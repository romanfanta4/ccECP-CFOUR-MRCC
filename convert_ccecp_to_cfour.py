#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
convert_ccecp_to_cfour.py -- canonical ccECP -> CFOUR/MRCC GENBAS + ECPDATA
=============================================================================

Convert the scalar-relativistic *canonical* ccECP pseudopotentials and their
associated Gaussian orbital basis sets from PseudopotentialLibrary.org into
CFOUR/MRCC library files:

    GENBAS   -- CFOUR Gaussian orbital basis-set library
    ECPDATA  -- CFOUR effective-core-potential library
    GENBAS.mrcc -- GENBAS with ECP entries appended for MRCC

Only ``recipes/<ELEMENT>/ccECP/`` is read.  Other recipe families
(``ccECP-soft``, ``eCEPP``, ``CEPP``, etc.) and spin-orbit/grid forms are
excluded deliberately.

The source logical names are retained. ECP labels remain exact (for example
``Fe:ccECP``), while orbital-basis labels receive ``-ccECP`` by default (for
example ``Fe:cc-pVQZ-ccECP``) to avoid collisions with installed basis
libraries. Use ``--basis-suffix ""`` only when exact source basis labels are
required and collisions have been checked. Optional ``*.upper`` companions
uppercase only the element symbol.

Requirements
------------
Python 3.8+ is sufficient for the converter itself.  ``basis_set_exchange`` is
optional but strongly recommended because it provides an independent basis and
ECP-convention cross-check.  Git is optional and is used only to record source
repository details in the generated README/summary.

Typical use
-----------
Clone the Pseudopotential Library, then run for example::

    git clone https://github.com/QMCPACK/pseudopotentiallibrary.git
    python3 -m pip install basis_set_exchange
    python3 convert_ccecp_to_cfour.py \
        --library ./pseudopotentiallibrary \
        --out ./ccECP_cfour_output

``--library`` may point either to the repository root or directly to its
``recipes`` directory.  When omitted, the script tries the environment variable
``PSEUDOPOTENTIAL_LIBRARY`` and a few locations next to the script/current
working directory.  Use ``--no-bse`` to run without the optional BSE checks.

Scientific conventions
----------------------
The converter preserves source decimal values and validates all analytic source
renderings that it knows how to parse.  ECPDATA is written as
``coefficient  N  alpha`` with radial factor ``r**(N-2)``; the local channel is
first, followed by the semi-local difference channels.  See the generated
README.md for the validation record and known program-compatibility notes.

"""

import argparse
import collections
import csv
import datetime
import json
import os
import re
import subprocess
import sys
from decimal import Decimal

VERSION = "1.2.0"

CANONICAL_FAMILY = "ccECP"

#: Suffix appended to orbital-basis labels so that ccECP-optimized basis sets
#: cannot be confused with installed all-electron or other basis sets carrying
#: the same conventional names.  A working-directory GENBAS can take precedence
#: over an installed library, so explicit disambiguation is safer for shared
#: inputs.  MRCC already distinguishes related sets with suffixes such as
#: `-PP`, `-DK`, and `-RI`.  Use `--basis-suffix ""` when exact source labels
#: are required and label collisions have been checked independently.
DEFAULT_BASIS_SUFFIX = "-ccECP"

#: The ECP label needs no suffix: `<El>:ccECP` does not collide with anything in
#: the CFOUR or MRCC libraries, and the family name must stay exactly `ccECP`.
ECP_LABEL_SUFFIX = ""


def basis_label_name(name, suffix):
    """Return the output orbital-basis label: source name plus optional suffix."""
    return name + suffix


#: Inserted before the label suffix for an auxiliary (density-fitting) set, so
#: `aug-cc-pCV5Z` becomes `aug-cc-pCV5Z-RI-ccECP`.  MRCC uses `-RI` for a
#: correlation-fitting set and `-RI-JK` for an SCF-fitting one; ORCA's AutoAux
#: /C sets are correlation-fitting, so `-RI` is the matching name.  MRCC already
#: ships entries such as `Ag:aug-cc-pV5Z-RI`, which is why the `-ccECP` suffix
#: still has to follow: `Ag:aug-cc-pCV5Z-RI-ccECP` collides with nothing.
RI_NAME_SUFFIX = "-RI"

#: file suffix of the library's auxiliary sets: `<El>.<basis>.AutoAuxC.orca`
AUX_FORMAT = "AutoAuxC.orca"


def ri_label_name(name, suffix):
    """Return the output auxiliary-basis label for orbital set `name`."""
    return name + RI_NAME_SUFFIX + suffix

FORBIDDEN_FAMILIES = (
    "ccECP-soft", "soft-ccECP-deprecated", "eCEPP", "CEPP", "BFD", "RRKJ", "TM",
)

AM_CHARS = "spdfghikl"

#: analytic representations this converter is able to read
BASIS_FORMATS = ("gamess", "nwchem", "gaussian")
ECP_FORMATS = ("plain", "gamess", "gaussian", "molpro", "nwchem")
#: tie-break preference when several renderings carry the same precision
#: tie-break order when several renderings carry the *same* precision:
#: prefer the library's native text file, then the widely used GAMESS
#: rendering, and only fall back to MOLPRO/Gaussian last
ECP_PREFERENCE = ("plain", "gamess", "nwchem", "molpro", "gaussian")
BASIS_PREFERENCE = ("gamess", "nwchem", "gaussian")

ELEMENT_Z = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9,
    "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16,
    "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23,
    "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
    "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36, "Rb": 37,
    "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44,
    "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50, "Sb": 51,
    "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Ce": 58,
    "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65,
    "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71, "Hf": 72,
    "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79,
    "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86,
}


#: Source-data observations worth reporting to users.  These are not assumed
#: to be errors: the converter preserves the repository numerical data exactly and does
#: not attempt to reinterpret, reoptimize, or repair basis sets.
SOURCE_DATA_NOTES = [
    dict(element="Ir", kind="basis", name="aug-cc-pVQZ",
         check="source_data_note", severity="review",
         detail="an f exponent is 0.577106 in aug-cc-pVQZ while the corresponding "
                "cc-pVQZ value is 1.577106; the value also breaks the descending "
                "order of the f "
                "set (0.577106 < 0.656915).  Identical in the .gamess, .gaussian, "
                ".nwchem and .molpro renderings.  The source value is preserved "
                "exactly as supplied."),
    dict(element="Ir", kind="basis", name="aug-cc-pCVQZ",
         check="source_data_note", severity="review",
         detail="same 0.577106 vs 1.577106 f exponent as Ir/aug-cc-pVQZ."),
    dict(element="Na", kind="basis", name="aug-cc-pVDZ/TZ/QZ",
         check="source_data_note", severity="review",
         detail="like Mg: the contracted kernels match cc-pVXZ, but the outermost "
                "uncontracted exponents differ rather than forming a strict "
                "primitive superset (cc-pVDZ s 0.120857 vs aug-cc-pVDZ s "
                "0.065057 and 0.026027). The source values are preserved exactly."),
    dict(element="Mg", kind="basis", name="aug-cc-pVDZ/TZ/QZ/5Z",
         check="source_data_note", severity="review",
         detail="the contracted s/p kernels are identical to cc-pVXZ but every "
                "uncontracted exponent differs (cc-pVDZ s 0.029784 vs "
                "aug-cc-pVDZ s 0.125638 and 0.050283); the aug- sets are not "
                "strict primitive supersets of the cc- sets for Mg.  The source "
                "values are preserved exactly."),
    dict(element="Eu", kind="basis", name="aug-cc-pCVTZ",
         check="source_data_note", severity="review",
         detail="the g space is {5.415869, 1.498739, 0.599496} while cc-pCVTZ has "
                "{5.415869, 1.498739, 0.563469}: the augmented set is not a strict "
                "primitive superset in the g space.  The source values are "
                "preserved exactly."),
]


#: Target-program compatibility notes.  These do not indicate a problem with
#: the converted numerical data, and the converter does not modify the ECP to
#: work around program-specific behavior.
KNOWN_PROGRAM_LIMITATIONS = [
    dict(element="H;He", kind="program", name="ccECP",
         check="mrcc_ncore_zero_unsupported", severity="program",
         detail="Executable validation with MRCC 25.1.1 and 26.1.1 shows that "
                "ECP entries with NCORE = 0 are not handled correctly.  This "
                "affects the canonical ccECPs of H and He, which keep all "
                "electrons explicit.  The converter writes the correct "
                "NCORE = 0 ECPDATA record without a numerical workaround.  "
                "Validate H and He independently before production use with "
                "these MRCC releases."),
]


#: MRCC 25.1.1 and 26.1.1 read the ECPDATA NCORE/LMAX record with this
#: fixed-column Fortran format:
#:
#:     read(gbasfile,"(11x,i3,11x,i1)") ncorecp(iatoms),lmax
#:
#: so NCORE must sit in columns 12-14 (I3) and LMAX in column 26 (I1).
MRCC_NCORE_COLUMNS = (12, 14)          # 1-based, inclusive
MRCC_LMAX_COLUMN = 26                  # 1-based


def format_ncore_lmax_line(ncore, lmax):
    """Format the CFOUR/MRCC ECPDATA NCORE/LMAX record.

    MRCC 25.1.1 and 26.1.1 read this record with

        read(gbasfile, "(11x,i3,11x,i1)") ncorecp, lmax

    so NCORE must occupy columns 12-14 and LMAX column 26.  This fixed-width
    layout is also valid CFOUR input.  The explicit checks below make the
    column contract part of the converter rather than relying on visual spacing.
    """
    if not (0 <= ncore <= 999):
        raise ConversionError("NCORE=%r does not fit the I3 field MRCC reads "
                              "from columns %d-%d" % ((ncore,) + MRCC_NCORE_COLUMNS))
    if not (0 <= lmax <= 9):
        raise ConversionError("LMAX=%r does not fit the I1 field MRCC reads "
                              "from column %d" % (lmax, MRCC_LMAX_COLUMN))
    line = "    NCORE =%3d    LMAX = %1d" % (ncore, lmax)
    # the whole point of this function is the column layout, so check it here
    lo, hi = MRCC_NCORE_COLUMNS
    assert int(line[lo - 1:hi]) == ncore, \
        "NCORE not in columns %d-%d of %r" % (lo, hi, line)
    assert int(line[MRCC_LMAX_COLUMN - 1]) == lmax, \
        "LMAX not in column %d of %r" % (MRCC_LMAX_COLUMN, line)
    return line


def mrcc_read_ncore_lmax(line):
    """Read a NCORE/LMAX record exactly as MRCC <= 26.1.1 does.

    Emulates `read(gbasfile,"(11x,i3,11x,i1)") ncorecp,lmax`, including the
    Fortran rule that an all-blank integer field reads as zero.  Used by the
    regression test for the fixed-column record.
    """
    def fortran_int(field):
        field = field.strip()
        return int(field) if field else 0

    lo, hi = MRCC_NCORE_COLUMNS
    return (fortran_int(line[lo - 1:hi]),
            fortran_int(line[MRCC_LMAX_COLUMN - 1:MRCC_LMAX_COLUMN]))


#: (element, NCORE, LMAX) covering one- and two-digit NCORE and LMAX 1..4
NCORE_LMAX_TEST_CASES = (("H", 0, 1), ("C", 2, 1), ("Fe", 10, 2),
                         ("Ru", 28, 3), ("Pt", 60, 4))


def selftest_ncore_lmax():
    """Regression-test the fixed-column NCORE/LMAX record.

    The test checks exact character positions and round-trips each record
    through an emulation of MRCC's fixed-format reader.
    """
    expected_ncore_field = {0: "  0", 2: "  2", 10: " 10", 28: " 28", 60: " 60"}
    results = []
    for elem, ncore, lmax in NCORE_LMAX_TEST_CASES:
        line = format_ncore_lmax_line(ncore, lmax)
        assert line[11:14] == expected_ncore_field[ncore], \
            "%s: columns 12-14 are %r, expected %r" % (
                elem, line[11:14], expected_ncore_field[ncore])
        assert line[25] == str(lmax), \
            "%s: column 26 is %r, expected %r" % (elem, line[25], str(lmax))
        assert mrcc_read_ncore_lmax(line) == (ncore, lmax), \
            "%s: MRCC-emulated read gave %r" % (elem, mrcc_read_ncore_lmax(line))
        results.append((elem, ncore, lmax, line, True))
    return results


class ConversionError(Exception):
    """Raised when a source entry cannot be converted unambiguously."""


# --------------------------------------------------------------------------- #
#  numeric helpers -- everything is kept as the *source decimal string* so that
#  no floating-point round-trip ever occurs; Decimal is used only to compare.
# --------------------------------------------------------------------------- #

def dec(s):
    return Decimal(str(s).replace("D", "E").replace("d", "e"))


def plain_decimal(s):
    """Return the value written in plain fixed-point form.

    The MOLPRO renderings use Fortran exponent notation
    (`2.32209171361153e+01`).  `Decimal` is an exact decimal type, so rewriting
    such a token as `23.2209171361153` changes no digit; it only avoids relying
    on how a given Fortran runtime parses `e+01` in list-directed input.  The
    conversion is asserted to be value-preserving.
    """
    d = dec(s)
    out = format(d, "f")
    assert dec(out) == d, "fixed-point rendering changed the value"
    return out


def norm_sig(s):
    """Significant digits of a written value, ignoring zero padding."""
    d = dec(s)
    if d == 0:
        return 0
    return len(d.normalize().as_tuple().digits)


#: every analytic rendering in the library carries at least this many
#: significant digits, so a value padded with trailing zeros is never credited
#: with less precision than this
MIN_WRITTEN_SIG = 6


def ulp(s):
    """Unit in the last *significant* written place of a value.

    Trailing zeros do not carry information -- `11.116996000000` in the
    `.gaussian` renderings really states six decimals -- so the unit in the
    last place is derived from the value with trailing zeros stripped, floored
    at MIN_WRITTEN_SIG significant digits.
    """
    d = dec(s)
    if d == 0:
        return Decimal(10) ** (1 - MIN_WRITTEN_SIG)
    t = d.normalize().as_tuple()
    nsig = max(len(t.digits), MIN_WRITTEN_SIG)
    magnitude = len(t.digits) + t.exponent            # floor(log10|d|) + 1
    return Decimal(10) ** (magnitude - nsig)


def consistent(a, b):
    """True when `a` and `b` are two renderings of the *same* number.

    Two decimal strings are consistent when they differ by no more than one
    unit in the last significant place of the coarser rendering, i.e. when the
    more precise value maps onto the coarser one under rounding *or*
    truncation.  (Truncation does occur in the library: Co's local channel is
    `25.00124115981202` in the MOLPRO rendering and `25.00124115` -- not
    ...16 -- in the GAMESS/NWChem/native renderings.)  Returns
    (consistent?, absolute difference, difference in units of that last-place
    unit).  The third value is a ULP fraction, not a relative error.
    """
    da, db = dec(a), dec(b)
    if da == db:
        return True, Decimal(0), Decimal(0)
    diff = abs(da - db)
    u = max(ulp(a), ulp(b))
    tol = u * Decimal("1.000000001")
    frac = diff / u
    return diff <= tol, diff, frac


def am_from_char(ch):
    ch = ch.lower()
    if ch not in AM_CHARS:
        raise ConversionError("unknown angular-momentum label %r" % ch)
    return AM_CHARS.index(ch)


def char_from_am(am):
    if am >= len(AM_CHARS):
        raise ConversionError("angular momentum %d out of table" % am)
    return AM_CHARS[am]


# --------------------------------------------------------------------------- #
#  Basis-set parsers.  Each returns a list of "source shells", one entry per
#  *contracted function* exactly as written in the source file:
#      {"am": int, "exps": [str, ...], "coefs": [str, ...]}
# --------------------------------------------------------------------------- #

def parse_gamess_basis(text, element):
    """GAMESS `$DATA` shell blocks.  Two header spellings occur in the library:

        `S  13`         upper case, nprim only
        `s 9 1.00`      lower case, nprim + scale factor (B and C files)

    Primitive lines are always `index exponent coefficient`.
    """
    shells, cur = [], None
    for lineno, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        if not line or line[0] in "!#$":
            continue
        tok = line.split()
        if re.fullmatch(r"[A-Za-z]", tok[0]):
            if len(tok) < 2:
                raise ConversionError("line %d: malformed shell header %r" % (lineno, line))
            try:
                nprim = int(tok[1])
            except ValueError:
                raise ConversionError("line %d: malformed shell header %r" % (lineno, line))
            if len(tok) >= 3 and dec(tok[2]) != Decimal(1):
                raise ConversionError("line %d: unsupported scale factor %r" % (lineno, tok[2]))
            if len(tok) > 3:
                raise ConversionError("line %d: unexpected tokens %r" % (lineno, line))
            if cur is not None:
                _close_shell(cur, shells)
            cur = {"am": am_from_char(tok[0]), "nprim": nprim, "exps": [], "coefs": []}
        elif re.fullmatch(r"[A-Za-z]{2,}", tok[0]):
            raise ConversionError("line %d: composite shell label %r not supported"
                                  % (lineno, tok[0]))
        else:
            if cur is None:
                raise ConversionError("line %d: primitive before any shell header" % lineno)
            if len(tok) != 3:
                raise ConversionError("line %d: expected 'idx exp coef', got %r" % (lineno, line))
            if int(tok[0]) != len(cur["exps"]) + 1:
                raise ConversionError("line %d: primitive index out of sequence" % lineno)
            cur["exps"].append(tok[1])
            cur["coefs"].append(tok[2])
    if cur is not None:
        _close_shell(cur, shells)
    if not shells:
        raise ConversionError("no shells found")
    return shells


def _close_shell(cur, shells):
    if len(cur["exps"]) != cur["nprim"]:
        raise ConversionError("shell declares %d primitives, %d given"
                              % (cur["nprim"], len(cur["exps"])))
    shells.append({"am": cur["am"], "exps": cur["exps"], "coefs": cur["coefs"]})


def parse_autoauxc(text, element):
    """ORCA `NewAuxCGTO` auxiliary basis, as shipped in `*.AutoAuxC.orca`.

    The shell body is the GAMESS spelling already handled by
    `parse_gamess_basis`; only the wrapper differs, being two lines at the top

        # Auxiliary/C basis set generated by ORCA AutoAux for Ag
        NewAuxCGTO Ag

    and a terminator at the bottom

        end;

    so the wrapper is removed and the body handed to the GAMESS parser.  No
    numerical field is touched here.

    The header element is checked against the file's element rather than
    trusted, since an auxiliary set carries no other identifying information.
    """
    lines = text.splitlines()
    body, seen_header, terminated = [], False, False
    for lineno, raw in enumerate(lines, 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if not seen_header:
            m = re.fullmatch(r"NewAuxC?GTO\s+([A-Za-z]{1,2})", line)
            if m is None:
                raise ConversionError(
                    "line %d: expected 'NewAuxCGTO <element>', got %r"
                    % (lineno, line))
            if m.group(1).capitalize() != element:
                raise ConversionError(
                    "line %d: header names element %r but the file is %r"
                    % (lineno, m.group(1), element))
            seen_header = True
            continue
        if line.rstrip(";").lower() == "end":
            terminated = True
            continue
        if terminated:
            raise ConversionError("line %d: content after the terminating 'end'"
                                  % lineno)
        body.append(raw)
    if not seen_header:
        raise ConversionError("no 'NewAuxCGTO <element>' header found")
    if not terminated:
        raise ConversionError("no terminating 'end' found")
    return parse_gamess_basis("\n".join(body) + "\n", element)


def aux_is_uncontracted(shells):
    """Is every shell a single primitive with coefficient exactly 1?

    That is what ORCA AutoAux emits, and what makes the CFOUR coefficient
    matrix an identity, exactly as MRCC's own `*-RI` entries are stored.  A set
    that is contracted instead still converts correctly; it is only reported,
    never altered.
    """
    for sh in shells:
        if len(sh["exps"]) != 1 or dec(sh["coefs"][0]) != Decimal(1):
            return False
    return True


def parse_nwchem_basis(text, element):
    """NWChem `<Elem> <L>` blocks (general-contraction columns supported)."""
    shells, exps, cols, am = [], [], [], None

    def flush(am, exps, cols):
        for j in range(len(cols)):
            shells.append({"am": am, "exps": list(exps), "coefs": list(cols[j])})

    for lineno, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        if not line or line[0] == "#":
            continue
        low = line.lower()
        if low in ("basis", "end", "spherical", "cartesian") or low.startswith("basis "):
            continue
        tok = line.split()
        if (len(tok) == 2 and re.fullmatch(r"[A-Za-z]{1,2}", tok[0])
                and re.fullmatch(r"[A-Za-z]", tok[1])):
            if tok[0].lower() != element.lower():
                raise ConversionError("line %d: header for %r in %s file"
                                      % (lineno, tok[0], element))
            if am is not None:
                flush(am, exps, cols)
            am, exps, cols = am_from_char(tok[1]), [], []
        else:
            if am is None:
                raise ConversionError("line %d: primitive before any shell header" % lineno)
            if len(tok) < 2:
                raise ConversionError("line %d: malformed primitive %r" % (lineno, line))
            if not cols:
                cols = [[] for _ in range(len(tok) - 1)]
            if len(tok) - 1 != len(cols):
                raise ConversionError("line %d: ragged contraction block" % lineno)
            exps.append(tok[0])
            for j, c in enumerate(tok[1:]):
                cols[j].append(c)
    if am is not None:
        flush(am, exps, cols)
    if not shells:
        raise ConversionError("no shells found")
    return shells


def parse_gaussian_basis(text, element):
    """Gaussian-94 style: `Fe 0` / `S 13 1.00` / rows / `****`."""
    shells = []
    lines = [l.strip() for l in text.splitlines()]
    i = 0
    while i < len(lines) and (not lines[i] or lines[i][0] == "!"):
        i += 1
    if i < len(lines):
        tok = lines[i].split()
        if len(tok) == 2 and tok[1] == "0":
            i += 1
    while i < len(lines):
        line = lines[i]
        i += 1
        if not line or line[0] == "!" or line.startswith("****"):
            continue
        tok = line.split()
        if not re.fullmatch(r"[A-Za-z]", tok[0]):
            raise ConversionError("unexpected line %r in gaussian basis" % line)
        am, nprim = am_from_char(tok[0]), int(tok[1])
        exps, coefs = [], []
        for _ in range(nprim):
            while i < len(lines) and not lines[i]:
                i += 1
            p = lines[i].split()
            i += 1
            if len(p) != 2:
                raise ConversionError("malformed gaussian primitive %r" % p)
            exps.append(p[0])
            coefs.append(p[1])
        shells.append({"am": am, "exps": exps, "coefs": coefs})
    if not shells:
        raise ConversionError("no shells found")
    return shells


BASIS_PARSERS = {
    "gamess": parse_gamess_basis,
    "nwchem": parse_nwchem_basis,
    "gaussian": parse_gaussian_basis,
}


# --------------------------------------------------------------------------- #
#  ECP parsers.  Each returns
#      {"ncore": int|None, "zeff": int|None, "lmax": int,
#       "channels": [(label, [(coef, N, alpha), ...]), ...]}
#  with channels[0] the local channel and channels[1:] l = 0 .. lmax-1.
#  Spin-orbit blocks are never returned.
# --------------------------------------------------------------------------- #

def chan_labels(lmax):
    loc = char_from_am(lmax)
    return [loc] + ["%s-%s" % (char_from_am(l), loc) for l in range(lmax)]


def parse_gamess_ecp(text, element):
    """`Fe-ccECP GEN <ncore> <lmax>`, then LMAX+1 blocks of
    (nterms, nterms x `coef N alpha`); local block first.  Carries no SO data."""
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    head = lines[0].split()
    if len(head) < 4 or head[1].upper() != "GEN":
        raise ConversionError("unrecognised gamess ECP header %r" % lines[0])
    ncore, lmax = int(head[2]), int(head[3])
    blocks, i = [], 1
    while i < len(lines):
        nterms = int(lines[i]); i += 1
        terms = []
        for _ in range(nterms):
            t = lines[i].split(); i += 1
            if len(t) != 3:
                raise ConversionError("malformed gamess ECP term %r" % t)
            terms.append((t[0], int(t[1]), t[2]))
        blocks.append(terms)
    if len(blocks) != lmax + 1:
        raise ConversionError("gamess ECP: LMAX=%d needs %d blocks, %d found"
                              % (lmax, lmax + 1, len(blocks)))
    return {"ncore": ncore, "zeff": None, "lmax": lmax,
            "channels": list(zip(chan_labels(lmax), blocks)), "n_so_discarded": 0}


def parse_gaussian_ecp(text, element):
    """`Fe 0` / `QMC <lmax> <ncore>` / (comment, nterms, nterms x `N alpha coef`)."""
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    i = 0
    if len(lines[i].split()) == 2 and lines[i].split()[1] == "0":
        i += 1
    head = lines[i].split(); i += 1
    if len(head) != 3:
        raise ConversionError("unrecognised gaussian ECP header %r" % lines[i - 1])
    lmax, ncore = int(head[1]), int(head[2])
    blocks = []
    while i < len(lines):
        if not re.fullmatch(r"\d+", lines[i]):
            i += 1
            continue
        nterms = int(lines[i]); i += 1
        terms = []
        for _ in range(nterms):
            t = lines[i].split(); i += 1
            if len(t) != 3:
                raise ConversionError("malformed gaussian ECP term %r" % t)
            terms.append((t[2], int(t[0]), t[1]))
        blocks.append(terms)
    if len(blocks) != lmax + 1:
        raise ConversionError("gaussian ECP: LMAX=%d needs %d blocks, %d found"
                              % (lmax, lmax + 1, len(blocks)))
    return {"ncore": ncore, "zeff": None, "lmax": lmax,
            "channels": list(zip(chan_labels(lmax), blocks)), "n_so_discarded": 0}


def parse_nwchem_ecp(text, element):
    """NWChem ECP; only the scalar `ecp ... end` section is read, a following
    `so ... end` (spin-orbit) section is discarded."""
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    n_so = 0
    for k, l in enumerate(lines):
        if l.lower() == "so":
            n_so = sum(1 for x in lines[k:]
                       if len(x.split()) == 2
                       and re.fullmatch(r"[A-Za-z]{1,2}", x.split()[0])
                       and re.fullmatch(r"[A-Za-z]", x.split()[1]))
            lines = lines[:k]
            break
    ncore, order, blocks, cur = None, [], [], None
    for raw in lines:
        if raw.lower() in ("ecp", "end"):
            continue
        tok = raw.split()
        if len(tok) == 3 and tok[1].lower() == "nelec":
            ncore = int(tok[2])
            continue
        if (len(tok) == 2 and re.fullmatch(r"[A-Za-z]{1,2}", tok[0])
                and re.fullmatch(r"[A-Za-z]{1,2}", tok[1])):
            cur = []
            order.append(tok[1].lower())
            blocks.append(cur)
            continue
        if cur is None:
            raise ConversionError("nwchem ECP term before any channel header: %r" % raw)
        if len(tok) != 3:
            raise ConversionError("malformed nwchem ECP term %r" % raw)
        cur.append((tok[2], int(tok[0]), tok[1]))
    if ncore is None:
        raise ConversionError("nwchem ECP: no `nelec` line")
    if not order or order[0] != "ul":
        raise ConversionError("nwchem ECP: first channel %r, expected `ul`"
                              % (order[0] if order else None))
    lmax = len(order) - 1
    for k, lab in enumerate(order[1:]):
        if am_from_char(lab) != k:
            raise ConversionError("nwchem ECP channel order anomaly: %r" % order)
    return {"ncore": ncore, "zeff": None, "lmax": lmax,
            "channels": list(zip(chan_labels(lmax), blocks)), "n_so_discarded": n_so}


def parse_molpro_ecp(text, element):
    """MOLPRO `ECP,<El>,<ncore>,<lmax>,<nso>`; the first LMAX+1 blocks are the
    scalar potential, the remaining `nso` blocks are spin-orbit and dropped."""
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    head = [x for x in re.split(r"[,\s;]+", lines[0]) if x]
    if head[0].lower() != "ecp":
        raise ConversionError("unrecognised molpro ECP header %r" % lines[0])
    ncore, lmax = int(head[2]), int(head[3])
    blocks, i = [], 1
    while i < len(lines):
        first = lines[i].split("!")[0].strip().rstrip(",;")
        nterms = int(first.split(",")[0]); i += 1
        terms = []
        for _ in range(nterms):
            t = [x.strip() for x in
                 lines[i].split("!")[0].strip().rstrip(",;").split(",")]
            i += 1
            if len(t) != 3:
                raise ConversionError("malformed molpro ECP term %r" % t)
            terms.append((t[2], int(t[0]), t[1]))
        blocks.append(terms)
    if len(blocks) < lmax + 1:
        raise ConversionError("molpro ECP: only %d of %d scalar blocks present"
                              % (len(blocks), lmax + 1))
    return {"ncore": ncore, "zeff": None, "lmax": lmax,
            "channels": list(zip(chan_labels(lmax), blocks[:lmax + 1])),
            "n_so_discarded": len(blocks) - (lmax + 1)}


def parse_plain_ccecp(text, element):
    """The library's native `<El>.ccECP` text form:

        Zeff  nchannels
        n_0 n_1 ... n_{lmax-1} n_local
        <blocks in the order l = 0, 1, ..., lmax-1, then the local channel>
    """
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    h = lines[0].split()
    zeff, nchan = int(dec(h[0])), int(h[1])
    counts = [int(x) for x in lines[1].split()]
    if len(counts) != nchan:
        raise ConversionError("plain ccECP: %d channels declared, %d counts"
                              % (nchan, len(counts)))
    lmax = nchan - 1
    blocks, i = [], 2
    for n in counts:
        terms = []
        for _ in range(n):
            t = lines[i].split(); i += 1
            if len(t) != 3:
                raise ConversionError("malformed plain ccECP term %r" % t)
            terms.append((t[2], int(t[0]), t[1]))
        blocks.append(terms)
    if i != len(lines):
        raise ConversionError("plain ccECP: %d trailing lines" % (len(lines) - i))
    ordered = [blocks[-1]] + blocks[:-1]        # local channel first
    return {"ncore": ELEMENT_Z[element] - zeff, "zeff": zeff, "lmax": lmax,
            "channels": list(zip(chan_labels(lmax), ordered)), "n_so_discarded": 0}


ECP_PARSERS = {
    "plain": parse_plain_ccecp,
    "gamess": parse_gamess_ecp,
    "gaussian": parse_gaussian_ecp,
    "nwchem": parse_nwchem_ecp,
    "molpro": parse_molpro_ecp,
}


# --------------------------------------------------------------------------- #
#  Comparison of independently parsed representations
# --------------------------------------------------------------------------- #

def compare_shells(a, b):
    """Compare two source-shell lists.  Returns (ok, identical, max_rel, msg)."""
    if len(a) != len(b):
        return False, False, None, "shell count %d vs %d" % (len(a), len(b))
    identical, worst = True, Decimal(0)
    for k, (x, y) in enumerate(zip(a, b)):
        if x["am"] != y["am"]:
            return False, False, None, "shell %d: am %d vs %d" % (k, x["am"], y["am"])
        if len(x["exps"]) != len(y["exps"]):
            return False, False, None, ("shell %d: nprim %d vs %d"
                                        % (k, len(x["exps"]), len(y["exps"])))
        for i in range(len(x["exps"])):
            for key in ("exps", "coefs"):
                ok, _, rel = consistent(x[key][i], y[key][i])
                if not ok:
                    return False, False, rel, ("shell %d prim %d %s: %s vs %s"
                                               % (k, i, key, x[key][i], y[key][i]))
                if rel:
                    identical = False
                    worst = max(worst, rel)
    return True, identical, worst, ""


def compare_ecp(a, b):
    """Compare two ECP dicts channel-by-channel (terms as sorted multisets).

    Returns (ok, identical, max_rel, msg).
    """
    if a["lmax"] != b["lmax"]:
        return False, False, None, "LMAX %d vs %d" % (a["lmax"], b["lmax"])
    if len(a["channels"]) != len(b["channels"]):
        return False, False, None, "channel count %d vs %d" % (
            len(a["channels"]), len(b["channels"]))
    if a["ncore"] is not None and b["ncore"] is not None and a["ncore"] != b["ncore"]:
        return False, False, None, "NCORE %d vs %d" % (a["ncore"], b["ncore"])
    identical, worst = True, Decimal(0)
    key = lambda t: (t[1], dec(t[2]), dec(t[0]))
    for (la, ta), (lb, tb) in zip(a["channels"], b["channels"]):
        if la != lb:
            return False, False, None, "channel label %s vs %s" % (la, lb)
        if len(ta) != len(tb):
            return False, False, None, ("channel %s: %d vs %d terms"
                                        % (la, len(ta), len(tb)))
        for x, y in zip(sorted(ta, key=key), sorted(tb, key=key)):
            if x[1] != y[1]:
                return False, False, None, ("channel %s: radial power %d vs %d"
                                            % (la, x[1], y[1]))
            for i in (0, 2):
                ok, _, rel = consistent(x[i], y[i])
                if not ok:
                    return False, False, rel, ("channel %s: %s vs %s"
                                               % (la, x[i], y[i]))
                if rel:
                    identical = False
                    worst = max(worst, rel)
    return True, identical, worst, ""


def ecp_precision(ecp):
    return sum(norm_sig(t[i]) for _, terms in ecp["channels"] for t in terms
               for i in (0, 2))


def basis_precision(shells):
    return sum(norm_sig(v) for sh in shells
               for key in ("exps", "coefs") for v in sh[key])


# --------------------------------------------------------------------------- #
#  GENBAS construction
# --------------------------------------------------------------------------- #

def build_genbas_blocks(shells):
    """Group source shells by angular momentum into CFOUR's layout.

    For each L, the primitive exponents are collected without numerical
    modification and then sorted in descending order to form the single
    exponent list CFOUR requires. Every source contracted function becomes one
    column of the coefficient matrix. Entries
    a given source shell does not use are `None` and are written as exact zero.
    Nothing is renormalised, recontracted or dropped; two primitives are
    identified only if their source decimals are numerically *exactly* equal.
    """
    byam = collections.OrderedDict()
    for sh in shells:
        byam.setdefault(sh["am"], []).append(sh)

    blocks, notes = [], []
    for am in sorted(byam):
        group = byam[am]
        order, index = [], {}
        for sh in group:
            here = set()
            for e in sh["exps"]:
                d = dec(e)
                if d in here:
                    raise ConversionError("L=%d: exponent %s repeated inside one "
                                          "contracted function" % (am, e))
                here.add(d)
                if d not in index:
                    index[d] = len(order)
                    order.append((d, e))
        ne, nc = len(order), len(group)
        mat = [[None] * nc for _ in range(ne)]
        for j, sh in enumerate(group):
            for e, c in zip(sh["exps"], sh["coefs"]):
                mat[index[dec(e)]][j] = c
        expd = [d for d, _ in order]
        # The union built above follows the order in which the source file
        # happens to list its shells, which for the core-valence sets puts tight
        # primitives after diffuse ones.  Reorder to descending exponent: this
        # permutes terms inside a sum and therefore changes no numerical value,
        # it only produces the canonical layout every code writes.
        was_sorted = all(expd[k] > expd[k + 1] for k in range(ne - 1))
        perm = sorted(range(ne), key=lambda i: -expd[i])
        order = [order[i] for i in perm]
        expd = [expd[i] for i in perm]
        mat = [mat[i] for i in perm]
        if not was_sorted:
            notes.append("L=%d: source listed primitives out of order, "
                         "re-sorted descending (value-preserving permutation)" % am)
        for k in range(ne - 1):
            if expd[k] != expd[k + 1] and abs(expd[k] - expd[k + 1]) \
                    < Decimal("1e-10") * max(abs(expd[k]), abs(expd[k + 1])):
                notes.append("L=%d NEAR-DEGENERATE primitives %s / %s"
                             % (am, order[k][1], order[k + 1][1]))
        blocks.append({"am": am, "exps": [t for _, t in order], "expd": expd,
                       "mat": mat, "resorted": not was_sorted,
                       "nprim_source": [len(t["exps"]) for t in group]})
    return blocks, notes


def format_genbas_entry(elem, name, blocks, comment, suffix=""):
    """Render one GENBAS entry in the CFOUR *new* format.

    Per the CFOUR specification the numeric records of an entry are read
    list-directed (Fortran `*`). Values are carried with exact Decimal arithmetic;
    scientific-notation tokens may be rewritten in fixed-point form by
    `plain_decimal()`, which is asserted to preserve the numerical value exactly.
    """
    blocks = [dict(b, exps=[plain_decimal(e) for e in b["exps"]],
                   mat=[[None if v is None else plain_decimal(v) for v in row]
                        for row in b["mat"]]) for b in blocks]
    exps = [e for b in blocks for e in b["exps"]]
    coefs = [c for b in blocks for row in b["mat"] for c in row if c is not None]
    ew = max(len(s) for s in exps) + 3
    cw = max([len(s) for s in coefs] + [10]) + 3
    per_e = max(1, min(5, 78 // ew))

    out = ["%s:%s" % (elem, basis_label_name(name, suffix)), comment, "",
           "%3d" % len(blocks),
           "".join("%5d" % b["am"] for b in blocks),
           "".join("%5d" % len(b["mat"][0]) for b in blocks),
           "".join("%5d" % len(b["exps"]) for b in blocks),
           ""]
    for b in blocks:
        for k in range(0, len(b["exps"]), per_e):
            out.append("".join("%*s" % (ew, s) for s in b["exps"][k:k + per_e]))
        out.append("")
        nc = len(b["mat"][0])
        per_c = max(1, min(nc, 78 // cw))
        for row in b["mat"]:
            vals = ["0.00000000" if v is None else v for v in row]
            for k in range(0, nc, per_c):
                out.append("".join("%*s" % (cw, v) for v in vals[k:k + per_c]))
        out.append("")
    return "\n".join(out) + "\n"


def parse_genbas_entry(text):
    """Independent reader for one GENBAS entry (round-trip verification)."""
    lines = text.splitlines()
    label = lines[0].strip()
    toks = []
    for l in lines[3:]:
        toks.extend(l.split())
    ns = int(toks[0]); p = 1
    am = [int(x) for x in toks[p:p + ns]]; p += ns
    nc = [int(x) for x in toks[p:p + ns]]; p += ns
    ne = [int(x) for x in toks[p:p + ns]]; p += ns
    blocks = []
    for k in range(ns):
        e = toks[p:p + ne[k]]; p += ne[k]
        flat = toks[p:p + ne[k] * nc[k]]; p += ne[k] * nc[k]
        blocks.append({"am": am[k], "exps": e,
                       "mat": [flat[r * nc[k]:(r + 1) * nc[k]] for r in range(ne[k])]})
    if p != len(toks):
        raise ConversionError("GENBAS round-trip: %d unconsumed tokens" % (len(toks) - p))
    return label, blocks


def compare_roundtrip(blocks, rblocks):
    if len(blocks) != len(rblocks):
        return False, "block count %d vs %d" % (len(blocks), len(rblocks))
    for b, r in zip(blocks, rblocks):
        if b["am"] != r["am"]:
            return False, "am %d vs %d" % (b["am"], r["am"])
        if len(b["exps"]) != len(r["exps"]):
            return False, "L=%d nprim %d vs %d" % (b["am"], len(b["exps"]), len(r["exps"]))
        for i, (x, y) in enumerate(zip(b["exps"], r["exps"])):
            if dec(x) != dec(y):
                return False, "L=%d exponent %d: %s vs %s" % (b["am"], i, x, y)
        for i in range(len(b["mat"])):
            if len(b["mat"][i]) != len(r["mat"][i]):
                return False, "L=%d row %d width" % (b["am"], i)
            for j in range(len(b["mat"][i])):
                want = Decimal(0) if b["mat"][i][j] is None else dec(b["mat"][i][j])
                if want != dec(r["mat"][i][j]):
                    return False, ("L=%d coefficient (%d,%d): %s vs %s"
                                   % (b["am"], i, j, b["mat"][i][j], r["mat"][i][j]))
    return True, ""


# --------------------------------------------------------------------------- #
#  ECPDATA rendering / re-reading
# --------------------------------------------------------------------------- #

def format_ecpdata_entry(elem, name, ecp, comment):
    chans = [(lab, [(plain_decimal(c), n, plain_decimal(a)) for c, n, a in ts])
             for lab, ts in ecp["channels"]]
    out = ["*", "%s:%s" % (elem, name), "# " + comment, "*",
           format_ncore_lmax_line(ecp["ncore"], ecp["lmax"])]
    w1 = max([len(t[0]) for _, ts in chans for t in ts] + [12]) + 2
    w2 = max([len(t[2]) for _, ts in chans for t in ts] + [12]) + 2
    for label, terms in chans:
        out.append(label)
        for coef, n, alpha in terms:
            out.append("%*s %4d %*s" % (w1, coef, n, w2, alpha))
    return "\n".join(out) + "\n"


def parse_ecpdata_entry(text):
    """Independent reader for one ECPDATA entry (round-trip verification)."""
    lines = [l for l in (x.rstrip() for x in text.splitlines()) if l.strip()]
    if lines[0].strip() != "*":
        raise ConversionError("ECPDATA entry does not start with '*'")
    label = lines[1].strip()
    if not lines[2].lstrip().startswith("#"):
        raise ConversionError("ECPDATA entry has no '#' comment line")
    if lines[3].strip() != "*":
        raise ConversionError("ECPDATA entry: missing second '*'")
    m = re.match(r"\s*NCORE\s*=\s*(\d+)\s+LMAX\s*=\s*(\d+)\s*$", lines[4])
    if not m:
        raise ConversionError("ECPDATA entry: bad NCORE/LMAX line %r" % lines[4])
    channels, cur = [], None
    for l in lines[5:]:
        t = l.split()
        if len(t) == 1 and re.fullmatch(r"[a-z](-[a-z])?", t[0]):
            cur = []
            channels.append((t[0], cur))
        else:
            if cur is None:
                raise ConversionError("ECPDATA entry: term before any channel %r" % l)
            if len(t) != 3:
                raise ConversionError("ECPDATA entry: bad term line %r" % l)
            cur.append((t[0], int(t[1]), t[2]))
    return label, {"ncore": int(m.group(1)), "zeff": None,
                   "lmax": int(m.group(2)), "channels": channels}


# --------------------------------------------------------------------------- #
#  Convention verification against CFOUR's own ECPDATA + BSE
# --------------------------------------------------------------------------- #

#: verbatim excerpts of the ECPDATA distributed with CFOUR
#: http://cfour.uni-mainz.de/cfour/uploads/Main.FormatOfECPDATAFile/ECPDATA
CFOUR_REFERENCE_ECPDATA = {
    "BR:ECP-28-MWB": """\
    NCORE = 28    LMAX = 3
f
   -8.16149293    2    2.7207
s-f
   61.51372099    2    5.0218
    9.02149299    2    2.5109
p-f
   53.87586402    2    4.2814
    4.62940227    2    2.1407
d-f
   20.84967744    2    2.8800
    2.96544431    2    1.4400
""",
    "BR:ECP-28 HAY & WADT": """\
  NCORE =    28           LMAX =     3
f
         -28.0000000    1         213.6143969
        -134.9268852    2          41.0585380
         -41.9271913    2           8.7086530
          -5.9336420    2           2.6074661
s-f
           3.0000000    0          54.1980682
          27.3430642    1          32.9053558
         118.8028847    2          13.6744890
          43.4354876    2           3.0341152
p-f
           5.0000000    0          54.2563340
          25.0504252    1          26.0095593
          92.6157463    2          28.2012995
          95.8249016    2           9.4341061
          26.2684983    2           2.5321764
d-f
           3.0000000    0          87.6328721
          22.5533557    1          61.7373377
         178.1241988    2          32.4385104
          76.9924162    2           8.7537199
           9.4818270    2           1.6633189
""",
}

#: (BSE basis name, element, maximum difference in units of the last
#: significant written place of the coarser representation).  compare_ecp()
#: reports this ULP-like quantity.  The Stuttgart/BSE rendering is rounded,
#: whereas the Hay-Wadt/LANL entry is expected to agree exactly.
CFOUR_REFERENCE_BSE = {
    "BR:ECP-28-MWB": ("stuttgart rlc ecp", "Br", "1.000000001"),
    "BR:ECP-28 HAY & WADT": ("lanl2dz", "Br", "0"),
}


def _parse_reference_block(text):
    lines = [l for l in text.splitlines() if l.strip()]
    m = re.match(r"\s*NCORE\s*=\s*(\d+)\s+LMAX\s*=\s*(\d+)", lines[0])
    channels, cur = [], None
    for l in lines[1:]:
        t = l.split()
        if len(t) == 1:
            cur = []
            channels.append((t[0], cur))
        else:
            cur.append((t[0], int(t[1]), t[2]))
    return {"ncore": int(m.group(1)), "zeff": None,
            "lmax": int(m.group(2)), "channels": channels}


def verify_conventions(log):
    results = []
    try:
        import basis_set_exchange as bse
    except ImportError:
        log("convention check: basis_set_exchange unavailable -- SKIPPED")
        return results
    for key, block in CFOUR_REFERENCE_ECPDATA.items():
        ref = _parse_reference_block(block)
        bs_name, elem, max_ulp = CFOUR_REFERENCE_BSE[key]
        try:
            data = bse.get_basis(bs_name, elements=[elem])["elements"][str(ELEMENT_Z[elem])]
        except Exception as exc:                                # pragma: no cover
            results.append((key, "ERROR", str(exc)))
            continue
        # some library files carry a dummy highest channel whose coefficients
        # are all exactly zero; drop it for the structural comparison only
        keep = [q for q in data["ecp_potentials"]
                if any(dec(c) != 0 for col in q["coefficients"] for c in col)]
        pots = sorted(keep, key=lambda x: x["angular_momentum"])
        pots = [pots[-1]] + pots[:-1]
        lmax = max(p["angular_momentum"][0] for p in keep)
        got = {"ncore": data["ecp_electrons"], "zeff": None, "lmax": lmax,
               "channels": []}
        for lab, p in zip(chan_labels(lmax), pots):
            got["channels"].append((lab, [(p["coefficients"][0][i],
                                           int(p["r_exponents"][i]),
                                           p["gaussian_exponents"][i])
                                          for i in range(len(p["r_exponents"]))]))
        ok, ident, worst, msg = compare_ecp(ref, got)
        good = ok and worst is not None and worst <= Decimal(max_ulp)
        results.append((key + "  vs BSE '" + bs_name + "'",
                        "PASS" if good else "FAIL",
                        "NCORE %s/%s LMAX %s/%s max_ulp_diff=%s %s"
                        % (ref["ncore"], got["ncore"], ref["lmax"], got["lmax"],
                           ("%.2e" % worst) if worst is not None else "n/a", msg)))
    for k, s, m in results:
        log("convention check [%s]: %s  (%s)" % (k, s, m))
    return results


def bse_read_shells(text, fmt, element):
    """Parse a source basis with BSE and return shells in our internal form."""
    import basis_set_exchange as bse
    bd = bse.readers.read.read_formatted_basis_str(text, basis_fmt=fmt,
                                                   validate=False,
                                                   as_component=False)
    els = bd["elements"]
    z = str(ELEMENT_Z[element])
    if z not in els:
        if len(els) != 1:
            raise ConversionError("BSE returned %d elements" % len(els))
        z = list(els)[0]
    shells = []
    for sh in els[z].get("electron_shells", []):
        for col in sh["coefficients"]:
            shells.append({"am": sh["angular_momentum"][0],
                           "exps": list(sh["exponents"]), "coefs": list(col)})
    return shells


# --------------------------------------------------------------------------- #
#  Library scanning
# --------------------------------------------------------------------------- #

ALL_SUFFIXES = ("gamess", "gaussian", "nwchem", "molpro", "dirac")


def scan_host_labels(path):
    """Collect `<El>:<name>` labels from an installed CFOUR/MRCC basis library.

    Accepts either a single GENBAS-style file or a directory of them, which is
    how MRCC 26 ships its library (`BASIS/<El>`).  Used only to count label
    collisions, so that the naming claims in the README are measured on the
    target system rather than asserted.
    """
    files = []
    if os.path.isdir(path):
        for name in sorted(os.listdir(path)):
            full = os.path.join(path, name)
            if os.path.isfile(full):
                files.append(full)
    elif os.path.isfile(path):
        files = [path]
    else:
        raise ConversionError("host basis library not found: %s" % path)
    labels = set()
    pat = re.compile(r"([A-Za-z]{1,2}):(\S+)\s*")
    for full in files:
        try:
            with open(full, errors="replace") as fh:
                for line in fh:
                    m = pat.fullmatch(line)
                    if m:
                        labels.add("%s:%s" % (m.group(1).capitalize(),
                                              m.group(2)))
        except OSError:
            continue
    return labels


def scan_library(recipes_dir):
    inv = collections.OrderedDict()
    excluded = collections.Counter()
    for elem in sorted(os.listdir(recipes_dir)):
        eldir = os.path.join(recipes_dir, elem)
        if not os.path.isdir(eldir):
            continue
        for fam in sorted(os.listdir(eldir)):
            if os.path.isdir(os.path.join(eldir, fam)) and fam != CANONICAL_FAMILY:
                excluded[fam] += 1
        ccdir = os.path.join(eldir, CANONICAL_FAMILY)
        if not os.path.isdir(ccdir):
            continue
        if elem not in ELEMENT_Z:
            raise ConversionError("unknown element directory %r" % elem)
        ecp_files, basis_files, skipped = {}, collections.defaultdict(dict), []
        aux_files = {}
        for f in sorted(os.listdir(ccdir)):
            if f in ("author.txt", "energies.txt"):
                continue
            if not f.startswith(elem + "."):
                skipped.append(f)
                continue
            rest = f[len(elem) + 1:]
            if ".SOREP." in "." + rest or ".AREP." in "." + rest:
                skipped.append(f)                    # grid-based AREP/SOREP data
                continue
            if rest.endswith("." + AUX_FORMAT):
                # `<El>.<basis>.AutoAuxC.orca`: the ORCA AutoAux /C auxiliary
                # set belonging to one orbital basis.  Keyed by that orbital
                # name so the two are matched up later.
                aux_files[rest[:-len("." + AUX_FORMAT)]] = f
                continue
            if rest == CANONICAL_FAMILY:
                ecp_files["plain"] = f
                continue
            m = re.fullmatch(r"(.+)\.([A-Za-z_0-9]+)", rest)
            if not m:
                skipped.append(f)
                continue
            logical, fmt = m.group(1), m.group(2)
            if fmt not in ALL_SUFFIXES:
                skipped.append(f)                    # xml / upf / rpt / *_deprecated
                continue
            if logical == CANONICAL_FAMILY:
                if fmt in ECP_PARSERS:
                    ecp_files[fmt] = f
                else:
                    skipped.append(f)                # .dirac
            else:
                basis_files[logical][fmt] = f
        author = ""
        apath = os.path.join(ccdir, "author.txt")
        if os.path.isfile(apath):
            with open(apath) as fh:
                author = fh.readline().strip()
        inv[elem] = {"element": elem, "Z": ELEMENT_Z[elem],
                     "recipe_dir": os.path.join("recipes", elem, CANONICAL_FAMILY),
                     "ecp_files": ecp_files,
                     "basis_files": {k: dict(v) for k, v in sorted(basis_files.items())},
                     "aux_files": dict(sorted(aux_files.items())),
                     "author": author, "skipped_files": skipped}
    return inv, excluded


def clean_reference(author_line):
    s = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r"\1 \2", author_line)
    return re.sub(r"\s+", " ", s).strip()


# --------------------------------------------------------------------------- #
#  Source selection
# --------------------------------------------------------------------------- #

def select_ecp(elem, ccdir, files, want, checks):
    """Parse every available analytic scalar ECP rendering, verify that they are
    mutually consistent, and return the highest-precision one."""
    parsed, errors = {}, {}
    for fmt in ECP_FORMATS:
        if fmt not in files:
            continue
        try:
            with open(os.path.join(ccdir, files[fmt])) as fh:
                parsed[fmt] = ECP_PARSERS[fmt](fh.read(), elem)
        except Exception as exc:
            errors[fmt] = str(exc)
    for fmt, msg in errors.items():
        checks.append(("parse_%s" % fmt, "FAIL", msg))
    if not parsed:
        raise ConversionError("no analytic ECP rendering could be parsed: %r" % errors)

    ranked = sorted(parsed, key=lambda f: (-ecp_precision(parsed[f]),
                                           ECP_PREFERENCE.index(f)))
    chosen = ranked[0] if want == "auto" else want
    if chosen not in parsed:
        raise ConversionError("requested ECP source .%s unavailable/unparsable" % chosen)

    ref = parsed[chosen]
    for fmt in ECP_FORMATS:
        if fmt == chosen or fmt not in parsed:
            continue
        ok, ident, worst, msg = compare_ecp(ref, parsed[fmt])
        checks.append(("crosscheck_%s" % fmt, "PASS" if ok else "FAIL",
                       msg if not ok else
                       ("identical" if ident else
                        "same numbers, %s written with fewer significant digits "
                        "(largest difference %.2f units in its last place, i.e. "
                        "%s)" % (fmt, worst,
                                 "rounded" if worst <= Decimal("0.5")
                                 else "truncated"))))
    best = ranked[0]
    if best == chosen:
        checks.append(("source_precision", "PASS",
                       ".%s is the most precise rendering available (%d "
                       "significant digits in total)"
                       % (chosen, ecp_precision(ref))))
    else:
        ok, ident, worst, msg = compare_ecp(parsed[best], ref)
        checks.append(("source_precision", "PASS",
                       ".%s carries %d significant digits, .%s carries %d; the "
                       "extra digits change no parameter by more than %.2f units "
                       "in the last place of .%s (see ECPDATA.maxprec)"
                       % (chosen, ecp_precision(ref), best,
                          ecp_precision(parsed[best]), worst or 0, chosen)))
    return chosen, ref, parsed, best


def select_basis(elem, ccdir, name, files, want, checks, bse_version):
    parsed, errors = {}, {}
    for fmt in BASIS_FORMATS:
        if fmt not in files:
            continue
        try:
            with open(os.path.join(ccdir, files[fmt])) as fh:
                parsed[fmt] = BASIS_PARSERS[fmt](fh.read(), elem)
        except Exception as exc:
            errors[fmt] = str(exc)
    for fmt, msg in errors.items():
        checks.append(("parse_%s" % fmt, "FAIL", msg))
    if not parsed:
        raise ConversionError("no analytic basis rendering could be parsed: %r" % errors)

    ranked = sorted(parsed, key=lambda f: (-basis_precision(parsed[f]),
                                           BASIS_PREFERENCE.index(f)))
    chosen = ranked[0] if want == "auto" else want
    if chosen not in parsed:
        raise ConversionError("requested basis source .%s unavailable/unparsable" % chosen)

    ref = parsed[chosen]
    for fmt in BASIS_FORMATS:
        if fmt == chosen or fmt not in parsed:
            continue
        ok, ident, worst, msg = compare_shells(ref, parsed[fmt])
        checks.append(("crosscheck_%s" % fmt, "PASS" if ok else "FAIL",
                       msg if not ok else
                       ("identical" if ident else
                        "same numbers, %s written with fewer significant digits "
                        "(largest difference %.2f units in its last place)"
                        % (fmt, worst))))
    # fully independent parse by the MolSSI Basis Set Exchange
    if bse_version and "gaussian" in files:
        try:
            with open(os.path.join(ccdir, files["gaussian"])) as fh:
                bsh = bse_read_shells(fh.read(), "gaussian94", elem)
            ok, ident, worst, msg = compare_shells(ref, bsh)
            checks.append(("crosscheck_bse", "PASS" if ok else "FAIL",
                           msg if not ok else
                           ("identical to the basis_set_exchange %s parse of the "
                            ".gaussian source" % bse_version if ident else
                            "agrees with the basis_set_exchange %s parse of the "
                            ".gaussian source to within its written precision "
                            "(%.2f ulp)" % (bse_version, worst))))
        except Exception as exc:
            checks.append(("crosscheck_bse", "FAIL", str(exc)))
    return chosen, ref, parsed


# --------------------------------------------------------------------------- #
#  Command-line path handling
# --------------------------------------------------------------------------- #

def _library_candidates():
    """Return plausible Pseudopotential Library locations in priority order.

    No candidate is trusted merely because of its name: ``resolve_library()``
    requires a real ``recipes`` directory before accepting it.  The environment
    variable is useful on clusters where the repository is installed centrally.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    cwd = os.path.abspath(os.getcwd())
    candidates = []
    env = os.environ.get("PSEUDOPOTENTIAL_LIBRARY")
    if env:
        candidates.append(env)
    candidates.extend([
        script_dir,                              # script copied into repository root
        os.path.join(script_dir, "pseudopotentiallibrary"),
        cwd,                                     # invoked from repository root
        os.path.join(cwd, "pseudopotentiallibrary"),
    ])
    # Preserve order while eliminating duplicate absolute paths.
    out = []
    for item in candidates:
        item = os.path.abspath(os.path.expanduser(item))
        if item not in out:
            out.append(item)
    return out


def resolve_library(user_path=None):
    """Resolve repository root and recipes directory.

    ``user_path`` may be either the Pseudopotential Library repository root or
    its ``recipes`` directory.  With no explicit path, a small deterministic
    set of local candidates is searched.
    """
    candidates = [user_path] if user_path else _library_candidates()
    tried = []
    for raw in candidates:
        if not raw:
            continue
        p = os.path.abspath(os.path.expanduser(raw))
        tried.append(p)
        if os.path.isdir(os.path.join(p, "recipes")):
            return p, os.path.join(p, "recipes")
        if os.path.basename(os.path.normpath(p)) == "recipes" and os.path.isdir(p):
            return os.path.dirname(p), p
    raise ConversionError(
        "could not find the Pseudopotential Library. Pass --library /path/to/"
        "pseudopotentiallibrary (or its recipes directory), or set "
        "PSEUDOPOTENTIAL_LIBRARY. Tried: %s" % ", ".join(tried))


def _path_is_inside(child, parent):
    """True if *child* lies inside *parent*; safe across Windows drive letters."""
    try:
        return os.path.commonpath([os.path.abspath(child), os.path.abspath(parent)]) \
            == os.path.abspath(parent)
    except ValueError:
        return False


# --------------------------------------------------------------------------- #
#  Main driver
# --------------------------------------------------------------------------- #

def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Convert canonical ccECP data to CFOUR/MRCC GENBAS + ECPDATA.",
        epilog=("Examples:\n"
                "  python3 %(prog)s --library ./pseudopotentiallibrary "
                "--out ./ccECP_cfour_output\n"
                "  PSEUDOPOTENTIAL_LIBRARY=/opt/pseudopotentiallibrary "
                "python3 %(prog)s --out ./output\n"
                "  python3 %(prog)s --library ./pseudopotentiallibrary/recipes "
                "--no-bse"),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--library", metavar="DIR", default=None,
        help=("Pseudopotential Library repository root or its recipes directory. "
              "If omitted, use $PSEUDOPOTENTIAL_LIBRARY or auto-detect a local "
              "pseudopotentiallibrary checkout."))
    ap.add_argument(
        "--out", metavar="DIR", default=None,
        help=("output directory. If omitted, use ./ccECP_cfour_output when the "
              "current directory is outside the source repository, otherwise use "
              "a sibling ccECP_cfour_output directory. Existing generated files "
              "with the same names are replaced."))
    ap.add_argument(
        "--basis-source", default="auto", choices=("auto",) + BASIS_FORMATS,
        help=("basis rendering to write; auto selects the most precise supported "
              "rendering after cross-checking all supported renderings (default: auto)"))
    ap.add_argument(
        "--ecp-source", default="gamess", choices=("auto",) + ECP_FORMATS,
        help=("ECP rendering to write after cross-checking all supported renderings "
              "(default: gamess; auto selects the most precise)"))
    ap.add_argument(
        "--basis-suffix", metavar="TEXT", default=DEFAULT_BASIS_SUFFIX,
        help=("suffix appended to orbital-basis labels to keep them distinct "
              "from the host program's all-electron sets of the same name "
              "(default: %(default)r; pass an empty string to use exact source "
              "basis names). The ECP label is never suffixed."))
    ap.add_argument(
        "--no-bse", action="store_true",
        help="skip optional MolSSI Basis Set Exchange cross-checks")
    ap.add_argument(
        "--no-by-element", action="store_true",
        help="do not write output/by_element/<El>/ convenience files")
    ap.add_argument(
        "--no-ri", action="store_true",
        help="skip the AutoAuxC auxiliary (RI) basis sets entirely")
    ap.add_argument(
        "--host-basis-library", metavar="PATH", default=None,
        help="installed CFOUR/MRCC basis library (a GENBAS-style file, or a "
             "directory of them such as MRCC's BASIS/). Only used to count "
             "label collisions for the naming report")
    ap.add_argument(
        "--extra-aux", metavar="DIR", default=None,
        help="directory of additional <El>.<basis>.AutoAuxC.orca files, used "
             "where the library ships none; each is recorded in "
             "ri_name_map.csv with source_origin=extra")
    ap.add_argument("--version", action="version", version="%(prog)s " + VERSION)
    args = ap.parse_args(argv)

    if not re.fullmatch(r"[A-Za-z0-9_.+\-]*", args.basis_suffix):
        ap.error("--basis-suffix may contain only letters, digits, '.', '_', '+', "
                 "and '-'; use an empty string for no suffix")

    try:
        library_root, recipes = resolve_library(args.library)
    except ConversionError as exc:
        ap.error(str(exc))
    args.library = library_root
    if args.out is None:
        cwd = os.path.abspath(os.getcwd())
        if _path_is_inside(cwd, library_root):
            outdir = os.path.join(os.path.dirname(library_root), "ccECP_cfour_output")
        else:
            outdir = os.path.join(cwd, "ccECP_cfour_output")
    else:
        outdir = os.path.abspath(os.path.expanduser(args.out))
    if _path_is_inside(outdir, library_root):
        ap.error("--out must be outside the Pseudopotential Library source tree; "
                 "the source repository is treated as read-only")

    # Validate the source tree before creating any output files.
    if not os.path.isdir(recipes):
        ap.error("recipes directory not found: %s" % recipes)
    try:
        inv_preview, excluded_preview = scan_library(recipes)
    except (OSError, ConversionError) as exc:
        ap.error("cannot read canonical ccECP recipes: %s" % exc)
    if not inv_preview:
        ap.error("no canonical recipes/<ELEMENT>/ccECP directories found under %s"
                 % recipes)

    os.makedirs(outdir, exist_ok=True)
    logfh = open(os.path.join(outdir, "conversion.log"), "w")

    def log(msg):
        print(msg)
        logfh.write(msg + "\n")
        logfh.flush()

    stamp = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    try:
        commit = subprocess.check_output(
            ["git", "-C", args.library, "log", "-1", "--format=%H"],
            text=True, stderr=subprocess.DEVNULL).strip()
        commit_date = subprocess.check_output(
            ["git", "-C", args.library, "log", "-1", "--format=%ci"],
            text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        commit = commit_date = "unknown"
    bse_version = None
    if not args.no_bse:
        try:
            import basis_set_exchange as bse
            bse_version = bse.version()
        except ImportError:
            bse_version = None

    log("convert_ccecp_to_cfour.py v%s" % VERSION)
    log("conversion date              : %s" % stamp)
    log("python                       : %s" % sys.version.split()[0])
    log("library root                 : %s" % args.library)
    log("recipes                      : %s" % recipes)
    log("output                       : %s" % outdir)
    log("library commit               : %s (%s)" % (commit, commit_date))
    log("basis_set_exchange           : %s" % (bse_version or "not used"))
    log("basis source selection       : %s" % args.basis_source)
    log("ECP   source selection       : %s" % args.ecp_source)
    log("")
    # The MRCC-compatible fixed-column NCORE/LMAX record is the one piece of
    # rigid formatting in the output, so prove it before writing anything.
    ncore_lmax_tests = selftest_ncore_lmax()
    log("MRCC fixed-column NCORE/LMAX record self-test:")
    for elem, ncore, lmax, line, _ in ncore_lmax_tests:
        log("    %-3s NCORE=%-3d LMAX=%d  cols 12-14=%r col 26=%r  -> MRCC "
            "reads %r  PASS"
            % (elem, ncore, lmax, line[11:14], line[25],
               mrcc_read_ncore_lmax(line)))
    log("")

    if args.no_bse:
        log("convention check: disabled by --no-bse -- SKIPPED")
        convention_results = []
    else:
        convention_results = verify_conventions(log)
    convention_failures = [r for r in convention_results if r[1] != "PASS"]
    log("")

    inv, excluded = inv_preview, excluded_preview
    log("canonical ccECP elements found : %d" % len(inv))
    log("recipe families present in the library but EXCLUDED by the ccECP rule:")
    for fam, n in sorted(excluded.items()):
        log("    %-28s (%d elements)" % (fam, n))
    log("")

    validation_rows, basis_map_rows, ecp_map_rows, naming_rows = [], [], [], []
    ri_entries, ri_map_rows, ri_file_of = [], [], {}
    anomalies = list(SOURCE_DATA_NOTES) + list(KNOWN_PROGRAM_LIMITATIONS)
    genbas_entries, ecpdata_entries, ecpdata_best_entries = [], [], []
    basis_file_of, ecp_file_of = {}, {}
    per_element, failures = {}, {"ecp": [], "basis": [], "ri": []}
    seen_labels = collections.Counter()
    src_used = {"ecp": collections.Counter(), "basis": collections.Counter(),
                "ecp_maxprec": collections.Counter()}

    # ------------------------------- ECPs ------------------------------ #
    for elem, meta in inv.items():
        z, ccdir = meta["Z"], os.path.join(recipes, elem, CANONICAL_FAMILY)
        ref_txt = clean_reference(meta["author"])
        checks = []
        try:
            chosen, ecp, parsed, best_fmt = select_ecp(
                elem, ccdir, meta["ecp_files"], args.ecp_source, checks)
        except Exception as exc:
            failures["ecp"].append((elem, str(exc)))
            validation_rows.append(dict(element=elem, kind="ecp",
                                        name=CANONICAL_FAMILY, check="source",
                                        status="FAIL", detail=str(exc)))
            log("ECP  %-3s FAILED: %s" % (elem, exc))
            continue

        ncore, lmax = ecp["ncore"], ecp["lmax"]
        nexpl = z - ncore
        loc = ecp["channels"][0][0]

        checks.append(("n_channels",
                       "PASS" if len(ecp["channels"]) == lmax + 1 else "FAIL",
                       "%d channels for LMAX=%d (local + l=0..%d)"
                       % (len(ecp["channels"]), lmax, lmax - 1)))

        # NCORE derived from the data, never from a filename: the r**-1 (N=1)
        # terms of the local channel must sum to Zeff = Z - NCORE
        ones = [t for t in ecp["channels"][0][1] if t[1] == 1]
        zsum = sum((dec(t[0]) for t in ones), Decimal(0))
        checks.append(("ncore_from_data",
                       "PASS" if zsum == Decimal(nexpl) else "FAIL",
                       "sum of local-channel N=1 coefficients = %s, Z-NCORE = %d "
                       "(%d term%s)" % (zsum, nexpl, len(ones),
                                        "" if len(ones) == 1 else "s")))
        if ecp["zeff"] is not None:
            checks.append(("zeff_header",
                           "PASS" if ecp["zeff"] == nexpl else "FAIL",
                           "native .ccECP header Zeff=%d" % ecp["zeff"]))
        # every parsed rendering must agree on NCORE
        ncset = {f: p["ncore"] for f, p in parsed.items() if p["ncore"] is not None}
        checks.append(("ncore_all_formats",
                       "PASS" if len(set(ncset.values())) == 1 else "FAIL",
                       "NCORE per rendering: %r" % ncset))

        fixed_record = format_ncore_lmax_line(ncore, lmax)
        fixed_ok = mrcc_read_ncore_lmax(fixed_record) == (ncore, lmax)
        checks.append(("mrcc_fixed_columns", "PASS" if fixed_ok else "FAIL",
                       "NCORE in columns 12-14=%r; LMAX in column 26=%r"
                       % (fixed_record[11:14], fixed_record[25])))

        rp = sorted({t[1] for _, ts in ecp["channels"] for t in ts})
        checks.append(("radial_powers", "PASS" if all(0 <= x <= 4 for x in rp) else "FAIL",
                       "N in %r, written through unchanged; radial factor r**(N-2)" % rp))
        checks.append(("no_r_minus_2_divergence",
                       "PASS" if 0 not in rp else "WARN",
                       "no N=0 (r**-2) terms" if 0 not in rp else "N=0 terms present"))

        so = sum(p.get("n_so_discarded", 0) for p in parsed.values())
        checks.append(("spin_orbit_excluded", "PASS",
                       "%d spin-orbit block(s) found in the source renderings and "
                       "discarded" % so))

        comment = ("canonical %s for %s: Z=%d NCORE=%d Zeff=%d LMAX=%d local=%s | "
                   "source %s | %s"
                   % (CANONICAL_FAMILY, elem, z, ncore, nexpl, lmax, loc,
                      os.path.join(meta["recipe_dir"], meta["ecp_files"][chosen]),
                      ref_txt))
        entry = format_ecpdata_entry(elem, CANONICAL_FAMILY, ecp, comment)
        try:
            rlabel, rt = parse_ecpdata_entry(entry)
            ok, ident, worst, msg = compare_ecp(ecp, rt)
            good = (ok and ident and rt["ncore"] == ncore and rt["lmax"] == lmax
                    and rlabel == "%s:%s" % (elem, CANONICAL_FAMILY))
            checks.append(("roundtrip", "PASS" if good else "FAIL",
                           msg or "written ECPDATA entry re-reads numerically exactly"))
        except Exception as exc:
            checks.append(("roundtrip", "FAIL", str(exc)))

        status = "PASS" if all(s.startswith("PASS") or s == "WARN"
                               for _, s, _ in checks) else "FAIL"
        for c, s, d in checks:
            validation_rows.append(dict(element=elem, kind="ecp",
                                        name=CANONICAL_FAMILY, check=c,
                                        status=s, detail=d))
        if status == "FAIL":
            failures["ecp"].append((elem, "; ".join(
                "%s=%s" % (c, s) for c, s, _ in checks
                if not (s.startswith("PASS") or s == "WARN"))))
            log("ECP  %-3s FAILED -- entry NOT written (%s)"
                % (elem, "; ".join("%s:%s" % (c, d) for c, s, d in checks
                                   if not (s.startswith("PASS") or s == "WARN"))))
            continue

        ecpdata_entries.append((elem, CANONICAL_FAMILY, entry))
        best_comment = ("canonical %s for %s: Z=%d NCORE=%d Zeff=%d LMAX=%d "
                        "local=%s | source %s | %s"
                        % (CANONICAL_FAMILY, elem, z, ncore, nexpl, lmax, loc,
                           os.path.join(meta["recipe_dir"],
                                        meta["ecp_files"][best_fmt]), ref_txt))
        ecpdata_best_entries.append(
            (elem, CANONICAL_FAMILY,
             format_ecpdata_entry(elem, CANONICAL_FAMILY, parsed[best_fmt],
                                  best_comment)))
        src_used["ecp_maxprec"][best_fmt] += 1
        label = "%s:%s" % (elem, CANONICAL_FAMILY)
        seen_labels[label] += 1
        src_used["ecp"][chosen] += 1

        ecp_map_rows.append(dict(
            element=elem, atomic_number=z,
            source_recipe_dir=meta["recipe_dir"],
            source_filename=meta["ecp_files"][chosen],
            source_format=chosen,
            all_source_formats=";".join(sorted(meta["ecp_files"])),
            library_family=CANONICAL_FAMILY, library_name=CANONICAL_FAMILY,
            cfour_name=CANONICAL_FAMILY, label_suffix=ECP_LABEL_SUFFIX,
            cfour_label=label,
            ncore=ncore, n_explicit_electrons=nexpl, lmax=lmax,
            local_channel=loc,
            nonlocal_channels=";".join(l for l, _ in ecp["channels"][1:]),
            n_terms_per_channel=";".join("%s:%d" % (l, len(t))
                                         for l, t in ecp["channels"]),
            radial_powers=";".join(str(x) for x in rp),
            precision_digits_per_format=";".join(
                "%s:%d" % (f, ecp_precision(p)) for f, p in sorted(parsed.items())),
            radial_power_convention="U_l(r) = sum_m c_m r**(N_m-2) exp(-alpha_m r**2)",
            local_channel_convention="local channel first, labelled l=LMAX; "
                                     "semi-local channels are U_l - U_LMAX",
            name_transformation="none: <Element>:<library name> verbatim",
            validation=status))
        naming_rows.append(dict(
            element=elem, source_recipe_dir=meta["recipe_dir"],
            source_filename=meta["ecp_files"][chosen], kind="ecp",
            library_family=CANONICAL_FAMILY, library_name=CANONICAL_FAMILY,
            cfour_name=CANONICAL_FAMILY, label_suffix=ECP_LABEL_SUFFIX,
            cfour_label=label, status="PASS", note=""))
        per_element[elem] = {
            "ecp": {"ncore": ncore, "lmax": lmax, "n_explicit_electrons": nexpl,
                    "local_channel": loc,
                    "nonlocal_channels": [l for l, _ in ecp["channels"][1:]],
                    "n_terms_per_channel": {l: len(t) for l, t in ecp["channels"]},
                    "radial_powers": rp,
                    "source_format": chosen,
                    "source_file": meta["ecp_files"][chosen],
                    "spin_orbit_blocks_discarded": so,
                    "channels": [[l, t] for l, t in ecp["channels"]],
                    "validation": status},
            "basis": {}}
        log("ECP  %-3s Z=%-3d NCORE=%-3d Zeff=%-3d LMAX=%d local=%-3s src=.%-8s %s"
            % (elem, z, ncore, nexpl, lmax, loc, chosen, status))

    # ------------------------------ basis sets -------------------------- #
    log("")
    for elem, meta in inv.items():
        z, ccdir = meta["Z"], os.path.join(recipes, elem, CANONICAL_FAMILY)
        ref_txt = clean_reference(meta["author"])
        pe = per_element.setdefault(elem, {"ecp": None, "basis": {}})
        for name, fmts in meta["basis_files"].items():
            checks = []
            try:
                chosen, shells, parsed = select_basis(elem, ccdir, name, fmts,
                                                      args.basis_source, checks,
                                                      bse_version)
                blocks, notes = build_genbas_blocks(shells)
            except Exception as exc:
                failures["basis"].append((elem, name, str(exc)))
                validation_rows.append(dict(element=elem, kind="basis", name=name,
                                            check="source", status="FAIL",
                                            detail=str(exc)))
                continue

            n_src = len(shells)
            n_prim = sum(len(s["exps"]) for s in shells)
            comp = ",".join("%d%s" % (len(b["mat"][0]), char_from_am(b["am"]))
                            for b in blocks)
            pcomp = ",".join("%d%s" % (len(b["exps"]), char_from_am(b["am"]))
                             for b in blocks)
            maxL = max(b["am"] for b in blocks)
            ncols = sum(len(b["mat"][0]) for b in blocks)
            used = sum(1 for b in blocks for row in b["mat"] for v in row if v is not None)

            checks.append(("contraction_count", "PASS" if ncols == n_src else "FAIL",
                           "%d GENBAS columns vs %d source contracted functions"
                           % (ncols, n_src)))
            checks.append(("primitive_count", "PASS" if used == n_prim else "FAIL",
                           "%d matrix entries populated vs %d source primitives"
                           % (used, n_prim)))
            checks.append(("angular_momenta",
                           "PASS" if all(blocks[k]["am"] < blocks[k + 1]["am"]
                                         for k in range(len(blocks) - 1)) else "FAIL",
                           "grouped by L ascending: %s"
                           % ",".join(char_from_am(b["am"]) for b in blocks)))
            degen = [n for n in notes if "NEAR-DEGENERATE" in n]
            checks.append(("exponent_ordering", "WARN" if degen else "PASS",
                           "; ".join(notes) if notes else
                           "source already listed every L in strictly descending "
                           "exponent order"))
            checks.append(("high_L_present", "PASS",
                           "highest angular momentum retained: %s (L=%d)"
                           % (char_from_am(maxL), maxL)))

            comment = ("ccECP %s for %s | source %s | %s"
                       % (name, elem, os.path.join(meta["recipe_dir"], fmts[chosen]),
                          ref_txt))
            entry = format_genbas_entry(elem, name, blocks, comment,
                                        args.basis_suffix)
            try:
                rlabel, rblocks = parse_genbas_entry(entry)
                ok, msg = compare_roundtrip(blocks, rblocks)
                good = ok and rlabel == "%s:%s" % (
                    elem, basis_label_name(name, args.basis_suffix))
                checks.append(("roundtrip_exact", "PASS" if good else "FAIL",
                               msg or "every exponent and coefficient re-reads "
                                      "with exact Decimal equality"))
            except Exception as exc:
                checks.append(("roundtrip_exact", "FAIL", str(exc)))

            status = "PASS" if all(s.startswith("PASS") or s == "WARN"
                                   for _, s, _ in checks) else "FAIL"
            for c, s, d in checks:
                validation_rows.append(dict(element=elem, kind="basis", name=name,
                                            check=c, status=s, detail=d))
            if status == "FAIL":
                failures["basis"].append((elem, name, "; ".join(
                    "%s=%s" % (c, s) for c, s, _ in checks
                    if not (s.startswith("PASS") or s == "WARN"))))
                continue

            genbas_entries.append((elem, name, entry))
            label = "%s:%s" % (elem,
                               basis_label_name(name, args.basis_suffix))
            seen_labels[label] += 1
            src_used["basis"][chosen] += 1
            basis_map_rows.append(dict(
                element=elem, atomic_number=z,
                source_recipe_dir=meta["recipe_dir"],
                source_filename=fmts[chosen], source_format=chosen,
                all_source_formats=";".join(sorted(fmts)),
                precision_digits_per_format=";".join(
                    "%s:%d" % (f, basis_precision(v)) for f, v in sorted(parsed.items())),
                library_family=CANONICAL_FAMILY, source_name=name,
                cfour_name=basis_label_name(name, args.basis_suffix),
                label_suffix=args.basis_suffix, cfour_label=label,
                n_contracted=ncols, n_primitives=n_prim,
                max_L=char_from_am(maxL),
                contracted_composition=comp, primitive_composition=pcomp,
                name_transformation=("source library name plus declared suffix %r"
                                     % args.basis_suffix),
                validation=status))
            naming_rows.append(dict(
                element=elem, source_recipe_dir=meta["recipe_dir"],
                source_filename=fmts[chosen], kind="basis",
                library_family=CANONICAL_FAMILY, library_name=name,
                cfour_name=basis_label_name(name, args.basis_suffix),
                label_suffix=args.basis_suffix,
                cfour_label=label, status="PASS", note=""))
            pe["basis"][name] = dict(source_format=chosen, source_file=fmts[chosen],
                                     n_contracted=ncols, n_primitives=n_prim,
                                     max_L=char_from_am(maxL),
                                     contracted_composition=comp,
                                     primitive_composition=pcomp,
                                     validation=status)
        log("BASIS %-3s %2d sets  %s" % (elem, len(pe["basis"]),
                                         " ".join(sorted(pe["basis"]))))

    # ------------------------- hierarchy checks ------------------------ #
    log("")
    hier = hierarchy_checks(inv, recipes, args.basis_source, validation_rows,
                            anomalies, log)

    # --------------------------- naming audit -------------------------- #
    for row in naming_rows:
        bad = [f for f in FORBIDDEN_FAMILIES
               if re.search(r"(^|[/.])%s([/.]|$)" % re.escape(f),
                            row["source_recipe_dir"] + "/" + row["source_filename"])]
        if bad:
            row["status"], row["note"] = "FAIL", "non-canonical family %r" % bad
        elif row["cfour_name"] != row["library_name"] + row["label_suffix"]:
            row["status"] = "FAIL"
            row["note"] = ("name changed %r -> %r, which is not the library "
                           "name plus the declared suffix %r"
                           % (row["library_name"], row["cfour_name"],
                              row["label_suffix"]))
        elif row["cfour_label"] != "%s:%s" % (row["element"], row["cfour_name"]):
            row["status"], row["note"] = "FAIL", "unexpected label form"
        else:
            row["status"] = "PASS"
            if row["label_suffix"]:
                row["note"] = ("source library name retained as the label stem; "
                               "declared disambiguating suffix %r appended: %s -> %s"
                               % (row["label_suffix"], row["library_name"],
                                  row["cfour_name"]))
            else:
                row["note"] = ("no transformation: label is exactly "
                               "<library element>:<library name>")

    dup = {k: v for k, v in seen_labels.items() if v > 1}

    # ------------------- auxiliary (RI) basis conversion ---------------- #
    #
    # Auxiliary sets are density-fitting bases belonging to one orbital basis.
    # They are converted through exactly the same path as the orbital sets --
    # same parser body, same block builder, same renderer -- so that no second
    # numerical code path exists.  Only the label differs, gaining `-RI` before
    # the suffix.  Nothing here is written into GENBAS or GENBAS.mrcc; the RI
    # entries go to their own files so the validated orbital library stays
    # byte-for-byte as it was.
    extra_aux = {}
    if args.extra_aux and not args.no_ri:
        adir = os.path.abspath(os.path.expanduser(args.extra_aux))
        if not os.path.isdir(adir):
            raise ConversionError("--extra-aux is not a directory: %s" % adir)
        for f in sorted(os.listdir(adir)):
            m = re.fullmatch(r"([A-Z][a-z]?)\.(.+)\.%s"
                             % re.escape(AUX_FORMAT), f)
            if m is None:
                continue
            extra_aux[(m.group(1), m.group(2))] = os.path.join(adir, f)
        log("extra auxiliary sets supplied  : %d from %s"
            % (len(extra_aux), adir))

    if not args.no_ri:
        for elem, meta in inv.items():
            ccdir = os.path.join(recipes, elem, CANONICAL_FAMILY)
            for name in sorted(meta["basis_files"]):
                origin, path = "library", None
                if name in meta["aux_files"]:
                    path = os.path.join(ccdir, meta["aux_files"][name])
                elif (elem, name) in extra_aux:
                    origin, path = "extra", extra_aux[(elem, name)]
                if path is None:
                    # The library ships no auxiliary set for this orbital
                    # basis.  Nothing is invented; it is recorded and skipped.
                    anomalies.append(dict(
                        element=elem, kind="ri_basis", name=name,
                        check="aux_present", status="INFO",
                        detail="no %s file in the library and none supplied "
                               "via --extra-aux; no auxiliary entry written"
                               % AUX_FORMAT))
                    validation_rows.append(dict(
                        element=elem, kind="ri_basis", name=name,
                        check="aux_present", status="INFO",
                        detail="absent upstream; not invented"))
                    continue

                checks = []
                try:
                    with open(path) as fh:
                        shells = parse_autoauxc(fh.read(), elem)
                except (ConversionError, OSError) as exc:
                    failures["ri"].append((elem, name, str(exc)))
                    validation_rows.append(dict(
                        element=elem, kind="ri_basis", name=name,
                        check="parse", status="FAIL", detail=str(exc)))
                    continue
                checks.append(("parse", "PASS",
                               "%d shells read from %s"
                               % (len(shells), os.path.basename(path))))
                unc = aux_is_uncontracted(shells)
                checks.append(("uncontracted", "PASS" if unc else "WARN",
                               "every shell one primitive with coefficient 1"
                               if unc else
                               "set is contracted; converted unchanged"))
                try:
                    blocks, notes = build_genbas_blocks(shells)
                except ConversionError as exc:
                    failures["ri"].append((elem, name, str(exc)))
                    validation_rows.append(dict(
                        element=elem, kind="ri_basis", name=name,
                        check="blocks", status="FAIL", detail=str(exc)))
                    continue
                for n in notes:
                    checks.append(("layout", "WARN", n))
                ri_name = ri_label_name(name, args.basis_suffix)
                comment = ("! ORCA AutoAux /C auxiliary basis for %s %s"
                           " | source %s | %s"
                           % (elem, name, os.path.basename(path),
                              "canonical ccECP library"
                              if origin == "library"
                              else "generated locally with ORCA AutoAux, "
                                   "absent from the library"))
                entry = format_genbas_entry(elem, ri_name, blocks, comment, "")
                try:
                    _, rblocks = parse_genbas_entry(entry)
                    good, msg = compare_roundtrip(blocks, rblocks)
                    checks.append(("roundtrip_exact", "PASS" if good else "FAIL",
                                   msg or "every exponent and coefficient "
                                          "re-reads with exact Decimal equality"))
                except Exception as exc:
                    checks.append(("roundtrip_exact", "FAIL", str(exc)))

                status = "PASS" if all(c[1].startswith("PASS") or c[1] == "WARN"
                                       for c in checks) else "FAIL"
                for c, st, d in checks:
                    validation_rows.append(dict(element=elem, kind="ri_basis",
                                                name=name, check=c, status=st,
                                                detail=d))
                if status == "FAIL":
                    failures["ri"].append((elem, name, "; ".join(
                        "%s=%s" % (c, st) for c, st, _ in checks
                        if not (st.startswith("PASS") or st == "WARN"))))
                    continue

                label = "%s:%s" % (elem, ri_name)
                ri_entries.append((elem, name, entry))
                seen_labels[label] += 1
                maxL = max(b["am"] for b in blocks)
                ri_map_rows.append(dict(
                    element=elem, atomic_number=ELEMENT_Z[elem],
                    source_recipe_dir=meta["recipe_dir"],
                    source_filename=os.path.basename(path),
                    source_format=AUX_FORMAT, source_origin=origin,
                    library_family=CANONICAL_FAMILY,
                    orbital_basis=name, source_name=name,
                    cfour_name=ri_name, ri_infix=RI_NAME_SUFFIX,
                    label_suffix=args.basis_suffix, cfour_label=label,
                    n_functions=len(shells),
                    n_primitives=sum(len(sh["exps"]) for sh in shells),
                    max_L=char_from_am(maxL),
                    composition=" ".join(
                        "%s%d" % (char_from_am(b["am"]), len(b["exps"]))
                        for b in blocks),
                    uncontracted="yes" if unc else "no",
                    name_transformation=("orbital name plus %r plus declared "
                                         "suffix %r"
                                         % (RI_NAME_SUFFIX, args.basis_suffix))))
                naming_rows.append(dict(
                    element=elem, source_recipe_dir=meta["recipe_dir"],
                    source_filename=os.path.basename(path), kind="ri_basis",
                    library_family=CANONICAL_FAMILY, library_name=name,
                    cfour_name=ri_name, label_suffix=args.basis_suffix,
                    cfour_label=label, status="PASS", note=""))
        dup = {k: v for k, v in seen_labels.items() if v > 1}

    # ---------------------------- write output ------------------------- #
    hdr = dict(stamp=stamp, commit=commit, commit_date=commit_date,
               bse=bse_version, basis_src=args.basis_source,
               ecp_src=args.ecp_source, basis_suffix=args.basis_suffix)
    gpath = os.path.join(outdir, "GENBAS")
    with open(gpath, "w") as fh:
        fh.write(file_header("GENBAS", hdr, "!"))
        fh.write("\n")
        for _, _, txt in genbas_entries:
            fh.write(txt)
    epath = os.path.join(outdir, "ECPDATA")
    with open(epath, "w") as fh:
        fh.write(file_header("ECPDATA", hdr, "#"))
        for _, _, txt in ecpdata_entries:
            fh.write(txt)
        fh.write("*\n")
    with open(os.path.join(outdir, "ECPDATA.maxprec"), "w") as fh:
        fh.write(file_header("ECPDATA (maximum-precision variant)", hdr, "#"))
        fh.write("# For every element the analytic rendering carrying the most\n"
                 "# significant digits was used here, instead of the uniform\n"
                 "# rendering used in ECPDATA.  All renderings were verified to\n"
                 "# agree to within the written precision of the coarser one.\n")
        for _, _, txt in ecpdata_best_entries:
            fh.write(txt)
        fh.write("*\n")
    with open(os.path.join(outdir, "GENBAS.mrcc"), "w") as fh:
        fh.write(file_header("GENBAS (MRCC variant: ECPs appended)", hdr, "!"))
        fh.write("\n")
        for _, _, txt in genbas_entries:
            fh.write(txt)
        fh.write("\n! Effective core potentials -- canonical ccECP\n")
        for _, _, txt in ecpdata_entries:
            fh.write(txt)
        fh.write("*\n")

    # Measured, not assumed: how many of our labels would clash with the
    # target system's own library, with and without the suffix.
    collide = {}
    if args.host_basis_library:
        host = scan_host_labels(
            os.path.abspath(os.path.expanduser(args.host_basis_library)))
        ours_orb = {"%s:%s" % (e, basis_label_name(n, args.basis_suffix))
                    for e, n, _ in genbas_entries}
        ours_ri = {"%s:%s" % (e, ri_label_name(n, args.basis_suffix))
                   for e, n, _ in ri_entries}
        bare_orb = {"%s:%s" % (e, n) for e, n, _ in genbas_entries}
        bare_ri = {"%s:%s%s" % (e, n, RI_NAME_SUFFIX)
                   for e, n, _ in ri_entries}
        collide = dict(host_labels=len(host),
                       orbital=len(host & ours_orb),
                       orbital_bare=len(host & bare_orb),
                       ri=len(host & ours_ri),
                       ri_bare=len(host & bare_ri))
        log("host basis library             : %s" % args.host_basis_library)
        log("  labels in host library       : %d" % collide["host_labels"])
        log("  collisions, orbital labels   : %d (%d without the suffix)"
            % (collide["orbital"], collide["orbital_bare"]))
        log("  collisions, auxiliary labels : %d (%d without the suffix)"
            % (collide["ri"], collide["ri_bare"]))

    ri_bases = []
    if ri_entries:
        # RI entries on their own, for a job that wants only the fitting sets.
        with open(os.path.join(outdir, "GENBAS.RI"), "w") as fh:
            fh.write(file_header("GENBAS (auxiliary/RI sets only)", hdr, "!"))
            fh.write("\n")
            for _, _, txt in ri_entries:
                fh.write(txt)
        # Orbital + RI, for CFOUR.
        with open(os.path.join(outdir, "GENBAS.withRI"), "w") as fh:
            fh.write(file_header("GENBAS (orbital + auxiliary/RI sets)",
                                 hdr, "!"))
            fh.write("\n")
            for _, _, txt in genbas_entries:
                fh.write(txt)
            fh.write("\n! Auxiliary (RI) basis sets -- ORCA AutoAux /C\n")
            for _, _, txt in ri_entries:
                fh.write(txt)
        # Orbital + RI + ECP in one file, which is what an MRCC job copies to
        # ./GENBAS when it uses dfbasis_cor.
        with open(os.path.join(outdir, "GENBAS.mrcc.withRI"), "w") as fh:
            fh.write(file_header("GENBAS (MRCC variant: orbital + auxiliary/RI"
                                 " + ECPs appended)", hdr, "!"))
            fh.write("\n")
            for _, _, txt in genbas_entries:
                fh.write(txt)
            fh.write("\n! Auxiliary (RI) basis sets -- ORCA AutoAux /C\n")
            for _, _, txt in ri_entries:
                fh.write(txt)
            fh.write("\n! Effective core potentials -- canonical ccECP\n")
            for _, _, txt in ecpdata_entries:
                fh.write(txt)
            fh.write("*\n")
        ri_bases = ["GENBAS.RI", "GENBAS.withRI", "GENBAS.mrcc.withRI"]
        write_csv(os.path.join(outdir, "ri_name_map.csv"), ri_map_rows)

    for base in ("GENBAS", "ECPDATA", "GENBAS.mrcc", "ECPDATA.maxprec") \
            + tuple(ri_bases):
        src_path = os.path.join(outdir, base)
        with open(src_path) as fh:
            body = fh.read()
        with open(src_path + ".upper", "w") as fh:
            fh.write(upcase_labels(body))

    write_csv(os.path.join(outdir, "basis_name_map.csv"), basis_map_rows)
    write_csv(os.path.join(outdir, "ecp_name_map.csv"), ecp_map_rows)
    write_csv(os.path.join(outdir, "naming_audit.csv"), naming_rows)
    write_csv(os.path.join(outdir, "validation_report.csv"), validation_rows)
    for r in validation_rows:
        if r["status"] in ("FAIL", "WARN"):
            anomalies.append(dict(element=r["element"], kind=r["kind"],
                                  name=r["name"], check=r["check"],
                                  severity="fail" if r["status"] == "FAIL" else "warn",
                                  detail=r["detail"]))
    write_csv(os.path.join(outdir, "anomalies.csv"), anomalies)

    inv_rows = []
    for elem, meta in inv.items():
        pe = per_element.get(elem) or {}
        ecp = pe.get("ecp") or {}
        names = sorted(pe.get("basis", {}))
        anomalies = []
        if meta["skipped_files"]:
            anomalies.append("library files present but not used: %s"
                             % ";".join(sorted(meta["skipped_files"])))
        if ecp.get("ncore") == 0:
            anomalies.append("NCORE=0: all electrons are explicit; the ccECP acts "
                             "purely as a shape-consistent potential")
        miss = sorted(set(meta["basis_files"]) - set(names))
        if miss:
            anomalies.append("basis sets in library but NOT converted: %s" % ";".join(miss))
        if not ecp:
            anomalies.append("ECP conversion FAILED")
        inv_rows.append(dict(
            element=elem, atomic_number=meta["Z"], recipe_dir=meta["recipe_dir"],
            library_family=CANONICAL_FAMILY,
            ecp_source_files=";".join("%s=%s" % kv for kv in
                                      sorted(meta["ecp_files"].items())),
            ecp_source_used=ecp.get("source_file"),
            ncore=ecp.get("ncore"), n_explicit_electrons=ecp.get("n_explicit_electrons"),
            lmax=ecp.get("lmax"), local_ecp_channel=ecp.get("local_channel"),
            nonlocal_ecp_channels=";".join(ecp.get("nonlocal_channels", [])),
            n_basis_sets=len(names), basis_names=";".join(names),
            basis_source_formats=";".join(sorted(
                {f for v in meta["basis_files"].values() for f in v})),
            ecp_validation=ecp.get("validation") or "FAIL",
            anomalies=" | ".join(anomalies)))
    write_csv(os.path.join(outdir, "inventory.csv"), inv_rows)

    with open(os.path.join(outdir, "inventory.json"), "w") as fh:
        json.dump({"generated": stamp, "converter_version": VERSION,
                   "library_commit": commit, "library_commit_date": commit_date,
                   "library_path": args.library,
                   "python": sys.version.split()[0],
                   "basis_set_exchange": bse_version,
                   "family": CANONICAL_FAMILY,
                   "selection_rule": "recipes/<ELEMENT>/ccECP/ only",
                   "excluded_recipe_families": dict(excluded),
                   "basis_source_selection": args.basis_source,
                   "ecp_source_selection": args.ecp_source,
                   "source_formats_used": {k: dict(v) for k, v in src_used.items()},
                   "convention_checks": [{"reference": k, "status": s, "detail": d}
                                          for k, s, d in convention_results],
                   "elements": {e: {"element": e, "Z": m["Z"],
                                    "recipe_dir": m["recipe_dir"],
                                    "reference": clean_reference(m["author"]),
                                    "ecp_files": m["ecp_files"],
                                    "basis_files": m["basis_files"],
                                    "unused_files": m["skipped_files"],
                                    "converted": per_element.get(e)}
                                 for e, m in inv.items()}},
                  fh, indent=1, sort_keys=True)

    if not args.no_by_element:
        grouped = collections.OrderedDict()
        for e, n, t in genbas_entries:
            grouped.setdefault(e, {"basis": "", "ecp": "", "files": []})["basis"] += t
        for e, n, t in ecpdata_entries:
            grouped.setdefault(e, {"basis": "", "ecp": "", "files": []})["ecp"] += t
        for e, n, t in ri_entries:
            grouped.setdefault(e, {"basis": "", "ecp": "", "files": []})
            grouped[e].setdefault("ri", "")
            grouped[e]["ri"] += t

        # one file per (element, object), using the source logical object name
        # with a `.cfour` suffix:
        #     Fe.ccECP.gamess   ->   Fe.ccECP.cfour     (one ECPDATA entry)
        #     Fe.cc-pVQZ.gamess ->   Fe.cc-pVQZ.cfour   (one GENBAS entry)
        for e, n, t in ecpdata_entries:
            dd = os.path.join(outdir, "by_element", e)
            os.makedirs(dd, exist_ok=True)
            fn = "%s.%s.cfour" % (e, n)
            with open(os.path.join(dd, fn), "w") as fh:
                fh.write(single_file_header(e, n, "ECPDATA", hdr))
                fh.write(t)
                fh.write("*\n")
            grouped[e]["files"].append(fn)
            ecp_file_of[(e, n)] = os.path.join("by_element", e, fn)
        for e, n, t in genbas_entries:
            dd = os.path.join(outdir, "by_element", e)
            os.makedirs(dd, exist_ok=True)
            fn = "%s.%s.cfour" % (e, n)
            with open(os.path.join(dd, fn), "w") as fh:
                fh.write(single_file_header(e, n, "GENBAS", hdr))
                fh.write("\n")
                fh.write(t)
            grouped[e]["files"].append(fn)
            basis_file_of[(e, basis_label_name(n, args.basis_suffix))] = (
                os.path.join("by_element", e, fn))
        #     Fe.cc-pVQZ.AutoAuxC.orca -> Fe.cc-pVQZ-RI.cfour  (one RI entry)
        for e, n, t in ri_entries:
            dd = os.path.join(outdir, "by_element", e)
            os.makedirs(dd, exist_ok=True)
            fn = "%s.%s%s.cfour" % (e, n, RI_NAME_SUFFIX)
            with open(os.path.join(dd, fn), "w") as fh:
                fh.write(single_file_header(e, n + RI_NAME_SUFFIX, "GENBAS",
                                            hdr))
                fh.write("\n")
                fh.write(t)
            grouped[e]["files"].append(fn)
            ri_file_of[(e, ri_label_name(n, args.basis_suffix))] = (
                os.path.join("by_element", e, fn))

        for e, d in grouped.items():
            dd = os.path.join(outdir, "by_element", e)
            os.makedirs(dd, exist_ok=True)
            if d["basis"]:
                with open(os.path.join(dd, "GENBAS"), "w") as fh:
                    fh.write("\n" + d["basis"])
            if d.get("ri"):
                with open(os.path.join(dd, "GENBAS.RI"), "w") as fh:
                    fh.write("\n" + d["ri"])
            if d["ecp"]:
                with open(os.path.join(dd, "ECPDATA"), "w") as fh:
                    fh.write(d["ecp"] + "*\n")
            with open(os.path.join(dd, "summary.json"), "w") as fh:
                json.dump(per_element.get(e), fh, indent=1, sort_keys=True)

        for r in ri_map_rows:
            r["cfour_file"] = ri_file_of.get((r["element"], r["cfour_name"]), "")
        for r in basis_map_rows:
            r["cfour_file"] = basis_file_of.get((r["element"], r["cfour_name"]), "")
        for r in ecp_map_rows:
            r["cfour_file"] = ecp_file_of.get((r["element"], r["cfour_name"]), "")
        for r in naming_rows:
            key = (r["element"], r["cfour_name"])
            r["cfour_file"] = (ecp_file_of if r["kind"] == "ecp"
                               else basis_file_of).get(key, "")
        write_csv(os.path.join(outdir, "basis_name_map.csv"), basis_map_rows)
        write_csv(os.path.join(outdir, "ecp_name_map.csv"), ecp_map_rows)
        if ri_map_rows:
            write_csv(os.path.join(outdir, "ri_name_map.csv"), ri_map_rows)
        write_csv(os.path.join(outdir, "naming_audit.csv"), naming_rows)

    summary = dict(
        elements_found=len(inv),
        ecpdata_entries=len(ecpdata_entries), genbas_entries=len(genbas_entries),
        ri_entries=len(ri_entries),
        ri_from_library=sum(1 for r in ri_map_rows
                            if r["source_origin"] == "library"),
        ri_from_extra=sum(1 for r in ri_map_rows
                          if r["source_origin"] == "extra"),
        ri_missing=sum(1 for r in validation_rows
                       if r["kind"] == "ri_basis" and r["check"] == "aux_present"),
        collisions=collide or None,
        distinct_basis_names=sorted({r["source_name"] for r in basis_map_rows}),
        distinct_basis_labels=sorted({r["cfour_name"] for r in basis_map_rows}),
        naming_failures=[r for r in naming_rows if r["status"] != "PASS"],
        duplicate_labels=dup, failed_ecp=failures["ecp"], failed_basis=failures["basis"],
        convention_checks=convention_results, hierarchy=hier,
        source_formats_used={k: dict(v) for k, v in src_used.items()},
        library_commit=commit, library_commit_date=commit_date, generated=stamp,
        library_path=args.library, output_path=outdir,
        python=sys.version.split()[0], bse=bse_version,
        basis_source_selection=args.basis_source, ecp_source_selection=args.ecp_source,
        excluded_recipe_families=dict(excluded),
        elements_without_ccECP=sorted(set(os.listdir(recipes)) - set(inv)
                                      - {f for f in os.listdir(recipes)
                                         if not os.path.isdir(os.path.join(recipes, f))}),
        validation_counts=dict(collections.Counter(r["status"] for r in validation_rows)))
    with open(os.path.join(outdir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1, default=str)

    write_readme(outdir, summary, inv, per_element, inv_rows, basis_map_rows,
                 ecp_map_rows, validation_rows)

    log("")
    log("=" * 78)
    log("Canonical ccECP elements found : %d" % len(inv))
    log("Successfully converted (ECP)   : %d" % len(ecpdata_entries))
    log("Failed (ECP)                   : %d %s"
        % (len(failures["ecp"]), failures["ecp"] or ""))
    log("GENBAS entries                 : %d" % len(genbas_entries))
    log("auxiliary (RI) entries         : %d" % len(ri_entries))
    if failures["ri"]:
        log("auxiliary (RI) FAILURES        : %d" % len(failures["ri"]))
        for e, n, why in failures["ri"][:10]:
            log("    %s %s: %s" % (e, n, why))
    log("Failed (basis)                 : %d %s"
        % (len(failures["basis"]), failures["basis"] or ""))
    log("Distinct basis-set names       : %d" % len(summary["distinct_basis_names"]))
    log("Duplicate ELEM:NAME labels     : %d" % len(dup))
    log("Naming-audit failures          : %d" % len(summary["naming_failures"]))
    log("Validation rows                : %d %r"
        % (len(validation_rows), summary["validation_counts"]))
    log("Convention-check failures      : %d" % len(convention_failures))
    log("NCORE/LMAX fixed-column tests  : %d/%d PASS (cols 12-14 and 26)"
        % (len(ncore_lmax_tests), len(ncore_lmax_tests)))
    log("Source formats used (ECP)      : %r" % dict(src_used["ecp"]))
    log("Source formats used (basis)    : %r" % dict(src_used["basis"]))
    log("=" * 78)
    logfh.close()
    return 0 if (not failures["ecp"] and not failures["basis"]
                 and not summary["naming_failures"] and not dup
                 and not convention_failures) else 1


# --------------------------------------------------------------------------- #
#  Cross-basis scientific consistency
# --------------------------------------------------------------------------- #

def hierarchy_checks(inv, recipes, basis_fmt, validation_rows, anomalies, log):
    """Diagnostic consistency checks between related basis-set families.

    The checks report whether an `aug-X` set contains its parent `X` primitives
    and whether a `cc-pCVnZ` set contains its `cc-pVnZ` parent.  Some published
    families are independently optimized rather than strict supersets, so a
    mismatch is reported as INFO for review and is not treated as a conversion
    failure.  Primitives are matched with `consistent()` to allow harmless
    differences in written precision between source renderings.
    """
    cache = {}

    def load(elem, name):
        key = (elem, name)
        if key in cache:
            return cache[key]
        fmts = inv[elem]["basis_files"].get(name)
        got = None
        if fmts:
            fmt = (basis_fmt if basis_fmt != "auto"
                   else ("gamess" if "gamess" in fmts else sorted(fmts)[0]))
            if fmt in fmts and fmt in BASIS_PARSERS:
                with open(os.path.join(recipes, elem, CANONICAL_FAMILY,
                                       fmts[fmt])) as fh:
                    got = BASIS_PARSERS[fmt](fh.read(), elem)
        cache[key] = got
        return got

    def expsets(shells):
        d = collections.defaultdict(list)
        for sh in shells:
            for e in sh["exps"]:
                d[sh["am"]].append(e)
        return d

    stats = collections.Counter()
    for elem in inv:
        for name in sorted(inv[elem]["basis_files"]):
            if name.startswith("aug-"):
                base, kind = name[4:], "aug"
            elif "cc-pCV" in name:
                base, kind = name.replace("cc-pCV", "cc-pV"), "cv"
            else:
                continue
            check = "aug_vs_parent" if kind == "aug" else "core_valence_vs_parent"
            if base not in inv[elem]["basis_files"]:
                validation_rows.append(dict(element=elem, kind="basis", name=name,
                                            check=check, status="NA",
                                            detail="parent %s is not in the "
                                                   "canonical library for this "
                                                   "element" % base))
                stats["NA"] += 1
                continue
            a, b = load(elem, name), load(elem, base)
            if a is None or b is None:
                stats["NA"] += 1
                continue
            ea, eb = expsets(a), expsets(b)
            na = collections.Counter(sh["am"] for sh in a)
            nb = collections.Counter(sh["am"] for sh in b)
            comp = lambda c: ",".join("%d%s" % (c[l], char_from_am(l))
                                      for l in sorted(c))

            missing, extreme_ok, notes = {}, True, []
            for l in sorted(eb):
                for e in eb[l]:
                    if not any(consistent(e, f)[0] for f in ea.get(l, ())):
                        missing.setdefault(char_from_am(l), []).append(e)
                if l not in ea:
                    extreme_ok = False
                    continue
                if kind == "aug":
                    if min(dec(x) for x in ea[l]) >= min(dec(x) for x in eb[l]):
                        extreme_ok = False
                        notes.append("FYI no primitive more diffuse than the "
                                     "parent's smallest in %s" % char_from_am(l))
                else:
                    if max(dec(x) for x in ea[l]) <= max(dec(x) for x in eb[l]):
                        notes.append("FYI no primitive tighter than the parent's "
                                     "largest in %s" % char_from_am(l))
            grew = all(na[l] >= nb[l] for l in nb)
            retains = not missing
            # NOTE: whether the *added* function is more diffuse (aug) or tighter
            # (core-valence) than every parent primitive is not a well-posed
            # criterion, because a parent contraction can already contain a very
            # diffuse or very tight primitive.  It is reported in the detail text
            # only and does not decide PASS/INFO.
            ok = retains and grew
            detail = ("%s [%s] keeps every primitive of %s [%s]"
                      % (name, comp(na), base, comp(nb)))
            if not ok:
                detail = ("%s [%s] vs %s [%s]: primitives of the parent that are "
                          "absent: %r%s%s"
                          % (name, comp(na), base, comp(nb), missing,
                             "" if grew else " | fewer functions in some L",
                             (" | " + "; ".join(notes)) if notes else ""))
            status = "PASS" if ok else "INFO"
            validation_rows.append(dict(element=elem, kind="basis", name=name,
                                        check=check, status=status, detail=detail))
            if not ok:
                anomalies.append(dict(element=elem, kind="basis", name=name,
                                      check=check, severity="review",
                                      detail=detail))
            stats[status] += 1
    log("cross-basis hierarchy checks: %r" % dict(stats))
    return dict(stats)


# --------------------------------------------------------------------------- #
#  Output helpers
# --------------------------------------------------------------------------- #

def single_file_header(elem, name, kind, hdr):
    """Header for a one-object file, e.g. `Fe.cc-pVQZ.cfour`."""
    cc = "!" if kind == "GENBAS" else "#"
    label_name = (basis_label_name(name, hdr.get("basis_suffix", ""))
                  if kind == "GENBAS" else name)
    return ("%s %s:%s -- canonical ccECP data from PseudopotentialLibrary.org\n"
            "%s CFOUR/MRCC %s format; append to your %s file, or use as is.\n"
            "%s library commit %s | generated %s by convert_ccecp_to_cfour.py v%s\n"
            % (cc, elem, label_name, cc, kind, kind, cc, hdr["commit"],
               hdr["stamp"], VERSION))


def file_header(kind, hdr, cc):
    text = """\
---------------------------------------------------------------------------
CFOUR/MRCC {kind}
canonical ccECP data from PseudopotentialLibrary.org
---------------------------------------------------------------------------
This file contains ONLY the canonical ccECP family, i.e. exactly the entries
found under  recipes/<ELEMENT>/ccECP/  in the Pseudopotential Library
repository.  No other recipe family (ccECP-soft, soft-ccECP-deprecated,
eCEPP, CEPP, BFD, RRKJ, TM, ...) is present, and no spin-orbit (SOREP) data
is present: only the scalar-relativistic canonical ccECP is converted.

Naming: ECP entries use  <Element>:ccECP.  Orbital-basis entries use
<Element>:<SOURCE BASIS NAME>{suffix}.  The optional basis-label suffix prevents
collisions with installed basis libraries; it does not alter any numerical data.
Optional *.upper companion files uppercase only the element symbol.

ECP convention:  U_l(r) = sum_m c_m r**(N_m - 2) exp(-alpha_m r**2)
                 columns are   coefficient   N   alpha
                 local channel first (labelled l = LMAX), then s-<l>, p-<l>...

Pseudopotential Library commit : {commit}  ({cdate})
generated                      : {stamp}
generator                      : convert_ccecp_to_cfour.py v{ver}
python                         : {py}
basis_set_exchange             : {bse}
basis source selection         : {bsrc}
ECP   source selection         : {esrc}

The orbital basis sets are pure spherical-harmonic sets (5d/7f/9g/...):
run CFOUR with SPHERICAL=ON and MRCC without cartesian d functions.
---------------------------------------------------------------------------
""".format(kind=kind, commit=hdr["commit"], cdate=hdr["commit_date"],
           stamp=hdr["stamp"], ver=VERSION, py=sys.version.split()[0],
           bse=hdr["bse"] or "not used", bsrc=hdr["basis_src"],
           esrc=hdr["ecp_src"], suffix=hdr.get("basis_suffix", ""))
    return "".join("%s %s\n" % (cc, l) if l else "%s\n" % cc
                   for l in text.splitlines())


def upcase_labels(text):
    """Return the same file with the element symbol of every entry label in
    upper case (`Ag:cc-pVQZ` -> `AG:cc-pVQZ`).

    Some CFOUR builds fold the atom string taken from the ZMAT to upper case
    before searching GENBAS/ECPDATA, and the GENBAS shipped with CFOUR is
    written that way.  This is the *only* naming transformation this converter
    ever performs, it is confined to the optional `*.upper` files, and it
    touches nothing but the element symbol in front of the colon.
    """
    return re.sub(r"(?m)^([A-Z][a-z]?):(\S+)\s*$",
                  lambda m: "%s:%s" % (m.group(1).upper(), m.group(2)), text)


def write_csv(path, rows):
    if not rows:
        open(path, "w").close()
        return
    keys = list(rows[0].keys())
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_readme(outdir, s, inv, per_element, inv_rows, basis_rows, ecp_rows,
                 val_rows):
    L = []
    A = L.append
    # Naming suffix actually used for this conversion.
    s_suffix = basis_rows[0]["label_suffix"] if basis_rows else ""
    A("# Canonical ccECP for CFOUR / MRCC")
    A("")
    A("`GENBAS` and `ECPDATA` generated from the **canonical `ccECP`** family of")
    A("[PseudopotentialLibrary.org](https://pseudopotentiallibrary.org).")
    A("")
    A("> **This output represents ONLY the canonical ccECP family from")
    A("> PseudopotentialLibrary.org.**  No `ccECP-soft`, `soft-ccECP-deprecated`,")
    A("> `eCEPP`, `CEPP`, `BFD`, `RRKJ` or `TM` data is included, and no")
    A("> grid-based AREP/SOREP files are included; spin-orbit SOREP data is")
    A("> excluded -- the conversion covers the scalar-relativistic ccECP only.")
    A("")
    A("## Quick start")
    A("")
    A("The converter is a standalone Python script. CFOUR and MRCC are **not**")
    A("required to generate the files. Python 3.8+ is sufficient; the optional")
    A("`basis_set_exchange` package is strongly recommended for independent")
    A("cross-checks.")
    A("")
    A("```bash")
    A("git clone https://github.com/QMCPACK/pseudopotentiallibrary.git")
    A("python3 -m venv .venv")
    A("source .venv/bin/activate")
    A("python -m pip install basis_set_exchange")
    A("python convert_ccecp_to_cfour.py \\")
    A("  --library ./pseudopotentiallibrary \\")
    A("  --out ./ccECP_cfour_output")
    A("```")
    A("")
    A("On Windows, run the same commands in WSL, Git Bash, or another Python")
    A("environment. `--library` may point to the repository root or directly to")
    A("its `recipes` directory. If it is omitted, the script checks")
    A("`$PSEUDOPOTENTIAL_LIBRARY` and local `pseudopotentiallibrary` directories.")
    A("Use `--no-bse` only when Basis Set Exchange is unavailable.")
    A("")
    A("By default, orbital-basis labels receive the suffix `%s` to distinguish" % s_suffix)
    A("the ccECP-optimized basis from installed basis sets with the same conventional")
    A("name. ECP labels remain exactly `ccECP`. See section 6 for details.")
    A("")
    A("The source repository is treated as read-only; choose an output directory")
    A("outside the Pseudopotential Library checkout.")
    A("")
    A("## 1. Source and validation details")
    A("")
    A("| item | value |")
    A("|---|---|")
    A("| Pseudopotential Library commit | `%s` |" % s["library_commit"])
    A("| commit date | %s |" % s["library_commit_date"])
    A("| conversion date | %s |" % s["generated"])
    A("| converter | `convert_ccecp_to_cfour.py` v%s |" % VERSION)
    A("| Python | %s |" % s["python"])
    A("| basis_set_exchange | %s (independent cross-check only) |" % (s["bse"] or "not used"))
    A("")
    A("## 2. Canonical ccECP selection rule")
    A("")
    A("A recipe is in scope **iff** it lives in `recipes/<ELEMENT>/ccECP/`.")
    A("Supported **basis** renderings are `.gamess`, `.nwchem` and `.gaussian`.")
    A("Supported **ECP** renderings are the native `<El>.ccECP` text plus")
    A("`.gamess`, `.nwchem`, `.gaussian` and `.molpro`. Molpro basis files and")
    A("all `.dirac`, `.xml`, `.upf`, `.rpt`, `*_deprecated`, `.AREP.*` and")
    A("`.SOREP.*` files are deliberately ignored. No numerical-grid or spin-orbit")
    A("data can therefore leak into the scalar-relativistic output.")
    A("")
    A("Recipe families that exist in the library and are **excluded**:")
    A("")
    for fam, n in sorted(s["excluded_recipe_families"].items()):
        A("* `%s` (%d elements)" % (fam, n))
    A("")
    if s["elements_without_ccECP"]:
        A("Element directories that exist but have **no** `ccECP` recipe: %s."
          % ", ".join("`%s`" % e for e in s["elements_without_ccECP"]))
        A("")
    A("## 3. Source representations used")
    A("")
    A("Selection is `--ecp-source %s` / `--basis-source %s`.  In `auto` mode every")
    A("available analytic rendering is parsed independently, all renderings must")
    A("be mutually *consistent* (their difference must not exceed one unit in")
    A("the last written place of the coarser one), and the rendering carrying the")
    A("largest number of significant digits is used.  The rendering actually used")
    A("is recorded per entry in `ecp_name_map.csv` / `basis_name_map.csv`.")
    L[-6] = L[-6] % (s["ecp_source_selection"], s["basis_source_selection"])
    A("")
    A("| kind | rendering | entries |")
    A("|---|---|---|")
    for kind in ("ecp", "basis"):
        for fmt, n in sorted(s["source_formats_used"][kind].items(),
                             key=lambda kv: -kv[1]):
            A("| %s | `.%s` | %d |" % (kind, fmt, n))
    A("")
    A("## 4. GENBAS conversion procedure")
    A("")
    A("1. Parse the source file into one record per *contracted function*,")
    A("   keeping every exponent and coefficient as its **source decimal string**")
    A("   (no float round-trip ever happens).")
    A("2. Group the contracted functions by angular momentum, s first.")
    A("3. For each L build the union of the primitive exponents, sort it in")
    A("   descending exponent order, and place each contracted function in one")
    A("   column of the coefficient matrix required by CFOUR.")
    A("   Two primitives are identified only when their decimals are *exactly*")
    A("   equal; matrix positions a given function does not use are exact zero.")
    A("4. Emit the CFOUR *new-format* entry")
    A("   (`ELEM:name` / comment / blank / NS / L / NC / NE / blank /")
    A("   per-shell exponents / blank / NE x NC coefficient matrix / blank).")
    A("   The numeric records of a GENBAS entry are read list-directed")
    A("   (Fortran `*`). Numerical values are preserved exactly; scientific")
    A("   notation may be rendered in equivalent fixed-point form. The new-format")
    A("   list-directed layout avoids the range and precision limits of the legacy")
    A("   fixed-width GENBAS representation.")
    A("")
    A("Nothing is normalised, recontracted, merged, deleted or refitted.")
    A("")
    A("## 5. ECPDATA conversion procedure")
    A("")
    A("1. Parse every available scalar rendering.  The MOLPRO `!*-so` blocks and")
    A("   the NWChem `so ... end` section are dropped explicitly.")
    A("2. Read `NCORE` and `LMAX` from the source header of each rendering and")
    A("   require them to agree across renderings.")
    A("3. Emit `local channel first`, labelled with the letter of `l = LMAX`,")
    A("   then `s-<lmax>`, `p-<lmax>`, ... for `l = 0 .. LMAX-1`.")
    A("4. Columns are `coefficient   N   alpha`, with `N` copied through")
    A("   unchanged from the source.")
    A("")
    A("### Radial-power convention")
    A("")
    A("```")
    A("U_l(r) = sum_m  c_m  r**(N_m - 2)  exp(-alpha_m r**2)")
    A("```")
    A("")
    A("The integer in the middle ECPDATA column is copied unchanged from the")
    A("GAMESS / Gaussian / NWChem / MOLPRO representations.  The convention is")
    A("checked by comparing the `ECPDATA` file shipped with CFOUR against the")
    A("same literature ECPs taken from the MolSSI Basis Set Exchange, whose")
    A("internal `r_exponents` are exactly those integers:")
    A("")
    A("| CFOUR ECPDATA entry | reference | result |")
    A("|---|---|---|")
    for k, st, d in s["convention_checks"]:
        A("| `%s` | %s | **%s** |" % (k.split("  vs ")[0],
                                      k.split("  vs ")[-1], st))
    A("")
    A("(Note: the CFOUR documentation page writes the formula as `r**N_m`; the")
    A("shipped data files show that the stored integer is the standard one, e.g.")
    A("`BR:ECP-28 HAY & WADT` carries a `-28.0  1  213.6143969` term whose only")
    A("physical reading is `-28 * r**-1 * exp(-213.6 r^2)`.)")
    A("")
    A("### Local-channel convention")
    A("")
    A("The library's `.gamess` rendering writes `LMAX+1` blocks: the local")
    A("potential first (`ul` in the NWChem rendering) followed by `l = 0 ..")
    A("LMAX-1`.  CFOUR expects exactly the same ordering, with the local block")
    A("labelled by the letter of `l = LMAX` and the semi-local blocks labelled")
    A("`<l>-<lmax>` because they represent `U_l - U_LMAX`.  The mapping is")
    A("therefore a direct pass-through.")
    A("")
    A("### NCORE is derived from the data, never from a filename")
    A("")
    A("For every element the coefficients of the `N = 1` (i.e. `r**-1`) terms of")
    A("the local channel are summed and required to equal `Z - NCORE`.  This")
    A("independent check passes for all converted elements (nitrogen, for")
    A("instance, splits `Zeff = 5` into two terms, `3.25 + 1.75`).")
    A("")
    A("## 6. Naming policy")
    A("")
    A("| Pseudopotential Library | CFOUR/MRCC label |")
    A("|---|---|")
    A("| `Fe.ccECP` | `Fe:ccECP` |")
    A("| `Fe.cc-pVQZ` | `Fe:cc-pVQZ%s` |" % s_suffix)
    A("| `Fe.aug-cc-pCV5Z` | `Fe:aug-cc-pCV5Z%s` |" % s_suffix)
    A("")
    A("The element symbol keeps its periodic-table spelling. The source basis")
    A("name is retained as the stem of each orbital-basis label. The **ECP**")
    A("label is exactly the family name, `ccECP`; nothing is renamed to")
    A("`CC-ECP`, `ccecp` or")
    A("`ccECP-cc-pVQZ`, and no prefix is added.")
    A("")
    if s_suffix:
        A("Orbital-basis labels carry the suffix **`%s`**." % s_suffix)
        A("")
        A("This is deliberate and is the one naming transformation in the")
        A("output. Installed MRCC/CFOUR basis libraries may already contain")
        A("entries with conventional names such as `cc-pVTZ` or `aug-cc-pCVTZ`.")
        A("A working-directory `GENBAS` can take precedence over an installed")
        A("library, so an explicit suffix makes the selected ccECP-optimized")
        A("basis unambiguous. MRCC already distinguishes related sets by suffix")
        A("(`-PP`, `-DK`, `-RI`), so `%s` follows that convention." % s_suffix)
        A("")
        A("The source basis name remains intact as the label stem: `cc-pVQZ` ->")
        A("`cc-pVQZ%s`." % s_suffix)
        A("`basis_name_map.csv` records `source_name`, `label_suffix` and")
        A("`cfour_name` for every entry, and `naming_audit.csv` verifies that")
        A("each label is exactly `<element>:<library name><suffix>`.")
        A("Run with `--basis-suffix \"\"` to use the exact source basis names")
        A("after verifying that they do not collide with the target installation.")
        A("")
    else:
        A("Orbital-basis labels use the exact source basis names")
        A("(`--basis-suffix \"\"`). Check the target installation for label")
        A("collisions, because a working-directory `GENBAS` can take precedence")
        A("over installed basis-library entries.")
    A("`naming_audit.csv` records the audit for every entry; failures: **%d**."
      % len(s["naming_failures"]))
    A("")
    A("For CFOUR builds that fold the ZMAT atom string to upper case before")
    A("searching the library, `GENBAS.upper` and `ECPDATA.upper` carry the same")
    A("data with the element symbol upper-cased (`AG:aug-cc-pCV5Z%s`).  Only the"
      % s_suffix)
    A("element symbol changes there; the text after the colon is untouched.")
    A("")
    A("## 7. Contents")
    A("")
    A("* **%d** canonical ccECP elements" % s["elements_found"])
    A("* **%d** `ECPDATA` entries (one `ccECP` per element)" % s["ecpdata_entries"])
    A("* **%d** `GENBAS` entries" % s["genbas_entries"])
    if s.get("ri_entries"):
        A("* **%d** auxiliary (RI) entries" % s["ri_entries"])
    A("* **%d** distinct basis-set names: %s"
      % (len(s["distinct_basis_names"]),
         ", ".join("`%s`" % n for n in s["distinct_basis_names"])))
    A("")
    A("### Elements")
    A("")
    A("| El | Z | NCORE | e- explicit | LMAX | local | non-local | # basis | basis names |")
    A("|---|---|---|---|---|---|---|---|---|")
    for r in inv_rows:
        A("| %s | %d | %s | %s | %s | %s | %s | %d | %s |"
          % (r["element"], r["atomic_number"], r["ncore"],
             r["n_explicit_electrons"], r["lmax"], r["local_ecp_channel"],
             r["nonlocal_ecp_channels"], r["n_basis_sets"],
             r["basis_names"].replace(";", " ")))
    A("")
    A("## 8. Anomalies and notes")
    A("")
    A("* `In` and `Sr` have element directories in the library but **no**")
    A("  `ccECP` recipe, so they are absent from this output.")
    A("* `Ne` has no `aug-` sets in the canonical library; `Li`, `Na` and `Pb`")
    A("  have no `6Z` and no core-valence sets.  Missing members of a series are")
    A("  **not** invented.")
    A("* `H` and `He` have `NCORE = 0`: the ccECP is a shape-consistent")
    A("  potential with all electrons explicit.  Their non-local `s` channel is")
    A("  a single term with coefficient exactly `0.0`, which is preserved.")
    A("* The `B` and `C` `cc-pV{D,T,Q,5}Z` `.gamess` files use lower-case")
    A("  Gaussian-style shell headers (`s 9 1.00`) instead of the usual `S 9`;")
    A("  the parser handles both spellings.")
    A("* `Ce` has a primitive with an exactly zero contraction coefficient in")
    A("  several `f` shells; it is preserved because the exponent is shared with")
    A("  the other `f` functions.")
    A("* The `.gaussian` renderings pad numbers with trailing zeros and, for")
    A("  nitrogen, round two orbital exponents (`9.3345971` -> `9.334597`);")
    A("  the `.gamess`/`.nwchem` renderings carry the extra digit and are used.")
    A("* Molpro **basis** renderings are not used by this converter. Basis `auto`")
    A("  selection compares the supported `.gamess`, `.nwchem` and `.gaussian`")
    A("  renderings. Molpro is supported for ECP cross-checking only.")
    A("* `cc-pV6Z`/`cc-pCV6Z` sets contain `i` functions (L = 6).  Make sure the")
    A("  CFOUR/MRCC build supports that angular momentum.")
    A("")
    A("## 9. MRCC ECP compatibility")
    A("")
    A("MRCC 25.1.1 and 26.1.1 read the `NCORE`/`LMAX` ECPDATA record with")
    A("the fixed Fortran format")
    A("")
    A("```fortran")
    A('read(gbasfile,"(11x,i3,11x,i1)") ncorecp(iatoms),lmax')
    A("```")
    A("")
    A("Accordingly, this converter always writes `NCORE` in columns 12-14 and")
    A("`LMAX` in column 26. `format_ncore_lmax_line()` is the single formatter")
    A("used for every ECPDATA entry, and `selftest_ncore_lmax()` verifies the")
    A("exact positions at startup. The layout is also valid CFOUR input.")
    A("")
    A("```text")
    A("         1         2         3")
    A("123456789012345678901234567890")
    for elem, ncore, lmax in NCORE_LMAX_TEST_CASES:
        A(format_ncore_lmax_line(ncore, lmax) + "   <- %s" % elem)
    A("```")
    A("")
    A("### `NCORE = 0` in MRCC")
    A("")
    A("Executable validation shows that MRCC 25.1.1 and 26.1.1 do not correctly")
    A("handle ECP entries with `NCORE = 0`. This affects the canonical ccECPs")
    A("of **H** and **He**, which keep all electrons explicit. The converter")
    A("still writes the correct `NCORE = 0` data and applies no numerical")
    A("workaround. Validate H and He independently before production use with")
    A("those MRCC releases. Canonical ccECP entries with `NCORE > 0` pass the")
    A("fixed-column ECP parsing checks.")
    A("")
    A("## 10. Validation")
    A("")
    counts = collections.Counter(r["status"] for r in val_rows)
    A("`validation_report.csv` holds one row per check (%d rows)." % len(val_rows))
    A("")
    A("| status | rows |")
    A("|---|---|")
    for k, v in sorted(counts.items()):
        A("| %s | %d |" % (k, v))
    A("")
    A("Checks performed for **every basis set**: source parse of all available")
    A("renderings and mutual consistency; independent parse by")
    A("`basis_set_exchange`; contracted-function count; primitive count;")
    A("angular-momentum grouping and ordering; exponent ordering and")
    A("near-degeneracy; highest-L retention; and a full")
    A("write -> independent re-read -> exact-decimal comparison round trip.")
    A("Plus cross-basis checks that every `aug-X` retains all primitives of `X`")
    A("and that every `cc-pCVnZ` retains all primitives of `cc-pVnZ`.")
    A("")
    A("Checks performed for **every ECP**: channel count vs `LMAX`; `NCORE`")
    A("agreement across all renderings; `NCORE` re-derived from the local")
    A("channel's `r**-1` coefficients; radial-power range; absence of `r**-2`")
    A("terms; spin-orbit exclusion; cross-rendering numerical agreement; and a")
    A("write -> independent re-read -> exact-decimal comparison round trip.")
    A("")
    A("Cross-basis hierarchy results: %r" % s["hierarchy"])
    A("")
    if s["failed_ecp"] or s["failed_basis"]:
        A("### Failures")
        A("")
        for e, m in s["failed_ecp"]:
            A("* ECP `%s`: %s" % (e, m))
        for t in s["failed_basis"]:
            A("* basis `%s/%s`: %s" % t)
    else:
        A("**No failures**: every canonical ccECP element and every canonical")
        A("ccECP basis set converted and validated.")
    A("")
    A("## 11. Files")
    A("")
    A("All paths below are relative to the output directory chosen with")
    A("`--out`. If omitted, the script chooses a safe `ccECP_cfour_output`")
    A("directory outside the source repository.")
    A("")
    A("```")
    A("GENBAS                 CFOUR basis-set library (all %d entries)"
      % s["genbas_entries"])
    A("ECPDATA                CFOUR ECP library (all %d entries)"
      % s["ecpdata_entries"])
    A("GENBAS.mrcc            basis data with ECP entries appended for MRCC;")
    A("                       copy/rename this file to GENBAS in an MRCC job")
    if s.get("ri_entries"):
        A("GENBAS.RI              auxiliary (RI) fitting sets only (%d entries)"
          % s["ri_entries"])
        A("GENBAS.withRI          orbital + auxiliary sets, for CFOUR")
        A("GENBAS.mrcc.withRI     orbital + auxiliary + ECP in one file; copy")
        A("                       this to GENBAS for an MRCC job that uses")
        A("                       dfbasis_cor")
        A("ri_name_map.csv        element / source file / library name / CFOUR")
        A("                       name for every auxiliary set")
    A("ECPDATA.maxprec        ECPs taken from the most precise rendering of")
    A("                       each element rather than from one uniform one.")
    A("                       For the 3d transition metals the .molpro ECP")
    A("                       rendering carries up to seven more significant")
    A("                       digits than .gamess, and for Ag, Au, Bi, I, Ir,")
    A("                       Pd and Te the native <El>.ccECP text does; this")
    A("                       file uses whichever is most precise, after")
    A("                       verifying all renderings agree.")
    A("*.upper                the same four files with upper-cased element")
    A("                       symbols, for CFOUR builds that need that")
    A("inventory.csv/.json    what was found in the library")
    A("basis_name_map.csv     element / source file / library name / CFOUR name")
    A("ecp_name_map.csv       element / source file / library name / CFOUR name")
    A("naming_audit.csv       per-entry naming audit")
    A("validation_report.csv  one row per validation check")
    A("anomalies.csv          every non-PASS check plus source-data and")
    A("                       target-program compatibility notes")
    A("summary.json           machine-readable run summary")
    A("conversion.log         full run log")
    A("by_element/<El>/       per-element GENBAS, ECPDATA and summary.json,")
    A("                       plus one file per logical object using the")
    A("                       source object name with a `.cfour` suffix:")
    A("                         Fe.ccECP.cfour      (one ECPDATA entry)")
    A("                         Fe.cc-pVQZ.cfour    (one GENBAS entry)")
    A("                         Fe.aug-cc-pCV5Z.cfour  ...")
    A("                       mirroring Fe.ccECP.gamess / Fe.cc-pVQZ.gamess.")
    A("                       MRCC uses the same two formats, so these files")
    A("                       serve both programs.")
    A("```")
    A("")
    A("## 12. Example CFOUR input (`ZMAT`)")
    A("")
    A("Copy `GENBAS` and `ECPDATA` into the CFOUR working directory.")
    A("")
    A("```")
    A("CO, CCSD(T)/cc-pVTZ with canonical ccECPs")
    A("C")
    A("O 1 R")
    A("")
    A("R=1.1280")
    A("")
    A("*CFOUR(CALC=CCSD(T),BASIS=SPECIAL,ECP=ON,SPHERICAL=ON")
    A("REF=RHF,SCF_CONV=9,CC_CONV=9,ABCDTYPE=AOBASIS)")
    A("")
    A("C:cc-pVTZ%s" % s_suffix)
    A("O:cc-pVTZ%s" % s_suffix)
    A("")
    A("C:ccECP")
    A("O:ccECP")
    A("")
    A("```")
    A("")
    A("The basis labels carry the `%s` suffix; the ECP labels do not." % s_suffix
      if s_suffix else "")
    A("")
    A("`BASIS=SPECIAL` makes CFOUR read the per-atom basis labels from the first")
    A("block and, because `ECP=ON`, the per-atom ECP labels from the second.")
    A("The labels are spelled exactly as in `GENBAS`/`ECPDATA`, i.e. with the")
    A("element symbol in periodic-table case (`C:`/`O:`). If your CFOUR build")
    A("folds atom labels to upper case, use `GENBAS.upper` and `ECPDATA.upper`")
    A("with the corresponding upper-case element labels.")
    A("Use `<Element>:NONE` in the ECP block for an atom that should carry no ECP.")
    A("")
    A("## 13. Example MRCC input (`MINP`)")
    A("")
    A("MRCC reads user basis sets *and* user ECPs from a single `GENBAS` file in")
    A("the working directory, so use the generated `GENBAS.mrcc`:")
    A("")
    A("```")
    A("cp GENBAS.mrcc ./GENBAS")
    A("```")
    A("")
    A("```")
    A("basis=special")
    A("cc-pVTZ%s" % s_suffix)
    A("cc-pVTZ%s" % s_suffix)
    A("")
    A("ecp=special")
    A("ccECP")
    A("ccECP")
    A("")
    A("calc=CCSD(T)")
    A("mem=4GB")
    A("scftype=RHF")
    A("gauss=spher")
    A("unit=angs")
    A("geom=xyz")
    A("2")
    A("")
    A("C    0.0000   0.0000   0.0000")
    A("O    0.0000   0.0000   1.1280")
    A("```")
    A("")
    A("The `basis=special` / `ecp=special` lists must follow the atom order of")
    A("the geometry, and the ECP labels must match the `GENBAS` labels exactly.")
    A("")
    A("For canonical ccECP H and He, read the `NCORE = 0` compatibility note")
    A("in section 9 before production use with MRCC 25.1.1 or 26.1.1.")
    A("")
    if s.get("ri_entries"):
        A("## 13b. Auxiliary (RI) basis sets")
        A("")
        A("The library ships an ORCA AutoAux `/C` auxiliary set beside most")
        A("orbital basis sets, as `<El>.<basis>.AutoAuxC.orca`.  These are")
        A("density-fitting sets for the correlation step, and **%d** of them are"
          % s["ri_entries"])
        A("converted here: %d taken from the library and %d supplied through"
          % (s.get("ri_from_library", 0), s.get("ri_from_extra", 0)))
        A("`--extra-aux` for orbital sets the library does not cover.")
        A("")
        A("They travel through the same parser body, the same block builder and")
        A("the same renderer as the orbital sets, so there is no second")
        A("numerical path.  Only the wrapper differs in the source: two lines")
        A("at the top (`# ...` and `NewAuxCGTO <El>`) and a terminating `end;`.")
        A("")
        A("### Naming")
        A("")
        A("The label is the orbital basis name, then `%s`, then the basis-label"
          % RI_NAME_SUFFIX)
        A("suffix:")
        A("")
        A("```")
        A("Ag:aug-cc-pCV5Z-RI%s" % s_suffix)
        A("```")
        A("")
        A("`-RI` is MRCC's own name for a correlation-fitting set, as against")
        A("`-RI-JK` for an SCF-fitting one; the AutoAux `/C` sets are")
        A("correlation-fitting, so `-RI` is the matching name.  The suffix after")
        A("it is still required: MRCC already ships entries such as")
        A("`Ag:aug-cc-pV5Z-RI`, so the bare name is already taken.")
        c = s.get("collisions")
        if c:
            A("")
            A("Measured against the installed library given to")
            A("`--host-basis-library` (%d labels): the auxiliary labels written"
              % c["host_labels"])
            A("here collide **%d** times, and would collide **%d** times if"
              % (c["ri"], c["ri_bare"]))
            A("`%s` were dropped." % s_suffix)
        A("")
        A("### What the sets look like")
        A("")
        A("Every AutoAux function is a single uncontracted primitive with")
        A("coefficient exactly 1, so each angular momentum becomes an identity")
        A("coefficient matrix.  That is precisely how MRCC stores its own")
        A("`*-RI` entries, so the layout needs no special handling.  Angular")
        A("momentum reaches `k` (l = 7) in the heavier sets, which is within")
        A("what MRCC accepts -- its own `Ag:cc-pV5Z-PP-RI` reaches `l`.")
        A("")
        A("### Use in MRCC")
        A("")
        A("```")
        A("cp GENBAS.mrcc.withRI ./GENBAS")
        A("```")
        A("")
        A("```")
        A("basis=atomtype")
        A("Fe:cc-pVTZ%s" % s_suffix)
        A("dfbasis_scf=none")
        A("dfbasis_cor=cc-pVTZ-RI%s" % s_suffix)
        A("ecp=atomtype")
        A("Fe:ccECP")
        A("calc=MP2")
        A("```")
        A("")
        A("Note that `dfbasis_cor` takes the label without the element prefix.")
        A("")
        A("### Orbital sets with no auxiliary partner")
        A("")
        A("Missing auxiliary sets are recorded, never reconstructed.  Any")
        A("orbital set without one is listed in `anomalies.csv` with")
        A("`check=aux_present`; there are **%d** such entries in this run."
          % s.get("ri_missing", 0))
        A("")
    A("## 14. Troubleshooting and reproducibility")
    A("")
    A("* **Library not found:** pass `--library /path/to/pseudopotentiallibrary`")
    A("  or set `PSEUDOPOTENTIAL_LIBRARY`. The option may also point directly")
    A("  to the `recipes/` directory.")
    A("* **`basis_set_exchange` unavailable:** install it with `python -m pip`")
    A("  `install basis_set_exchange`, or use `--no-bse` to skip only the")
    A("  independent BSE checks.")
    A("* **Basis-label collision:** keep the default `%s` suffix, or use" % s_suffix)
    A("  `--basis-suffix \"\"` only after checking the target basis library.")
    A("* **CFOUR cannot find a mixed-case label:** try the corresponding")
    A("  `*.upper` files. Only the element symbol changes.")
    A("* **Very high angular momentum:** some 6Z sets contain `i` functions")
    A("  (`L = 6`). Confirm that the target CFOUR/MRCC build supports them.")
    A("* **Reproducibility:** retain `summary.json`, `conversion.log`, the mapping")
    A("  CSV files, and the Pseudopotential Library commit hash.")
    A("")
    A("## 15. Citation and source data")
    A("")
    A("The converter does not redefine or refit ccECP data. Scientific credit")
    A("belongs to the original ccECP/Pseudopotential Library sources. Each")
    A("recipe's `author.txt` is carried into generated comments and should be")
    A("consulted for the appropriate literature citation. Source-data observations")
    A("and non-fatal consistency notes are recorded in `anomalies.csv`; numerical")
    A("source values are preserved rather than silently modified.")
    A("")
    with open(os.path.join(outdir, "README.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")


if __name__ == "__main__":
    sys.exit(main())
