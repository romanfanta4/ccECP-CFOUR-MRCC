#!/usr/bin/env python3
"""Convert GAMESS-US ECP and orbital-basis files below a recipes tree to ORCA.

Run from the pseudopotential library's ``recipes`` directory:

    python3 convert_gamess_ecp_to_orca.py

The script recursively inspects ``*.gamess`` files.  Files whose first
significant line has the form ``Element-name GEN ncore lmax`` become ORCA
``NewECP`` blocks.  Orbital basis files are wrapped in the ORCA-readable
GAMESS-US form ``$DATA / Element / shells / $END``.  For example:

* ``Fe/ccECP/Fe.ccECP.gamess`` becomes ``Fe/ccECP/Fe.ccECP.orca``.
* ``Fe/ccECP/Fe.aug-cc-pVDZ.gamess`` becomes
  ``Fe/ccECP/Fe.aug-cc-pVDZ.orca``.

Existing outputs are skipped when they are identical to the expected
conversion.  A differing existing output is preserved unless ``--overwrite``
is supplied.  Files are written atomically, so an interrupted run can be
restarted safely.

Two unrelated file formats share the ``.orca`` extension, because the two
kinds of input do:

* An ECP conversion is a bare ``NewECP ... end`` block.  ORCA reads it when the
  block sits inside the input's ``%basis ... end`` section.  Checked against
  ORCA 6.1.0: inlining the block works, while ``%include`` of the file and
  ``ReadFragECP`` of the file are both rejected, so the block has to be pasted
  in rather than referenced.
* A basis conversion is a GAMESS-US ``$DATA ... $END`` file, which ORCA reads
  through ``%basis  GTOName = "<file>"  end``.

Numerical fields travel as the strings found in the source and are never
reparsed into floats on the way out, so no digit of the library data is lost.
ORCA's own ``ECP BASIS IN INPUT FORMAT`` echo reproduces the source exponents,
coefficients and radial powers exactly.

Both GAMESS ECP spellings of a group header are read: the bare primitive count
used by the ccECP library, and the labelled form written by GAMESS-US and the
Basis Set Exchange (``3     ----- d-ul potential -----``).

Shell headers are emitted in the one spelling ORCA accepts, ``S 9``.  Boron and
carbon ship their ``cc-pV{D,T,Q,5}Z`` sets in the Gaussian spelling ``s 9 1.00``
instead, which makes ORCA hang while reading the basis; see
normalize_shell_header.

Every parsed ECP is checked against ``Z - n_core`` through the r**-1 terms of
its local channel.  For a potential whose header says ``ccECP`` a mismatch is
fatal; for any other family it is only reported, because the invariant follows
from the ccECP form rather than from the file format.  See coulomb_tail_error.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import os
from pathlib import Path
import re
import stat
import sys
import tempfile


# Conventional angular-momentum labels omit J.
ANGULAR_MOMENTA = "spdfghikl"

ELEMENTS = (
    "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
    "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr",
    "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
    "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po",
    "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs",
    "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
)
ATOMIC_NUMBERS = {symbol.upper(): z for z, symbol in enumerate(ELEMENTS) if symbol}

# The header names the element and, after a dash, the potential itself
# ("Fe-ccECP GEN 10 2").  That name is captured because it selects how strictly
# the Coulomb tail is checked; see coulomb_tail_error.
ECP_HEADER = re.compile(
    r"([A-Za-z]{1,2})(?:-([^\s]+))?\s+GEN\s+(\d+)\s+(\d+)",
    flags=re.IGNORECASE,
)

# A shell header names the shell type and its primitive count, optionally
# followed by a scale factor in the Gaussian spelling ("s 9 1.00").
SHELL_HEADER = re.compile(r"(\s*)([A-Za-z]{1,2})(\s+)(\d+)(?:\s+(\S+))?\s*")

# A group header is a primitive count, optionally followed by descriptive text.
# The ccECP library writes the count alone; GAMESS-US and the Basis Set
# Exchange append a label, as in "3     ----- d-ul potential -----".  Both are
# accepted.  A label that is itself entirely numeric is refused, so a primitive
# line can never be misread as a group header.
GROUP_HEADER = re.compile(r"(\d+)(?:\s+(\S.*))?")


@dataclass(frozen=True)
class ECPPrimitive:
    """One radial term of a channel: c * r**(power - 2) * exp(-alpha * r**2).

    ``coefficient`` and ``exponent`` stay as the source strings, so the output
    carries the library digits unchanged.  Only ``power`` is interpreted, and
    it is written to ORCA unaltered because both codes print the same integer
    for the same radial factor.
    """

    coefficient: str
    power: int
    exponent: str


@dataclass(frozen=True)
class GamessECP:
    element: str
    #: potential name as spelled in the header, e.g. ``ccECP``; ``""`` if the
    #: header carried no ``-<name>`` suffix.
    name: str
    n_core: int
    local_l: int
    # GAMESS order: local channel first, then s, p, ..., local_l - 1.  ORCA
    # wants the local channel last; format_orca_ecp does that permutation.
    # The semi-local groups are differences U_l - U_local in both codes, so no
    # arithmetic is needed -- only a reordering.
    groups: tuple[tuple[ECPPrimitive, ...], ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Recursively convert GAMESS-US GEN ECP files to ORCA NewECP blocks "
            "and wrap GAMESS orbital bases for direct ORCA use."
        )
    )
    parser.add_argument(
        "folder",
        nargs="?",
        type=Path,
        default=Path("."),
        help="recipes directory or one GAMESS ECP/basis file (default: current directory)",
    )
    parser.add_argument(
        "--pattern",
        default="*.gamess",
        help="recursive input glob (default: *.gamess)",
    )
    parser.add_argument(
        "--extension",
        default=".orca",
        help="output extension, including the leading dot (default: .orca)",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="replace existing outputs that differ from the expected conversion",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="show conversions without writing files",
    )
    return parser.parse_args()


def clean_ecp_lines(text: str) -> list[str]:
    """Remove blank lines, comments, and optional GAMESS group delimiters."""
    lines: list[str] = []
    for raw in text.splitlines():
        line = raw.split("!", 1)[0].split("#", 1)[0].strip()
        if not line or line.upper() in {"$ECP", "$END"}:
            continue
        lines.append(line)
    return lines


def normalize_element(value: str) -> str:
    key = value.upper()
    if key not in ATOMIC_NUMBERS:
        raise ValueError(f"unknown element symbol {value!r}")
    return ELEMENTS[ATOMIC_NUMBERS[key]]


def parse_number(token: str) -> float:
    return float(token.replace("D", "E").replace("d", "e"))


def is_number(token: str) -> bool:
    try:
        parse_number(token)
    except ValueError:
        return False
    return True


def parse_group_header(line: str) -> int | None:
    """Return the primitive count on a GAMESS ECP group header, else ``None``.

    Trailing descriptive text is allowed, because GAMESS-US writes it, but only
    when it is not itself a list of numbers.  Without that restriction a
    primitive line such as ``16 1 23.2`` would be indistinguishable from a
    count of 16 followed by a label.
    """
    match = GROUP_HEADER.fullmatch(line)
    if match is None:
        return None
    label = match.group(2)
    if label is not None and all(is_number(field) for field in label.split()):
        return None
    return int(match.group(1))


def looks_like_gen_ecp(text: str) -> bool:
    lines = clean_ecp_lines(text)
    return bool(lines and ECP_HEADER.fullmatch(lines[0]))


def basis_content_lines(text: str) -> list[str]:
    """Return significant basis lines without optional GAMESS wrappers/header."""
    lines: list[str] = []
    for raw in text.splitlines():
        stripped = raw.strip()
        if stripped.upper() in {"$DATA", "$END"}:
            continue
        significant = raw.split("!", 1)[0].split("#", 1)[0].strip()
        if significant:
            lines.append(significant)

    # A complete GAMESS-US basis may already contain an element header such as
    # "Fe" or "IRON".  Shell headers always also contain a primitive count, so
    # a leading alphabetic-only line is unambiguously an element header here.
    if lines and re.fullmatch(r"[A-Za-z]+", lines[0]):
        lines.pop(0)
    return lines


def looks_like_orbital_basis(text: str) -> bool:
    """Detect a GAMESS shell block from its leading ``<type> <nprim>`` line.

    ``SP`` and ``L`` are GAMESS spellings of a fused s+p shell.  They are
    recognized so such a file is not silently ignored, but they are never
    interpreted, because shell data is copied through untouched.
    """
    lines = basis_content_lines(text)
    if not lines:
        return False
    return re.fullmatch(
        r"(?:S|P|D|F|G|H|I|K|L|SP)\s+\d+(?:\s+.*)?",
        lines[0],
        flags=re.IGNORECASE,
    ) is not None


def parse_gamess_ecp(text: str, source: str) -> GamessECP:
    """Parse ``Element-name GEN ncore lmax`` GAMESS-US ECP syntax."""
    lines = clean_ecp_lines(text)
    if not lines:
        raise ValueError(f"{source}: empty file")

    header = ECP_HEADER.fullmatch(lines[0])
    if header is None:
        raise ValueError(
            f"{source}: expected 'Element-name GEN ncore lmax'; got {lines[0]!r}"
        )

    element = normalize_element(header.group(1))
    name = header.group(2) or ""
    n_core = int(header.group(3))
    local_l = int(header.group(4))

    # n_core == 0 is legitimate: the ccECP potentials for H and He keep every
    # electron explicit.  ORCA 6.1.0 handles that correctly, unlike some other
    # programs, so it is allowed through here.
    if n_core < 0 or n_core >= ATOMIC_NUMBERS[element.upper()]:
        raise ValueError(f"{source}: invalid core-electron count {n_core} for {element}")
    if local_l >= len(ANGULAR_MOMENTA):
        raise ValueError(f"{source}: unsupported local angular momentum {local_l}")

    cursor = 1
    groups: list[tuple[ECPPrimitive, ...]] = []
    expected_groups = local_l + 1

    for group_index in range(expected_groups):
        count = parse_group_header(lines[cursor]) if cursor < len(lines) else None
        if count is None:
            raise ValueError(
                f"{source}: missing primitive count for group "
                f"{group_index + 1}/{expected_groups}"
            )
        cursor += 1
        if count < 1:
            raise ValueError(f"{source}: ECP groups cannot be empty")

        primitives: list[ECPPrimitive] = []
        for _ in range(count):
            if cursor >= len(lines):
                raise ValueError(f"{source}: ECP group {group_index + 1} is truncated")

            primitive_line = lines[cursor]
            cursor += 1
            fields = primitive_line.split()
            if len(fields) != 3:
                raise ValueError(
                    f"{source}: expected 'coefficient power exponent'; "
                    f"got {primitive_line!r}"
                )

            coefficient, power_token, exponent = fields
            try:
                # The two parse_number results are discarded on purpose: they
                # only prove the fields are numeric.  Keeping the original
                # strings is what avoids any rounding of the source data.
                parse_number(coefficient)
                power = int(power_token)
                parse_number(exponent)
            except ValueError as exc:
                raise ValueError(f"{source}: invalid primitive {primitive_line!r}") from exc

            primitives.append(ECPPrimitive(coefficient, power, exponent))

        groups.append(tuple(primitives))

    if cursor != len(lines):
        raise ValueError(f"{source}: unexpected trailing content {lines[cursor]!r}")

    ecp = GamessECP(element, name, n_core, local_l, tuple(groups))
    message = coulomb_tail_error(ecp)
    if message is not None and is_ccecp(ecp):
        raise ValueError(f"{source}: {message}")
    return ecp


def is_ccecp(ecp: GamessECP) -> bool:
    """Does the header claim this potential is a ccECP?"""
    return "ccecp" in ecp.name.lower()


def coulomb_tail_error(ecp: GamessECP) -> str | None:
    """Check the local channel against ``Z - n_core``, or return ``None``.

    A ccECP writes its local channel as

        -Zeff/r + (Zeff/r) * exp(-alpha * r**2) + ...

    and GAMESS supplies the bare ``-Zeff/r`` implicitly, so the file's own
    r**-1 terms -- the ones with ``power == 1`` -- must add up to ``+Zeff``.
    This is the strongest check available on a parsed ECP: it fails if the
    group order, the column assignment or the core count was misread.

    They have to be *summed*, not taken one at a time: nitrogen's ccECP splits
    the tail over two terms, 3.25 + 1.75 = 5.

    The invariant is specific to how ccECP regularizes the tail, and is not a
    property of the GAMESS ECP format.  LANL2DZ/Cu (-10 against Zeff = 19) and
    CRENBL/Fe (-8.99 against Zeff = 16) are both legitimate and both violate
    it, so this is only fatal for a potential whose header says ``ccECP``; for
    anything else main() downgrades it to a warning.
    """
    z_effective = ATOMIC_NUMBERS[ecp.element.upper()] - ecp.n_core
    tail = sum(parse_number(primitive.coefficient)
               for primitive in ecp.groups[0] if primitive.power == 1)
    if abs(tail - z_effective) <= 1e-6 * max(1.0, abs(z_effective)):
        return None
    return (
        f"local-channel r**-1 coefficients sum to {tail:.10g}, "
        f"expected Z - n_core = {z_effective}"
    )


def orca_number(token: str) -> str:
    """Normalize Fortran D exponents to notation accepted by ORCA."""
    return token.replace("D", "E").replace("d", "e")


def format_orca_ecp(ecp: GamessECP, source_name: str) -> str:
    """Reorder GAMESS groups and columns into an ORCA ``NewECP`` block.

    Three conventions have to agree for a pure reordering to be correct.  All
    three were confirmed by rendering the same ECP in both formats with the
    MolSSI Basis Set Exchange (LANL2DZ/Cu and CRENBL/Fe):

    * the trailing integer carries the same meaning in both codes, so ``power``
      is copied across rather than shifted by two;
    * ORCA orders the columns as index, exponent, coefficient, power, which is
      neither the GAMESS column order nor a simple reversal of it;
    * ORCA lists the semi-local groups first and the local group last, named by
      ``lmax``, which is the opposite of the GAMESS placement.

    Nothing here recomputes a coefficient, so a failure of any of these
    assumptions would show up as reordered output, not as altered numbers.
    """
    local_label = ANGULAR_MOMENTA[ecp.local_l]
    lines = [
        f"# Converted from {source_name}",
        "# GAMESS: coefficient power exponent",
        "# ORCA:   index exponent coefficient power",
        f"NewECP {ecp.element}",
        f"  N_core {ecp.n_core}",
        f"  lmax {local_label}",
        "",
    ]

    # GAMESS group order: local, s, p, ...
    # ORCA group order:    s, p, ..., local
    ordered_groups = [
        (ANGULAR_MOMENTA[l_value], ecp.groups[l_value + 1])
        for l_value in range(ecp.local_l)
    ]
    ordered_groups.append((local_label, ecp.groups[0]))

    for label, primitives in ordered_groups:
        lines.append(f"  {label} {len(primitives)}")
        for index, primitive in enumerate(primitives, start=1):
            lines.append(
                f"    {index} {orca_number(primitive.exponent)} "
                f"{orca_number(primitive.coefficient)} {primitive.power}"
            )
        lines.append("")

    lines.append("end")
    return "\n".join(lines) + "\n"


def infer_basis_element(source: Path) -> str:
    """Infer the element from ``Element.*.gamess`` or the recipes path."""
    prefix = re.match(r"([A-Za-z]{1,2})(?=[._-])", source.name)
    if prefix is not None:
        try:
            return normalize_element(prefix.group(1))
        except ValueError:
            pass

    for parent in source.parents:
        try:
            return normalize_element(parent.name)
        except ValueError:
            continue
    raise ValueError(
        f"{source}: cannot infer element from filename or parent directories"
    )


def normalize_shell_header(line: str) -> str:
    """Rewrite a shell header into the only spelling ORCA's GAMESS reader takes.

    ORCA needs an upper-case shell letter and exactly two fields.  Anything
    else makes it hang while reading the basis, with no error message: tested
    with ORCA 6.1.0, ``s 9 1.00``, ``S 9 1.00`` and ``s 9`` all fail, and only
    ``S 9`` works.

    This matters because ``B`` and ``C`` ship their ``cc-pV{D,T,Q,5}Z`` sets in
    the Gaussian spelling ``s 9 1.00`` rather than the GAMESS ``S 9`` used by
    every other element -- 8 of the library's 901 basis files.

    Only the header spelling changes.  A scale factor multiplies the exponents
    of its shell by its square, so one that is not 1 is refused rather than
    discarded; every occurrence in the library is exactly ``1.00``.
    """
    match = SHELL_HEADER.fullmatch(line)
    if match is None:
        return line
    indent, letter, gap, count, scale = match.groups()
    if scale is None and not letter.islower():
        return line          # already GAMESS form; do not even touch spacing
    if scale is not None and parse_number(scale) != 1.0:
        raise ValueError(
            f"shell header {line.strip()!r} carries scale factor {scale}, "
            f"which would rescale the exponents; refusing to drop it"
        )
    return f"{indent}{letter.upper()}{gap}{count}"


def format_orca_basis(text: str, element: str) -> str:
    """Normalize a shell-only or wrapped GAMESS basis for ORCA's GTO reader.

    Primitive lines are copied verbatim; exponents and coefficients are never
    reformatted, renumbered or re-sorted, which is why a source-data anomaly
    stays visible in the output instead of being smoothed over.  Only the
    optional ``$DATA``/``$END`` wrapper, the element header and the shell-header
    spelling are rewritten -- see normalize_shell_header for why the last of
    those is unavoidable.

    ``basis_content_lines`` strips the same element header for detection
    purposes; the two need to stay in step.
    """
    original_lines = text.splitlines()
    body = [
        line
        for line in original_lines
        if line.strip().upper() not in {"$DATA", "$END"}
    ]

    # Remove an existing symbol/full-name element header.  Ignore comments and
    # blank lines while locating it, but retain those lines in the body.
    for index, line in enumerate(body):
        significant = line.split("!", 1)[0].split("#", 1)[0].strip()
        if not significant:
            continue
        if re.fullmatch(r"[A-Za-z]+", significant):
            del body[index]
        break

    while body and not body[0].strip():
        body.pop(0)
    while body and not body[-1].strip():
        body.pop()

    body = [normalize_shell_header(line) for line in body]
    return f"$DATA\n{element}\n\n" + "\n".join(body) + "\n$END\n"


def output_mode(path: Path) -> int:
    """Mode a plain ``open(path, "w")`` would have left on ``path``.

    ``mkstemp`` hardcodes 0600 and ``os.replace`` preserves it, so an output
    written without this would be owner-only, unlike the world-readable
    ``.gamess`` sources beside it.  A file that already exists keeps the mode
    it already had, so ``--overwrite`` never changes permissions.
    """
    try:
        return stat.S_IMODE(path.stat().st_mode)
    except OSError:
        mask = os.umask(0)
        os.umask(mask)
        return 0o666 & ~mask


def atomic_write(path: Path, text: str) -> None:
    """Replace ``path`` only after a complete temporary file has been written."""
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        text=True,
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(temporary, output_mode(path))
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def discover_inputs(root: Path, pattern: str) -> list[Path]:
    if root.is_file():
        return [root]
    return sorted(path for path in root.rglob(pattern) if path.is_file())


def display_path(path: Path, root: Path) -> str:
    if root.is_dir():
        try:
            return str(path.relative_to(root))
        except ValueError:
            pass
    return str(path)


def main() -> int:
    args = parse_args()
    root = args.folder.expanduser().resolve()

    if not root.exists():
        print(f"ERROR: path does not exist: {root}", file=sys.stderr)
        return 2
    if not args.extension.startswith(".") or args.extension == ".":
        print("ERROR: --extension must start with a dot", file=sys.stderr)
        return 2
    if args.extension.lower() == ".gamess":
        print("ERROR: output extension cannot be .gamess", file=sys.stderr)
        return 2

    inputs = discover_inputs(root, args.pattern)
    ecp_converted = 0
    basis_converted = 0
    skipped_complete = 0
    ignored_unrecognized = 0
    warned = 0
    failed = 0

    for source in inputs:
        source_label = display_path(source, root)
        try:
            text = source.read_text(encoding="utf-8")
        except OSError as exc:
            failed += 1
            print(f"FAIL  {source_label}: {exc}", file=sys.stderr)
            continue

        try:
            output = source.with_suffix(args.extension)
            if looks_like_gen_ecp(text):
                kind = "ECP"
                ecp = parse_gamess_ecp(text, source_label)
                # A ccECP that broke the invariant already raised inside the
                # parser.  Reaching here with a message means some other
                # family, where the tail is genuinely written differently, so
                # it is reported without failing the conversion.
                message = coulomb_tail_error(ecp)
                if message is not None:
                    warned += 1
                    print(f"WARN  {source_label}: {message}", file=sys.stderr)
                expected = format_orca_ecp(ecp, source.name)
            elif looks_like_orbital_basis(text):
                kind = "BASIS"
                element = infer_basis_element(source)
                expected = format_orca_basis(text, element)
            else:
                ignored_unrecognized += 1
                print(f"IGNORE {source_label}: not a recognized GEN ECP or orbital basis")
                continue

            if output.exists():
                current = output.read_text(encoding="utf-8")
                if current == expected:
                    skipped_complete += 1
                    print(
                        f"SKIP  {kind:<5} {source_label}: "
                        f"{output.name} already complete"
                    )
                    continue
                if not args.overwrite:
                    raise FileExistsError(
                        f"{output.name} exists but differs; use --overwrite to replace it"
                    )

            if args.dry_run:
                print(f"READY {kind:<5} {source_label} -> {output.name}")
            else:
                atomic_write(output, expected)
                print(f"DONE  {kind:<5} {source_label} -> {output.name}")

            if kind == "ECP":
                ecp_converted += 1
            else:
                basis_converted += 1

        except (OSError, ValueError) as exc:
            failed += 1
            print(f"FAIL  {source_label}: {exc}", file=sys.stderr)

    print(
        f"Inspected {len(inputs)} GAMESS file(s): converted {ecp_converted} ECP, "
        f"converted {basis_converted} basis, already complete "
        f"{skipped_complete}, unrecognized {ignored_unrecognized}, "
        f"warned {warned}, failed {failed}."
    )
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
