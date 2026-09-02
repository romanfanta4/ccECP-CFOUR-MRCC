#!/usr/bin/env python3
"""Check converted ccECP ORCA files against the library's own reference energies.

Run from this directory:

    python3 test_orca_conversion.py             # every case
    python3 test_orca_conversion.py Fe Cu       # only these
    python3 test_orca_conversion.py --rerun     # ignore finished cases

A case that already produced an energy is reused rather than repeated, so an
interrupted run can simply be started again.

For each element it builds an ORCA input that uses only converted files --
``<El>.cc-pV5Z.orca`` for the orbital basis and ``<El>.ccECP.orca`` for the
potential -- runs a UKS/PBE atomic calculation, and compares the total energy
with the ``Pyscf/PBE`` value in that element's ``energies.txt``.  Those
reference values are UKS/PBE/cc-pV5Z (see NOTE 2 there), so the same method and
basis are used here and the numbers are directly comparable.

Agreement at the 1e-4 Ha level is the expected outcome: the reference is quoted
to five decimals and the two programs use different DFT integration grids.  A
misread ECP would show up as a whole-Hartree error, and a misread contraction
as several mHa, so this is a sharp test of the conversion despite the loose
tolerance.

Nothing here writes to the recipes tree; all scratch stays in this directory.
"""

import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RECIPES = os.path.dirname(HERE)
ORCA = os.environ.get(
    "ORCA_EXE", os.path.expanduser("~/codes/orca_6_1_0_avx2/orca"))
BASIS = "cc-pV5Z"

#: element -> spin multiplicity of the neutral ground state
CASES = [("C", 3), ("O", 3), ("Si", 3), ("Fe", 5),
         ("Cu", 2), ("Br", 2), ("I", 2), ("Au", 2)]

#: NORI keeps the four-centre integrals exact, so no auxiliary basis has to be
#: invented for an ECP basis that ships without one.  DEFGRID3 and
#: VeryTightSCF push grid and convergence error below the comparison level.
TEMPLATE = """! UKS PBE NORI DEFGRID3 VeryTightSCF PAL8
%maxcore 3000
%scf
  MaxIter 500
end
%basis
  GTOName = "{basis_file}"
{ecp}
end
* xyz 0 {mult}
  {element}  0.0  0.0  0.0
*
"""


def reference(element):
    """The Pyscf/PBE value from the element's energies.txt, or None."""
    path = os.path.join(RECIPES, element, "ccECP", "energies.txt")
    if not os.path.isfile(path):
        return None
    with open(path) as handle:
        for line in handle:
            match = re.match(r"\s*Pyscf/PBE\S*\s+(-?\d+\.\d+)", line)
            if match:
                return float(match.group(1))
    return None


def harvest(log):
    """Pull energy, core count and electron count out of an ORCA output."""
    if not os.path.isfile(log):
        return None, None, None, "not run"
    text = open(log, errors="replace").read()
    energy = ncore = nelectron = None
    match = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", text)
    if match:
        energy = float(match[-1])
    match = re.search(r"replacing\s+(\d+)\s+core electrons", text)
    if match:
        ncore = int(match.group(1))
    match = re.search(r"Number of Electrons\s+NEL\s+\.+\s+(\d+)", text)
    if match:
        nelectron = int(match.group(1))
    if energy is None:
        note = "no energy"
        for line in text.splitlines():
            if "Error" in line or "ABORT" in line:
                note = line.strip()[:60]
                break
        return None, ncore, nelectron, note
    return energy, ncore, nelectron, ""


def run(element, mult, rerun=False):
    source = os.path.join(RECIPES, element, "ccECP")
    basis_file = "%s.%s.orca" % (element, BASIS)
    ecp_file = "%s.ccECP.orca" % element
    for name in (basis_file, ecp_file):
        if not os.path.isfile(os.path.join(source, name)):
            return None, None, None, "missing %s" % name

    work = os.path.join(HERE, element)
    # A finished case is reused, so a long run can be resumed after an
    # interruption without repeating the atoms that already converged.
    if not rerun:
        done = harvest(os.path.join(work, "t.out"))
        if done[0] is not None:
            return done[0], done[1], done[2], "reused"
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    shutil.copy(os.path.join(source, basis_file), work)

    with open(os.path.join(source, ecp_file)) as handle:
        # The converted ECP is a bare NewECP block, so it is pasted into
        # %basis; ORCA does not accept it as a referenced file.
        ecp = "".join("  " + l for l in handle.readlines())

    with open(os.path.join(work, "t.inp"), "w") as handle:
        handle.write(TEMPLATE.format(basis_file=basis_file, ecp=ecp.rstrip(),
                                     mult=mult, element=element))

    log = os.path.join(work, "t.out")
    with open(log, "w") as handle:
        try:
            subprocess.run([ORCA, "t.inp"], cwd=work, stdout=handle,
                           stderr=subprocess.STDOUT, timeout=7200)
        except subprocess.TimeoutExpired:
            return None, None, None, "timed out"

    return harvest(log)


def main():
    if not os.path.isfile(ORCA):
        print("ORCA not found at %s; set ORCA_EXE" % ORCA, file=sys.stderr)
        return 2

    argv = sys.argv[1:]
    if "-h" in argv or "--help" in argv:
        print(__doc__.strip())
        return 0
    rerun = "--rerun" in argv
    wanted = [a for a in argv if not a.startswith("-")]
    cases = [c for c in CASES if not wanted or c[0] in wanted]
    unknown = sorted(set(wanted) - {c[0] for c in CASES})
    if unknown:
        print("not a configured case: %s" % " ".join(unknown), file=sys.stderr)
        return 2

    print("Converted ccECP files vs the library's Pyscf/PBE reference")
    print("basis %s, UKS/PBE, ORCA at %s" % (BASIS, ORCA))
    print()
    print("%-3s %5s %4s %14s %14s %11s  %s"
          % ("El", "NCORE", "e-", "ORCA (Ha)", "reference (Ha)", "difference",
             "verdict"))

    failures = 0
    for element, mult in cases:
        ref = reference(element)
        energy, ncore, nelectron, note = run(element, mult, rerun=rerun)
        if energy is None:
            print("%-3s %5s %4s %14s %14s %11s  FAILED: %s"
                  % (element, ncore if ncore is not None else "-",
                     nelectron if nelectron is not None else "-", "-",
                     "%.5f" % ref if ref is not None else "-", "-", note))
            failures += 1
            continue
        if ref is None:
            print("%-3s %5d %4d %14.6f %14s %11s  no reference"
                  % (element, ncore, nelectron, energy, "-", "-"))
            continue
        difference = energy - ref
        if abs(difference) < 1e-4:
            verdict = "agrees" + (" (reused)" if note == "reused" else "")
        elif abs(difference) < 1e-3:
            verdict = "agrees (1e-3)"
        else:
            verdict = "DIFFERS"
            failures += 1
        print("%-3s %5d %4d %14.6f %14.5f %11.2e  %s"
              % (element, ncore, nelectron, energy, ref, difference, verdict))
        sys.stdout.flush()

    print()
    sys.stdout.flush()
    print("%d of %d cases outside tolerance" % (failures, len(cases)))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
