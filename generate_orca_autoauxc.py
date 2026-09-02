#!/usr/bin/env python3
"""Generate ORCA AutoAux correlation-fitting (AuxC) bases in batch.

Run this script from the pseudopotential library's ``recipes`` directory.  It
visits every ``recipes/Element/Recipe`` directory whose recipe name contains
``ccECP`` (case-insensitive), including variants such as ``ccECP_He_core``,
``ccECP.S``, and ``ccECP-soft``.  Within each recipe it identifies the ECP
GAMESS file from its ``Element-name GEN ncore lmax`` header, treats the other
``*.gamess`` files as orbital bases, converts the ECP to an ORCA ``NewECP``
block, runs AutoAux, and extracts the printed NewAuxCGTO block.

No SCF or MP2 calculation is performed: ORCA is run with NoIter solely as a
basis parser and AutoAux generator.

Shell headers are emitted in the one spelling ORCA accepts, ``S 9``.  Boron and
carbon ship their ``cc-pV{D,T,Q,5}Z`` sets in the Gaussian spelling
``s 9 1.00`` instead, on which ORCA hangs while reading the basis; see
``normalize_shell_header``.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import re
import shlex
import subprocess
import sys
from pathlib import Path


ANGULAR_MOMENTA = "spdfghikl"

# Index in this tuple is the atomic number.  Only the atomic number is needed:
# it lets the script choose a formally valid singlet/doublet for the NoIter
# atom after subtracting the ECP core.  No ground-state atomic configuration is
# implied or required.
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


@dataclass(frozen=True)
class ECPPrimitive:
    coefficient: str
    power: int
    exponent: str


@dataclass(frozen=True)
class CcECP:
    element: str
    n_core: int
    local_l: int
    # GAMESS order: local channel first, then s, p, ..., local_l - 1.
    groups: tuple[tuple[ECPPrimitive, ...], ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Visit every recipes/Element/Recipe directory whose name contains "
            "ccECP, run ORCA AutoAux for every orbital GAMESS basis, and "
            "extract each Auxiliary/C block."
        )
    )
    parser.add_argument(
        "folder",
        nargs="?",
        default=".",
        type=Path,
        help=(
            "Pseudopotential-library recipes directory, or one recipe "
            "directory whose name contains ccECP (default: current directory)."
        ),
    )
    parser.add_argument(
        "--pattern",
        default="*.gamess",
        help=(
            "Glob used inside each ccECP-family directory; the detected ECP file "
            "is excluded automatically (default: *.gamess)."
        ),
    )
    parser.add_argument(
        "--orca",
        default="orca",
        help="ORCA executable or command (default: orca).",
    )
    parser.add_argument(
        "--output-template",
        default="{stem}.AutoAuxC.orca",
        help=(
            "Output filename template. {name}, {stem}, and {element} are "
            "available (default: {stem}.AutoAuxC.orca)."
        ),
    )
    parser.add_argument(
        "--jobs-dir",
        default="autoauxc_jobs",
        help=(
            "Subdirectory created within each ccECP-family recipe for ORCA inputs "
            "and complete outputs (default: autoauxc_jobs)."
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Regenerate an AuxC file even when its output already exists.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Prepare all ORCA job files but do not run ORCA.",
    )
    return parser.parse_args()


def significant_lines(text: str) -> list[str]:
    """Return nonblank, noncomment lines for light format inspection."""
    answer = []
    for line in text.splitlines():
        stripped = line.split("!", 1)[0].split("#", 1)[0].strip()
        if not stripped:
            continue
        answer.append(stripped)
    return answer


def normalize_element(value: str) -> str:
    key = value.strip().upper()
    if key not in ATOMIC_NUMBERS:
        raise ValueError(f"unknown element symbol: {value!r}")
    return ELEMENTS[ATOMIC_NUMBERS[key]]


def parse_number(token: str) -> float:
    return float(token.replace("D", "E").replace("d", "e"))


def clean_ecp_lines(text: str) -> list[str]:
    lines = []
    for raw in text.splitlines():
        line = raw.split("!", 1)[0].split("#", 1)[0].strip()
        if not line or line.upper() in {"$ECP", "$END"}:
            continue
        lines.append(line)
    return lines


def parse_ccecp(text: str, source_name: str = "ccECP input") -> CcECP:
    """Parse a ccECP in GAMESS ``Element-name GEN ncore lmax`` format."""
    lines = clean_ecp_lines(text)
    if not lines:
        raise ValueError(f"{source_name}: empty ECP file")

    header = re.fullmatch(
        r"([A-Za-z]{1,2})(?:-[^\s]+)?\s+GEN\s+(\d+)\s+(\d+)",
        lines[0],
        flags=re.IGNORECASE,
    )
    if header is None:
        raise ValueError(
            f"{source_name}: expected 'Element-name GEN ncore lmax', got "
            f"{lines[0]!r}"
        )

    element = normalize_element(header.group(1))
    n_core = int(header.group(2))
    local_l = int(header.group(3))
    if local_l >= len(ANGULAR_MOMENTA):
        raise ValueError(
            f"{source_name}: local angular momentum {local_l} is unsupported"
        )
    if n_core < 0 or n_core >= ATOMIC_NUMBERS[element.upper()]:
        raise ValueError(
            f"{source_name}: invalid N_core={n_core} for {element}"
        )

    cursor = 1
    groups: list[tuple[ECPPrimitive, ...]] = []
    expected_groups = local_l + 1
    for group_index in range(expected_groups):
        if cursor >= len(lines) or not re.fullmatch(r"\d+", lines[cursor]):
            raise ValueError(
                f"{source_name}: missing primitive count for ECP group "
                f"{group_index + 1} of {expected_groups}"
            )
        count = int(lines[cursor])
        cursor += 1
        if count < 1:
            raise ValueError(f"{source_name}: ECP groups cannot be empty")

        primitives = []
        for primitive_index in range(count):
            if cursor >= len(lines):
                raise ValueError(
                    f"{source_name}: ECP group {group_index + 1} is truncated"
                )
            fields = lines[cursor].split()
            cursor += 1
            if len(fields) != 3:
                raise ValueError(
                    f"{source_name}: expected 'coefficient power exponent', "
                    f"got {lines[cursor - 1]!r}"
                )
            coefficient, power_token, exponent = fields
            try:
                parse_number(coefficient)
                power = int(power_token)
                parse_number(exponent)
            except ValueError as exc:
                raise ValueError(
                    f"{source_name}: invalid primitive {lines[cursor - 1]!r}"
                ) from exc
            primitives.append(ECPPrimitive(coefficient, power, exponent))
        groups.append(tuple(primitives))

    if cursor != len(lines):
        raise ValueError(
            f"{source_name}: unexpected trailing content: {lines[cursor]!r}"
        )
    return CcECP(element, n_core, local_l, tuple(groups))


def read_ccecp(path: Path) -> CcECP:
    return parse_ccecp(path.read_text(encoding="utf-8"), path.name)


def orca_number(token: str) -> str:
    """ORCA accepts E notation; normalize Fortran D notation when present."""
    return token.replace("D", "E").replace("d", "e")


def format_orca_ecp(ecp: CcECP) -> str:
    """Convert GAMESS ccECP ordering and primitive columns to ORCA syntax."""
    local_label = ANGULAR_MOMENTA[ecp.local_l]
    lines = [
        f"      NewECP {ecp.element}",
        f"        N_core {ecp.n_core}",
        f"        lmax {local_label}",
        "",
    ]

    # GAMESS: local group, s group, p group, ...
    # ORCA:    s group, p group, ..., local group
    ordered = [
        (ANGULAR_MOMENTA[l_value], ecp.groups[l_value + 1])
        for l_value in range(ecp.local_l)
    ]
    ordered.append((local_label, ecp.groups[0]))

    for label, primitives in ordered:
        lines.append(f"        {label} {len(primitives)}")
        for index, primitive in enumerate(primitives, start=1):
            lines.append(
                f"          {index} {orca_number(primitive.exponent)} "
                f"{orca_number(primitive.coefficient)} {primitive.power}"
            )
        lines.append("")
    lines.append("      end")
    return "\n".join(lines)


SHELL_HEADER = re.compile(r"(\s*)([A-Za-z]{1,2})(\s+)(\d+)(?:\s+(\S+))?\s*")


def normalize_shell_header(line: str) -> str:
    """Rewrite a shell header into the only spelling ORCA's GAMESS reader takes.

    ORCA needs an upper-case shell letter and exactly two fields.  Anything
    else makes it hang while reading the basis, with no error message: with
    ORCA 6.1.0, ``s 9 1.00``, ``S 9 1.00`` and ``s 9`` all fail and only
    ``S 9`` works.

    ``B`` and ``C`` ship their ``cc-pV{D,T,Q,5}Z`` sets in the Gaussian
    spelling ``s 9 1.00``, which is why those 8 sets otherwise get no
    auxiliary basis.

    Only the spelling changes.  A scale factor multiplies the exponents of its
    shell by its square, so one that is not 1 is refused rather than dropped;
    every occurrence in the library is exactly ``1.00``.  A header already in
    GAMESS form is returned untouched, spacing included.
    """
    match = SHELL_HEADER.fullmatch(line)
    if match is None:
        return line
    indent, letter, gap, count, scale = match.groups()
    if scale is None and not letter.islower():
        return line
    if scale is not None and parse_number(scale) != 1.0:
        raise ValueError(
            f"shell header {line.strip()!r} carries scale factor {scale}, "
            f"which would rescale the exponents; refusing to drop it"
        )
    return f"{indent}{letter.upper()}{gap}{count}"


def prepare_gamess_basis(source: Path, element: str) -> str:
    """Add the element header required by ORCA's external GTO reader.

    The project files may begin directly with shell data (for example,
    ``S 13``).  ORCA needs an element identifier before those shells. Existing
    $DATA/$END wrappers and symbol/full-name element headers are normalized.
    """
    original = source.read_text(encoding="utf-8")
    lines = original.splitlines()

    # Remove wrappers so that exactly one clean wrapper can be written.
    inner = [
        line
        for line in lines
        if line.strip().upper() not in {"$DATA", "$END"}
    ]
    sig = significant_lines("\n".join(inner))
    if not sig:
        raise ValueError("file contains no basis data")

    # Remove a pre-existing symbol or full-name element header.  A shell header
    # such as ``F 1`` has two tokens and therefore cannot be mistaken for it.
    first_sig = sig[0].rstrip(";")
    if re.fullmatch(r"[A-Za-z]+", first_sig):
        removed = False
        cleaned = []
        for line in inner:
            content = line.split("!", 1)[0].split("#", 1)[0].strip().rstrip(";")
            if not removed and content == first_sig:
                removed = True
                continue
            cleaned.append(line)
        inner = cleaned

    inner = [normalize_shell_header(line) for line in inner]
    body = "\n".join(inner).strip()
    shell_match = re.search(
        r"(?m)^\s*[SPDFGHIKL]\s+\d+\s*$", body
    )
    if shell_match is None:
        raise ValueError("no GAMESS-style shell header such as 'S 13' was found")

    return f"$DATA\n{element}\n\n{body}\n$END\n"


def make_orca_input(element: str, ecp: CcECP) -> str:
    if ecp.element != element:
        raise ValueError(
            f"orbital basis is for {element}, but ECP header is for {ecp.element}"
        )
    valence_electrons = ATOMIC_NUMBERS[element.upper()] - ecp.n_core
    multiplicity = 1 if valence_electrons % 2 == 0 else 2
    ecp_block = format_orca_ecp(ecp)
    return f"""\
! NoIter PrintBasis

%maxcore 1000
%base "autoauxc"

%basis
      GTOName = "orbital.gamess"

      NewAuxCGTO {element}
        "AutoAux"
      end

{ecp_block}
end

* xyz 0 {multiplicity}
{element} 0.0 0.0 0.0
*
"""


def extract_auxc_block(output: str, element: str = "Fe") -> str:
    """Extract the AuxC basis in ORCA input form from a complete ORCA output."""
    marker = re.search(
        r"AUXILIARY/C BASIS SET IN INPUT FORMAT",
        output,
        flags=re.IGNORECASE,
    )
    # ORCA prints a dashed underline immediately after the section title, so
    # search everything following the title rather than stopping at dashes.
    search_text = output[marker.end() :] if marker else output

    start = re.search(
        rf"(?mi)^\s*NewAuxCGTO\s+{re.escape(element)}\s*$", search_text
    )
    if start is None:
        raise ValueError(
            f"could not find 'NewAuxCGTO {element}' in the ORCA output"
        )

    tail = search_text[start.start() :]
    end = re.search(r"(?mi)^\s*end\s*;?\s*$", tail)
    if end is None:
        raise ValueError("found NewAuxCGTO, but not its closing 'end;' line")

    block = tail[: end.end()].strip()
    if re.search(r"(?mi)^\s*[SPDFGHIKL]\s+\d+\s*$", block) is None:
        raise ValueError("the extracted NewAuxCGTO block contains no shells")

    return (
        f"# Auxiliary/C basis set generated by ORCA AutoAux for {element}\n"
        f"{block}\n"
    )


def safe_job_name(stem: str) -> str:
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", stem)
    return name.strip("._") or "basis"


def is_ccecp_recipe_name(name: str) -> bool:
    return "ccecp" in name.lower()


def discover_ccecp_directories(root: Path) -> list[Path]:
    """Return Element/Recipe directories whose recipe name contains ccECP."""
    if is_ccecp_recipe_name(root.name):
        return [root]
    return sorted(
        path
        for path in root.glob("*/*")
        if path.is_dir() and is_ccecp_recipe_name(path.name)
    )


def find_ccecp_file(
    gamess_files: list[Path], element: str
) -> tuple[Path, CcECP]:
    """Identify the one ECP file by parsing content rather than its filename."""
    matches: list[tuple[Path, CcECP]] = []
    wrong_elements: list[tuple[Path, str]] = []
    for path in gamess_files:
        try:
            ecp = read_ccecp(path)
        except (OSError, ValueError):
            continue
        if ecp.element == element:
            matches.append((path, ecp))
        else:
            wrong_elements.append((path, ecp.element))

    if not matches:
        detail = ""
        if wrong_elements:
            found = ", ".join(
                f"{path.name} ({found_element})"
                for path, found_element in wrong_elements
            )
            detail = f"; parsable ECP files for other elements: {found}"
        raise ValueError(
            f"no GAMESS ECP with a {element}-... GEN header was found{detail}"
        )
    if len(matches) > 1:
        names = ", ".join(path.name for path, _ in matches)
        raise ValueError(f"multiple GAMESS ECP files were found: {names}")
    return matches[0]


def run_one(
    source: Path,
    recipe_dir: Path,
    jobs_root: Path,
    orca_command: str,
    output_template: str,
    overwrite: bool,
    dry_run: bool,
    element: str,
    ecp: CcECP,
) -> str:
    output_name = output_template.format(
        name=source.name,
        stem=source.stem,
        element=element,
    )
    output_path = recipe_dir / output_name
    if output_path.resolve() == source.resolve():
        raise ValueError("output template would overwrite the source basis file")
    if output_path.exists() and not overwrite and not dry_run:
        return f"SKIP  {source.name}: {output_path.name} already exists"

    job_dir = jobs_root / safe_job_name(source.stem)
    job_dir.mkdir(parents=True, exist_ok=True)
    (job_dir / "orbital.gamess").write_text(
        prepare_gamess_basis(source, element), encoding="utf-8"
    )
    (job_dir / "ccECP.orca").write_text(
        format_orca_ecp(ecp) + "\n", encoding="utf-8"
    )
    input_path = job_dir / "autoauxc.inp"
    input_path.write_text(make_orca_input(element, ecp), encoding="utf-8")

    if dry_run:
        return f"READY {source.name}: {input_path}"

    command = shlex.split(orca_command) + [input_path.name]
    completed = subprocess.run(
        command,
        cwd=job_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    full_output_path = job_dir / "autoauxc.out"
    full_output_path.write_text(completed.stdout, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"ORCA exited with status {completed.returncode}; see {full_output_path}"
        )

    block = extract_auxc_block(completed.stdout, element)
    output_path.write_text(block, encoding="utf-8")
    return f"DONE  {source.name} -> {output_path.name}"


def main() -> int:
    args = parse_args()
    root = args.folder.expanduser().resolve()
    if not root.is_dir():
        print(f"ERROR: not a directory: {root}", file=sys.stderr)
        return 2

    recipe_dirs = discover_ccecp_directories(root)
    if not recipe_dirs:
        print(
            f"ERROR: no Element/Recipe directories containing 'ccECP' found "
            f"below {root}",
            file=sys.stderr,
        )
        return 2

    failed = 0
    processed = 0
    skipped = 0
    for recipe_dir in recipe_dirs:
        try:
            element = normalize_element(recipe_dir.parent.name)
        except ValueError as exc:
            failed += 1
            print(f"FAIL  {recipe_dir}: {exc}", file=sys.stderr)
            continue

        all_gamess = sorted(
            path for path in recipe_dir.glob("*.gamess") if path.is_file()
        )
        if not all_gamess:
            skipped += 1
            print(f"SKIP  {recipe_dir}: no GAMESS files")
            continue

        try:
            ecp_path, ecp = find_ccecp_file(all_gamess, element)
        except (OSError, ValueError) as exc:
            failed += 1
            print(f"FAIL  {recipe_dir}: {exc}", file=sys.stderr)
            continue

        sources = sorted(
            path
            for path in recipe_dir.glob(args.pattern)
            if path.is_file() and path.resolve() != ecp_path.resolve()
        )
        if not sources:
            skipped += 1
            print(
                f"SKIP  {recipe_dir}: found {ecp_path.name}, but no orbital "
                f"files matching {args.pattern!r}"
            )
            continue

        jobs_root = recipe_dir / args.jobs_dir
        jobs_root.mkdir(parents=True, exist_ok=True)
        for source in sources:
            processed += 1
            try:
                print(
                    run_one(
                        source=source,
                        recipe_dir=recipe_dir,
                        jobs_root=jobs_root,
                        orca_command=args.orca,
                        output_template=args.output_template,
                        overwrite=args.overwrite,
                        dry_run=args.dry_run,
                        element=element,
                        ecp=ecp,
                    )
                )
            except (OSError, RuntimeError, ValueError, KeyError) as exc:
                failed += 1
                print(f"FAIL  {source}: {exc}", file=sys.stderr)

    print(
        f"Visited {len(recipe_dirs)} ccECP-family recipe directories; "
        f"processed {processed} orbital basis file(s); skipped {skipped} "
        f"empty/ECP-only recipe(s); {failed} failure(s)"
    )
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
