# Canonical ccECP for CFOUR / MRCC

`GENBAS` and `ECPDATA` generated from the **canonical `ccECP`** family of
[PseudopotentialLibrary.org](https://pseudopotentiallibrary.org).

> **This output represents ONLY the canonical ccECP family from
> PseudopotentialLibrary.org.**  No `ccECP-soft`, `soft-ccECP-deprecated`,
> `eCEPP`, `CEPP`, `BFD`, `RRKJ` or `TM` data is included, and no
> grid-based AREP/SOREP files are included. Spin-orbit SOREP data is excluded;
> the conversion covers the scalar-relativistic ccECP only.

## Quick start

The accompanying `convert_ccecp_to_cfour.py` is a standalone Python script.
CFOUR and MRCC are **not required to generate the files**. Python 3.8+ is
sufficient for the converter itself. The optional `basis_set_exchange` package
is strongly recommended because it provides independent parsing and ECP
convention checks. Git is optional and is used only to record repository
source.

Example:
```bash
git clone https://github.com/QMCPACK/pseudopotentiallibrary.git
python3 -m venv .venv
source .venv/bin/activate
python -m pip install basis_set_exchange
python convert_ccecp_to_cfour.py \
  --library ./pseudopotentiallibrary \
  --out ./ccECP_cfour_output
```

`--library` may point either to the repository root or directly to its
`recipes` directory. If omitted, the converter checks
`$PSEUDOPOTENTIAL_LIBRARY`, the script directory, the current directory, and a
`pseudopotentiallibrary` subdirectory in those locations. Use `--no-bse` only
when Basis Set Exchange is unavailable. Run `python convert_ccecp_to_cfour.py
--help` for all options.

By default, orbital-basis labels receive the suffix `-ccECP` to distinguish the
ccECP-optimized basis from installed basis sets with the same conventional name.
ECP labels remain exactly `ccECP`. See section 6 for details.

The source Pseudopotential Library is treated as read-only. The output directory
must be outside that source tree. On Windows, the same commands can be run in
WSL or another Python environment.

## 1. Source and Conversion Details

| item | value |
|---|---|
| Pseudopotential Library commit | `0e81cb49cbf817558d9a4930dea8b461cc1ca2bd` |
| commit date | 2026-06-29 14:40:38 -0400 |
| converter | `convert_ccecp_to_cfour.py` 1.1.2 |
| Python | 3.12.3 |
| basis_set_exchange | 0.12 (independent cross-check only) |

## 2. Canonical ccECP selection rule

A recipe is in scope **iff** it lives in `recipes/<ELEMENT>/ccECP/`.
Supported **basis-set** renderings are `.gamess`, `.nwchem`, and `.gaussian`.
Supported **ECP** renderings are the native `<El>.ccECP` text plus `.gamess`,
`.nwchem`, `.gaussian`, and `.molpro`. Molpro basis files and all `.dirac`,
`.xml`, `.upf`, `.rpt`, `*_deprecated`, `.AREP.*`, and `.SOREP.*` files are
deliberately ignored. No numerical-grid or spin-orbit data can therefore leak
into the scalar-relativistic output.

Recipe families that exist in the library and are **excluded**:

* `CEPP` (14 elements)
* `RRKJ_PRB_93_075143` (10 elements)
* `TM_PRB_93_075143` (10 elements)
* `ccECP-soft` (7 elements)
* `ccECP.S` (2 elements)
* `ccECP_28_core` (2 elements)
* `ccECP_36_core` (1 elements)
* `ccECP_46_core` (1 elements)
* `ccECP_He_core` (8 elements)
* `ccECP_reg` (2 elements)
* `eCEPP` (15 elements)
* `soft-ccECP-deprecated` (4 elements)

Element directories that exist but have **no** `ccECP` recipe: `In`, `Sr`.

## 3. Source representations used

Selection is `--ecp-source gamess` / `--basis-source auto`.  In `auto` mode every
available analytic rendering is parsed independently, all renderings must
be mutually *consistent* (their difference must not exceed one unit in
the last written place of the coarser one), and the rendering carrying the
largest number of significant digits is used.  The rendering actually used
is recorded per entry in `ecp_name_map.csv` / `basis_name_map.csv`.

| kind | rendering | entries |
|---|---|---|
| ecp | `.gamess` | 65 |
| basis | `.gamess` | 841 |

The validated conversion selected the GAMESS rendering for all 841 basis
sets. Molpro basis renderings are not parsed by this converter; Molpro is used
only as an available ECP cross-check.

## 4. GENBAS conversion procedure

1. Parse the source file into one record per *contracted function*,
   keeping every exponent and coefficient as its **source decimal string**
   (no float round-trip ever happens).
2. Group the contracted functions by angular momentum, s first.
3. For each L build the union of the primitive exponents, sort it in descending
   exponent order, and place each contracted function in one column of the
   coefficient matrix required by CFOUR.
   Two primitives are identified only when their decimals are *exactly*
   equal; matrix positions a given function does not use are exact zero.
4. Emit the CFOUR *new-format* entry
   (`ELEM:name` / comment / blank / NS / L / NC / NE / blank /
   per-shell exponents / blank / NE x NC coefficient matrix / blank).
   The numeric records of a GENBAS entry are read list-directed
   (Fortran `*`). Numerical values are preserved exactly; scientific notation
   may be rendered in equivalent fixed-point form. The new-format list-directed
   layout avoids the range and precision limits of the legacy fixed-width
   GENBAS representation.

Nothing is normalised, recontracted, merged, deleted or refitted.

## 5. ECPDATA conversion procedure

1. Parse every available scalar rendering.  The MOLPRO `!*-so` blocks and
   the NWChem `so ... end` section are dropped explicitly.
2. Read `NCORE` and `LMAX` from the source header of each rendering and
   require them to agree across renderings.
3. Emit `local channel first`, labelled with the letter of `l = LMAX`,
   then `s-<lmax>`, `p-<lmax>`, ... for `l = 0 .. LMAX-1`.
4. Columns are `coefficient   N   alpha`, with `N` copied through
   unchanged from the source.

### Radial-power convention

```
U_l(r) = sum_m  c_m  r**(N_m - 2)  exp(-alpha_m r**2)
```

The integer in the middle ECPDATA column is copied unchanged from the
GAMESS / Gaussian / NWChem / MOLPRO ECP representations. In the documented
conversion, Basis Set Exchange was enabled and the convention was independently
cross-checked against literature ECPs in CFOUR's shipped `ECPDATA` and the
MolSSI Basis Set Exchange.

Two independent CFOUR/BSE convention checks were recorded. The Hay-Wadt
`LANL2DZ` comparison agrees exactly. The Stuttgart `ECP-28-MWB` comparison
agrees within 0.44 units in the last written place, consistent with rounding in
the BSE copy. When the converter is run with `--no-bse`, this optional section
is reported as skipped rather than emitting an empty results table.

The format convention is: the middle integer is copied
unchanged from the analytic source representations and corresponds to the
standard radial factor `r**(N-2)`. For example, CFOUR's shipped Hay-Wadt Br
entry contains `-28.0  1  213.6143969`, which represents an `r**-1` term.

### Local-channel convention

The library's `.gamess` rendering writes `LMAX+1` blocks: the local
potential first (`ul` in the NWChem rendering) followed by `l = 0 ..
LMAX-1`.  CFOUR expects exactly the same ordering, with the local block
labelled by the letter of `l = LMAX` and the semi-local blocks labelled
`<l>-<lmax>` because they represent `U_l - U_LMAX`.  The mapping is
therefore a direct pass-through.

### NCORE is derived from the data, never from a filename

For every element the coefficients of the `N = 1` (i.e. `r**-1`) terms of
the local channel are summed and required to equal `Z - NCORE`.  This
independent check passes for all converted elements (nitrogen, for
instance, splits `Zeff = 5` into two terms, `3.25 + 1.75`).

## 6. Naming policy

| Pseudopotential Library | CFOUR/MRCC label |
|---|---|
| `Fe.ccECP` | `Fe:ccECP` |
| `Fe.cc-pVQZ` | `Fe:cc-pVQZ-ccECP` |
| `Fe.aug-cc-pCV5Z` | `Fe:aug-cc-pCV5Z-ccECP` |
| `Fe.aug-cc-pCV5Z.AutoAuxC.orca` | `Fe:aug-cc-pCV5Z-RI-ccECP` |

The element symbol keeps its periodic-table spelling, and the source basis-set
name is retained as the stem of each orbital-basis label. The **ECP** label is
exactly the family name, `ccECP`. Nothing is renamed to `CC-ECP` or `ccecp`,
and no prefix is added.

Orbital-basis labels carry the suffix **`-ccECP`**. This is deliberate, and it
is the one naming transformation in the output.

Installed MRCC/CFOUR basis libraries can already contain entries with the same
conventional names, such as `cc-pVTZ` or `aug-cc-pCVTZ`, but referring to a
different basis. In a validation against the MRCC 26.1.1 `BASIS` library, 290
of the 841 source `(element, basis)` labels in this ccECP data set were already
defined. Because a working-directory `GENBAS` can take precedence over the
installed library, an unsuffixed label can silently change which basis a shared
input selects. MRCC already distinguishes related sets by suffixes such as
`-PP`, `-DK`, and `-RI`, so `-ccECP` follows the same disambiguation pattern.

With the default suffix, the generated ccECP basis labels are distinct from the
corresponding unsuffixed installed names.

Auxiliary (density-fitting) sets take `-RI` between the basis name and the
suffix, giving `<element>:<basis>-RI-ccECP`. `-RI` is MRCC's own name for a
correlation-fitting set, as against `-RI-JK` for an SCF-fitting one, and the
converted sets are ORCA AutoAux `/C` correlation-fitting sets. The `-ccECP`
suffix is still needed after it: MRCC already ships entries such as
`Ag:aug-cc-pV5Z-RI`. Measured against the MRCC 26.1.1 `BASIS` library (5261
labels), the auxiliary labels written here collide **0** times, and would
collide **187** times were `-ccECP` dropped. Pass `--host-basis-library` to
repeat that measurement against any installation.

The source basis name remains intact as the label stem: `cc-pVQZ` becomes
`cc-pVQZ-ccECP`.
`basis_name_map.csv` records `source_name`, `label_suffix` and `cfour_name` for
every entry, and `naming_audit.csv` verifies that each label is exactly
`<element>:<library name><suffix>`; failures: **0**.

Run the converter with `--basis-suffix ""` to use the exact source basis names
after checking the target installation for label collisions.

For CFOUR builds that fold the ZMAT atom string to upper case before searching
the library, `GENBAS.upper` and `ECPDATA.upper` carry the same data with the
element symbol upper-cased (`AG:aug-cc-pCV5Z-ccECP`); nothing after the colon
changes.

## 7. Contents

* **65** canonical ccECP elements
* **65** `ECPDATA` entries (one `ccECP` per element)
* **841** `GENBAS` entries
* **841** auxiliary (RI) entries, one per orbital basis set
* **20** distinct source basis-set names: `aug-cc-pCV5Z`, `aug-cc-pCV6Z`, `aug-cc-pCVDZ`, `aug-cc-pCVQZ`, `aug-cc-pCVTZ`, `aug-cc-pV5Z`, `aug-cc-pV6Z`, `aug-cc-pVDZ`, `aug-cc-pVQZ`, `aug-cc-pVTZ`, `cc-pCV5Z`, `cc-pCV6Z`, `cc-pCVDZ`, `cc-pCVQZ`, `cc-pCVTZ`, `cc-pV5Z`, `cc-pV6Z`, `cc-pVDZ`, `cc-pVQZ`, `cc-pVTZ`

### Elements

| El | Z | NCORE | e- explicit | LMAX | local | non-local | # basis | basis names |
|---|---|---|---|---|---|---|---|---|
| Ag | 47 | 28 | 19 | 3 | f | s-f;p-f;d-f | 16 | aug-cc-pCV5Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Al | 13 | 10 | 3 | 2 | d | s-d;p-d | 10 | aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Ar | 18 | 10 | 8 | 2 | d | s-d;p-d | 10 | aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| As | 33 | 28 | 5 | 3 | f | s-f;p-f;d-f | 10 | aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Au | 79 | 60 | 19 | 4 | g | s-g;p-g;d-g;f-g | 16 | aug-cc-pCV5Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| B | 5 | 2 | 3 | 1 | p | s-p | 10 | aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Ba | 56 | 46 | 10 | 4 | g | s-g;p-g;d-g;f-g | 16 | aug-cc-pCV5Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Be | 4 | 2 | 2 | 1 | p | s-p | 10 | aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Bi | 83 | 78 | 5 | 4 | g | s-g;p-g;d-g;f-g | 10 | aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Br | 35 | 28 | 7 | 3 | f | s-f;p-f;d-f | 10 | aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| C | 6 | 2 | 4 | 1 | p | s-p | 10 | aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Ca | 20 | 10 | 10 | 2 | d | s-d;p-d | 20 | aug-cc-pCV5Z aug-cc-pCV6Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCV6Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Cd | 48 | 28 | 20 | 3 | f | s-f;p-f;d-f | 16 | aug-cc-pCV5Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Ce | 58 | 46 | 12 | 4 | g | s-g;p-g;d-g;f-g | 12 | aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pVDZ cc-pVQZ cc-pVTZ |
| Cl | 17 | 10 | 7 | 2 | d | s-d;p-d | 10 | aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Co | 27 | 10 | 17 | 2 | d | s-d;p-d | 16 | aug-cc-pCV5Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Cr | 24 | 10 | 14 | 2 | d | s-d;p-d | 16 | aug-cc-pCV5Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Cs | 55 | 46 | 9 | 3 | f | s-f;p-f;d-f | 16 | aug-cc-pCV5Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Cu | 29 | 10 | 19 | 2 | d | s-d;p-d | 16 | aug-cc-pCV5Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Eu | 63 | 46 | 17 | 4 | g | s-g;p-g;d-g;f-g | 12 | aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pVDZ cc-pVQZ cc-pVTZ |
| F | 9 | 2 | 7 | 1 | p | s-p | 10 | aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Fe | 26 | 10 | 16 | 2 | d | s-d;p-d | 16 | aug-cc-pCV5Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Ga | 31 | 28 | 3 | 3 | f | s-f;p-f;d-f | 10 | aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Gd | 64 | 46 | 18 | 4 | g | s-g;p-g;d-g;f-g | 12 | aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pVDZ cc-pVQZ cc-pVTZ |
| Ge | 32 | 28 | 4 | 3 | f | s-f;p-f;d-f | 10 | aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| H | 1 | 0 | 1 | 1 | p | s-p | 10 | aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| He | 2 | 0 | 2 | 1 | p | s-p | 10 | aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| I | 53 | 46 | 7 | 3 | f | s-f;p-f;d-f | 10 | aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Ir | 77 | 60 | 17 | 4 | g | s-g;p-g;d-g;f-g | 16 | aug-cc-pCV5Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| K | 19 | 10 | 9 | 2 | d | s-d;p-d | 20 | aug-cc-pCV5Z aug-cc-pCV6Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCV6Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Kr | 36 | 28 | 8 | 3 | f | s-f;p-f;d-f | 10 | aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| La | 57 | 46 | 11 | 4 | g | s-g;p-g;d-g;f-g | 12 | aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pVDZ cc-pVQZ cc-pVTZ |
| Li | 3 | 2 | 1 | 1 | p | s-p | 8 | aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Mg | 12 | 10 | 2 | 2 | d | s-d;p-d | 10 | aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Mn | 25 | 10 | 15 | 2 | d | s-d;p-d | 16 | aug-cc-pCV5Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Mo | 42 | 28 | 14 | 3 | f | s-f;p-f;d-f | 16 | aug-cc-pCV5Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| N | 7 | 2 | 5 | 1 | p | s-p | 10 | aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Na | 11 | 10 | 1 | 2 | d | s-d;p-d | 8 | aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Nb | 41 | 28 | 13 | 3 | f | s-f;p-f;d-f | 16 | aug-cc-pCV5Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Ne | 10 | 2 | 8 | 1 | p | s-p | 5 | cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Ni | 28 | 10 | 18 | 2 | d | s-d;p-d | 16 | aug-cc-pCV5Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| O | 8 | 2 | 6 | 1 | p | s-p | 10 | aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| P | 15 | 10 | 5 | 2 | d | s-d;p-d | 10 | aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Pb | 82 | 78 | 4 | 4 | g | s-g;p-g;d-g;f-g | 8 | aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Pd | 46 | 28 | 18 | 3 | f | s-f;p-f;d-f | 16 | aug-cc-pCV5Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Pt | 78 | 60 | 18 | 4 | g | s-g;p-g;d-g;f-g | 16 | aug-cc-pCV5Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Rb | 37 | 28 | 9 | 3 | f | s-f;p-f;d-f | 20 | aug-cc-pCV5Z aug-cc-pCV6Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCV6Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Re | 75 | 60 | 15 | 4 | g | s-g;p-g;d-g;f-g | 16 | aug-cc-pCV5Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Rh | 45 | 28 | 17 | 3 | f | s-f;p-f;d-f | 16 | aug-cc-pCV5Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Ru | 44 | 28 | 16 | 3 | f | s-f;p-f;d-f | 16 | aug-cc-pCV5Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| S | 16 | 10 | 6 | 2 | d | s-d;p-d | 10 | aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Sb | 51 | 46 | 5 | 3 | f | s-f;p-f;d-f | 10 | aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Sc | 21 | 10 | 11 | 2 | d | s-d;p-d | 16 | aug-cc-pCV5Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Se | 34 | 28 | 6 | 3 | f | s-f;p-f;d-f | 10 | aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Si | 14 | 10 | 4 | 2 | d | s-d;p-d | 10 | aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Sn | 50 | 46 | 4 | 3 | f | s-f;p-f;d-f | 10 | aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Ta | 73 | 60 | 13 | 4 | g | s-g;p-g;d-g;f-g | 16 | aug-cc-pCV5Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Tb | 65 | 46 | 19 | 4 | g | s-g;p-g;d-g;f-g | 12 | aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pVDZ cc-pVQZ cc-pVTZ |
| Te | 52 | 46 | 6 | 3 | f | s-f;p-f;d-f | 10 | aug-cc-pV5Z aug-cc-pV6Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pV5Z cc-pV6Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Ti | 22 | 10 | 12 | 2 | d | s-d;p-d | 16 | aug-cc-pCV5Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| V | 23 | 10 | 13 | 2 | d | s-d;p-d | 16 | aug-cc-pCV5Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| W | 74 | 60 | 14 | 4 | g | s-g;p-g;d-g;f-g | 16 | aug-cc-pCV5Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Y | 39 | 28 | 11 | 3 | f | s-f;p-f;d-f | 16 | aug-cc-pCV5Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Zn | 30 | 10 | 20 | 2 | d | s-d;p-d | 16 | aug-cc-pCV5Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |
| Zr | 40 | 28 | 12 | 3 | f | s-f;p-f;d-f | 16 | aug-cc-pCV5Z aug-cc-pCVDZ aug-cc-pCVQZ aug-cc-pCVTZ aug-cc-pV5Z aug-cc-pVDZ aug-cc-pVQZ aug-cc-pVTZ cc-pCV5Z cc-pCVDZ cc-pCVQZ cc-pCVTZ cc-pV5Z cc-pVDZ cc-pVQZ cc-pVTZ |

## 8. Anomalies and notes

* `In` and `Sr` have element directories in the library but **no**
  `ccECP` recipe, so they are absent from this output.
* `Ne` has no `aug-` sets in the canonical library; `Li`, `Na` and `Pb`
  have no `6Z` and no core-valence sets.  Missing members of a series are
  **not** invented.
* `H` and `He` have `NCORE = 0`: the ccECP is a shape-consistent
  potential with all electrons explicit.  Their non-local `s` channel is
  a single term with coefficient exactly `0.0`, which is preserved.
* The `B` and `C` `cc-pV{D,T,Q,5}Z` `.gamess` files use lower-case
  Gaussian-style shell headers (`s 9 1.00`) instead of the usual `S 9`;
  the parser handles both spellings.
* `Ce` has a primitive with an exactly zero contraction coefficient in
  several `f` shells; it is preserved because the exponent is shared with
  the other `f` functions.
* The `.gaussian` renderings pad numbers with trailing zeros and, for
  nitrogen, round two orbital exponents (`9.3345971` -> `9.334597`);
  the `.gamess`/`.nwchem` renderings carry the extra digit and are used.
* Molpro **basis-set** renderings are not used by this converter. Basis
  `auto` selection compares the supported `.gamess`, `.nwchem`, and `.gaussian`
  renderings. Molpro is supported for ECP cross-checking.
* `cc-pV6Z`/`cc-pCV6Z` sets contain `i` functions (L = 6).  Make sure the
  CFOUR/MRCC build supports that angular momentum.

## 9. MRCC ECP compatibility

MRCC 25.1.1 and 26.1.1 read the `NCORE`/`LMAX` ECPDATA record using the fixed
Fortran format

```fortran
read(gbasfile,"(11x,i3,11x,i1)") ncorecp(iatoms),lmax
```

Therefore `NCORE` must occupy columns 12-14 and `LMAX` column 26. The converter
uses one formatter, `format_ncore_lmax_line()`, for every ECPDATA entry and
checks the exact columns with `selftest_ncore_lmax()` on every run. This layout
is also valid CFOUR input.

```text
         1         2         3
123456789012345678901234567890
    NCORE =  0    LMAX = 1   <- H
    NCORE =  2    LMAX = 1   <- C
    NCORE = 10    LMAX = 2   <- Fe
    NCORE = 28    LMAX = 3   <- Ru
    NCORE = 60    LMAX = 4   <- Pt
```

### MRCC validation

Representative MRCC 26.1.1 and PySCF UKS/PBE/cc-pVTZ calculations using the
same converted basis/ECP data give:

| El | NCORE | LMAX | MRCC (Ha) | PySCF full (Ha) | MRCC-PySCF |
|---|---:|---:|---:|---:|---:|
| Li | 2 | 1 | -0.201399 | -0.201399 | -3.6e-09 |
| Be | 2 | 1 | -0.995253 | -0.995253 | -2.2e-11 |
| B | 2 | 1 | -2.602291 | -2.602326 | +3.5e-05 |
| C | 2 | 1 | -5.409314 | -5.409318 | +4.8e-06 |
| N | 2 | 1 | -9.755863 | -9.755863 | +5.4e-09 |
| O | 2 | 1 | -15.883639 | -15.883673 | +3.5e-05 |
| F | 2 | 1 | -24.193257 | -24.193304 | +4.8e-05 |
| Ne | 2 | 1 | -35.022840 | -35.022840 | -7.4e-07 |
| Si | 10 | 2 | -3.760722 | -3.760722 | +1.2e-08 |
| Cr | 10 | 2 | -86.658120 | -86.658120 | -2.6e-10 |
| Pt | 60 | 4 | -119.662468 | -119.662468 | -7.2e-08 |

The small residual differences for B, O, and F are consistent with numerical
DFT grid/convergence differences between programs. Open-shell atoms can also
converge to different SCF solutions, so total-energy disagreements should be
diagnosed at the SCF level before being attributed to ECP parsing.

### `NCORE = 0` in MRCC

Executable validation shows that MRCC 25.1.1 and 26.1.1 do not correctly handle
ECP entries with `NCORE = 0`. This affects the canonical ccECPs of **H** and
**He**, which keep all electrons explicit:

| El | NCORE | LMAX | MRCC (Ha) | PySCF full (Ha) | MRCC-PySCF |
|---|---:|---:|---:|---:|---:|
| H | 0 | 1 | -0.499604 | -0.499866 | +2.6e-04 |
| He | 0 | 1 | -2.886897 | -2.892406 | +5.5e-03 |

The converter writes the correct `NCORE = 0` ECPDATA entries without a
program-specific numerical workaround. Validate H and He independently before
production use with these MRCC releases.

## 10. Validation

`validation_report.csv` contains one row per conversion/consistency check, and
`summary.json` records the status counts for the current run.

Checks performed for **every basis set**: source parse of all available
renderings and mutual consistency; contracted-function count; primitive count;
angular-momentum grouping and ordering; exponent ordering and near-degeneracy;
highest-L retention; and a full write -> independent re-read -> exact-decimal
comparison round trip. Plus cross-basis checks that every `aug-X` retains all
primitives of `X` and that every `cc-pCVnZ` retains all primitives of
`cc-pVnZ`.

For the documented conversion, Basis Set Exchange was enabled, so supported
Gaussian basis renderings were also parsed independently with
`basis_set_exchange`. Runs using `--no-bse` omit that optional check and state
so explicitly in the generated README.

Checks performed for **every ECP**: channel count vs `LMAX`; `NCORE`
agreement across all renderings; `NCORE` re-derived from the local
channel's `r**-1` coefficients; radial-power range; absence of `r**-2`
terms; spin-orbit exclusion; cross-rendering numerical agreement; and a
write -> independent re-read -> exact-decimal comparison round trip.

Cross-basis hierarchy diagnostics include 46 `INFO` cases where related source
basis families are not strict primitive supersets. These are preserved and are not
treated as conversion failures.

**No failures**: every canonical ccECP element and every canonical
ccECP basis set converted and validated.

### Reproduction of the published ccECP atomic energies

As an independent executable check, PySCF was used to read the generated
`GENBAS` and `ECPDATA` data and reproduce the UKS/PBE/cc-pV5Z atomic energies
reported in `recipes/<El>/ccECP/energies.txt`.

| El | Z | 2S | nao | ours (Ha) | published (Ha) | difference |
|---|---|---|---|---|---|---|
| H | 1 | 1 | 55 | -0.499956 | -0.49996 | +0.000004 |
| He | 2 | 0 | 55 | -2.892573 | -2.89257 | -0.000003 |
| Li | 3 | 1 | 79 | -0.201445 | -0.20144 | -0.000005 |
| Be | 4 | 0 | 90 | -0.995244 | -0.99524 | -0.000004 |
| B | 5 | 1 | 90 | -2.602589 | -2.60258 | -0.000009 |
| C | 6 | 2 | 90 | -5.409826 | -5.40983 | +0.000004 |
| N | 7 | 3 | 90 | -9.756575 | -9.75658 | +0.000005 |
| O | 8 | 2 | 90 | -15.884596 | -15.88454 | -0.000056 |
| F | 9 | 1 | 90 | -24.194576 | -24.19457 | -0.000006 |
| Ne | 10 | 0 | 90 | -35.023754 | -35.02375 | -0.000004 |
| Na | 11 | 1 | 90 | -0.191616 | -0.19162 | +0.000004 |
| Mg | 12 | 0 | 90 | -0.823519 | -0.82352 | +0.000001 |
| Al | 13 | 1 | 90 | -1.940214 | -1.94021 | -0.000004 |
| Si | 14 | 2 | 90 | -3.761088 | -3.76109 | +0.000002 |
| P | 15 | 3 | 90 | -6.448556 | -6.44856 | +0.000004 |
| S | 16 | 2 | 90 | -10.080050 | -10.08005 | -0.000000 |
| Cl | 17 | 1 | 90 | -14.895447 | -14.89545 | +0.000003 |
| Ar | 18 | 0 | 90 | -21.018323 | -21.01833 | +0.000007 |
| K | 19 | 1 | 103 | -28.198240 | -28.19824 | +0.000000 |
| Ca | 20 | 0 | 103 | -36.660000 | -36.66000 | -0.000000 |
| Sc | 21 | 1 | 148 | -46.508241 | -46.50819 | -0.000051 |
| Ti | 22 | 2 | 148 | -58.061315 | -58.06135 | +0.000035 |
| V | 23 | 3 | 148 | -71.417668 | -71.41763 | -0.000038 |
| Cr | 24 | 6 | 148 | -86.658635 | -86.65864 | +0.000005 |
| Mn | 25 | 5 | 148 | -103.884620 | -103.88462 | +0.000000 |
| Fe | 26 | 4 | 148 | -123.407373 | -123.40717 | -0.000203 |
| Co | 27 | 3 | 148 | -145.197658 | -145.19741 | -0.000248 |
| Ni | 28 | 2 | 148 | -169.457379 | -169.45740 | +0.000021 |
| Cu | 29 | 1 | 148 | -196.469980 | -196.46999 | +0.000010 |
| Zn | 30 | 0 | 148 | -226.433855 | -226.43387 | +0.000015 |
| Ga | 31 | 1 | 90 | -2.050668 | -2.05067 | +0.000002 |
| Ge | 32 | 2 | 90 | -3.758982 | -3.75899 | +0.000008 |
| As | 33 | 3 | 90 | -6.185996 | -6.18600 | +0.000004 |
| Se | 34 | 2 | 90 | -9.319983 | -9.31999 | +0.000007 |
| Br | 35 | 1 | 90 | -13.341596 | -13.34161 | +0.000014 |
| Kr | 36 | 0 | 90 | -18.486193 | -18.48621 | +0.000017 |
| Mo | 42 | 6 | 148 | -68.020779 | -68.02078 | +0.000001 |
| Pd | 46 | 0 | 148 | -127.338204 | -127.33827 | +0.000066 |
| Ag | 47 | 1 | 148 | -146.934262 | -146.93439 | +0.000128 |
| Te | 52 | 2 | 90 | -8.184423 | -8.18443 | +0.000007 |
| I | 53 | 1 | 90 | -11.442870 | -11.44288 | +0.000010 |
| W | 74 | 4 | 148 | -67.345790 | -67.34571 | -0.000080 |
| Ir | 77 | 3 | 148 | -104.497338 | -104.49734 | +0.000002 |
| Au | 79 | 1 | 148 | -135.962162 | -135.96224 | +0.000078 |
| Bi | 83 | 3 | 90 | -5.466003 | -5.46602 | +0.000017 |

Across 45 atoms, the largest absolute deviation is **0.000248 Ha** (Co).

## 11. Files

All paths below are relative to the output directory. Pass `--out DIR` to
choose it explicitly; otherwise the converter selects a safe `ccECP_cfour_output`
directory outside the source repository.

```
GENBAS                 CFOUR basis-set library (all 841 entries)
ECPDATA                CFOUR ECP library (all 65 entries)
GENBAS.mrcc            basis data with ECP entries appended for MRCC;
                       copy/rename this file to GENBAS in an MRCC job
GENBAS.RI              auxiliary (RI) fitting sets only (all 841 entries)
GENBAS.withRI          orbital + auxiliary sets, for CFOUR
GENBAS.mrcc.withRI     orbital + auxiliary + ECP in one file; copy/rename
                       this to GENBAS for an MRCC job using dfbasis_cor
ri_name_map.csv        element / source file / library name / CFOUR name
                       for every auxiliary set
ECPDATA.maxprec        ECPs taken from the most precise rendering of each
                       element rather than from a single uniform one. For
                       the 3d transition metals the .molpro ECP rendering
                       carries up to seven more significant digits than
                       .gamess, and for Ag, Au, Bi, I, Ir, Pd and Te the
                       native <El>.ccECP text does; this file uses
                       whichever is most precise, after verifying that all
                       renderings agree.
*.upper                the same four files with upper-cased element
                       symbols, for CFOUR builds that need that
inventory.csv/.json    what was found in the library
basis_name_map.csv     element / source file / library name / CFOUR name
ecp_name_map.csv       element / source file / library name / CFOUR name
naming_audit.csv       per-entry naming audit
validation_report.csv  one row per validation check
anomalies.csv          every non-PASS check plus source-data and
                       target-program compatibility notes
summary.json           machine-readable run summary
conversion.log         full run log
by_element/<El>/       per-element GENBAS, ECPDATA and summary.json,
                       plus one file per logical basis/ECP object using the
                       source object name with a `.cfour` suffix:
                         Fe.ccECP.cfour      (one ECPDATA entry)
                         Fe.cc-pVQZ.cfour    (one GENBAS entry)
                         Fe.aug-cc-pCV5Z.cfour  ...
                       mirroring Fe.ccECP.gamess / Fe.cc-pVQZ.gamess.
                       MRCC uses the same two formats, so these files
                       serve both programs.
```

## 12. Example CFOUR input (`ZMAT`)

Copy `GENBAS` and `ECPDATA` into the CFOUR working directory.

```
CO, CCSD(T)/cc-pVTZ with canonical ccECPs
C
O 1 R

R=1.1280

*CFOUR(CALC=CCSD(T),BASIS=SPECIAL,ECP=ON,SPHERICAL=ON
REF=RHF,SCF_CONV=9,CC_CONV=9,ABCDTYPE=AOBASIS)

C:cc-pVTZ-ccECP
O:cc-pVTZ-ccECP

C:ccECP
O:ccECP

```

Basis labels carry the `-ccECP` suffix (section 6); ECP labels do not.

`BASIS=SPECIAL` makes CFOUR read the per-atom basis labels from the first
block and, because `ECP=ON`, the per-atom ECP labels from the second.
Use `<Element>:NONE` in the ECP block for an atom that should carry no ECP.

The labels are spelled exactly as in `GENBAS`/`ECPDATA`, with periodic-table
element spelling in the standard files. If your CFOUR build folds atom labels
to upper case, use `GENBAS.upper` and `ECPDATA.upper` with the corresponding
upper-case element labels.

## 13. Example MRCC input (`MINP`)

MRCC reads user basis sets *and* user ECPs from a single `GENBAS` file in
the working directory, so use the generated `GENBAS.mrcc`:

```
cp GENBAS.mrcc ./GENBAS
```

```
basis=special
cc-pVTZ-ccECP
cc-pVTZ-ccECP

ecp=special
ccECP
ccECP

calc=CCSD(T)
mem=4GB
scftype=RHF
gauss=spher
unit=angs
geom=xyz
2

C    0.0000   0.0000   0.0000
O    0.0000   0.0000   1.1280
```

The `basis=special` / `ecp=special` lists must follow the atom order of
the geometry, and the ECP labels must match the `GENBAS` labels exactly.

The generated files use the MRCC fixed-column `NCORE`/`LMAX` record described
in section 9. For canonical ccECP H and He, read the `NCORE = 0` note before
using MRCC 25.1.1 or 26.1.1.

## 14. Auxiliary (RI) basis sets

The Pseudopotential Library ships an ORCA AutoAux `/C` auxiliary set beside
most orbital basis sets, as `<El>.<basis>.AutoAuxC.orca`. These are
density-fitting sets for the correlation step. All **841** are converted:
**833** taken from the library, and **8** generated locally with ORCA AutoAux
for the orbital sets the library does not cover (`B` and `C`
`cc-pV{D,T,Q,5}Z`, whose Gaussian-style shell headers make ORCA hang while
reading the basis, which is the likely reason they are absent upstream). The
generated ones are marked `source_origin=extra` in `ri_name_map.csv` and are
labelled as locally generated in their entry comment.

Auxiliary sets travel through the same parser body, block builder and renderer
as the orbital sets, so there is no second numerical code path. Only the source
wrapper differs: two lines at the top (`# ...` and `NewAuxCGTO <El>`) and a
terminating `end;`.

### Regenerating the auxiliary sets

`generate_orca_autoauxc.py` produces them from scratch. Run it from the
library's `recipes` directory:

```
python3 generate_orca_autoauxc.py --orca /path/to/orca
```

It walks every `recipes/<El>/<recipe>` directory whose recipe name contains
`ccECP`, identifies the ECP GAMESS file from its `Element-name GEN ncore lmax`
header, treats the remaining `*.gamess` files as orbital bases, and runs ORCA
with `! NoIter PrintBasis` purely as a basis parser and AutoAux generator. No
SCF or MP2 step is involved, so the multiplicity it writes is a formal
singlet/doublet and no ground-state configuration is implied; AutoAux derives
the fitting set from the orbital basis, not from a density.

Checked against the library: regenerating N `cc-pV{D,T,Q,5}Z` and
`aug-cc-pV{D,T}Z` reproduces the shipped `.AutoAuxC.orca` files
byte-for-byte. The 8 locally generated B and C sets are likewise reproduced
byte-for-byte by this script, even though it chooses a different multiplicity
than the run that first produced them, which confirms the multiplicity has no
effect here.

`autoauxc_generated/` holds all **841** auxiliary sources in one place so the
conversion can be driven by `--extra-aux` alone: 833 copied from the library
and the 8 generated locally. `autoauxc_generated/MANIFEST.csv` records the
origin, composition and a SHA-256 prefix for each, so which file came from
where stays explicit. The library's other `ccECP_*` recipe families carry a
further 232 AutoAuxC files; those are outside the canonical-ccECP scope of
this conversion and are not included.

Every AutoAux function is a single uncontracted primitive with coefficient
exactly 1, so each angular momentum becomes an identity coefficient matrix.
That is how MRCC stores its own `*-RI` entries, so the layout needs no special
handling. Angular momentum reaches `k` (l = 7) in the heavier sets, which is
within what MRCC accepts; its own `Ag:cc-pV5Z-PP-RI` reaches `l`.

Orbital sets without an auxiliary partner are recorded, never reconstructed:
any such entry appears in `anomalies.csv` with `check=aux_present`. With the
8 generated sets supplied, there are **0** in this output.

### Use in MRCC

```
cp GENBAS.mrcc.withRI ./GENBAS
```

```
basis=atomtype
Fe:cc-pVTZ-ccECP
dfbasis_scf=none
dfbasis_cor=cc-pVTZ-RI-ccECP
ecp=atomtype
Fe:ccECP
calc=MP2
core=0
mult=5
scftype=UHF
```

`dfbasis_cor` takes the label **without** the element prefix, unlike the
`basis`/`ecp` blocks under `atomtype`.

### Executable validation

MP2 was run twice per atom with the same orbital basis and ECP, once with
exact four-index integrals and once density-fitted with the converted set.
MRCC 26.1.1, cc-pVTZ:

| El | MP2 correlation, exact | MP2 correlation, DF | difference |
|---|---:|---:|---:|
| C | -0.06768956 | -0.06768933 | +2.4e-07 |
| Ne | -0.25115075 | -0.25115408 | -3.3e-06 |
| Si | -0.05843479 | -0.05844050 | -5.7e-06 |
| Fe | -0.62370956 | -0.62374618 | -3.7e-05 |
| Cu | -0.90836095 | -0.90842748 | -6.7e-05 |

MRCC echoes the resolved set, for example
`Fe cc-pvtz-ri-ccecp [ 15s 13p 12d 11f 10g 4h 3i ]`, which matches the
`composition` column of `ri_name_map.csv` exactly.

## 15. Troubleshooting and reproducibility

* **Library not found:** pass `--library /path/to/pseudopotentiallibrary` or set
  `PSEUDOPOTENTIAL_LIBRARY`. The option may also point directly to `recipes/`.
* **`basis_set_exchange` unavailable:** install it with `python -m pip install
  basis_set_exchange`, or use `--no-bse` if only the built-in parser/round-trip
  checks are required.
* **Basis-label collision:** keep the default `-ccECP` suffix, or use
  `--basis-suffix ""` only after checking the target basis library.
* **CFOUR cannot find a mixed-case label:** try the corresponding `*.upper`
  library file. Only the element symbol changes; the basis/ECP label text after
  the colon is unchanged.
* **Very high angular momentum:** `cc-pV6Z`/`cc-pCV6Z` may contain `i` functions.
  Confirm that the target CFOUR/MRCC build supports them.
* **Reproducibility:** retain `summary.json`, `conversion.log`, the mapping CSVs,
  and the Pseudopotential Library commit hash. These identify exactly which
  source data and rendering were used for every entry.

## 16. Citation and source data

The converter does not redefine or refit ccECP data. Scientific credit belongs
to the original ccECP/Pseudopotential Library sources. Each recipe's
`author.txt` is carried into generated comments and should be consulted for the
appropriate literature citation. Source-data observations and non-fatal consistency notes are reported in
`anomalies.csv`. Numerical source values are preserved exactly rather than
silently modified.

