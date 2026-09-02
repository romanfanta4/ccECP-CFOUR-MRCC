# ORCA conversion test

Checks that the `.orca` files produced by `../convert_gamess_ecp_to_orca.py`
reproduce the reference energies that ship with the library.

## Running

```
python3 test_orca_conversion.py
```

ORCA is taken from `~/codes/orca_6_1_0_avx2/orca`; override with `ORCA_EXE`.
Put the ORCA directory on `PATH` as well, since ORCA calls its own
sub-executables. Each case gets a subdirectory here holding its `t.inp` and
`t.out`; nothing is written into the recipes tree.

## What it does

For every case it builds an ORCA input that uses **only converted files**:

* `<El>.cc-pV5Z.orca` for the orbital basis, via `%basis GTOName = "..."`
* `<El>.ccECP.orca` for the potential, pasted into `%basis` as a `NewECP`
  block

The `NewECP` block has to be pasted rather than referenced: ORCA 6.1.0 rejects
both `%include` of the file and `ReadFragECP` of the file.

The total energy is then compared with the `Pyscf/PBE` value in that element's
`energies.txt`. Those references are UKS/PBE/cc-pV5Z (their NOTE 2), so the
same method and basis are used here and the numbers are directly comparable.

Settings: `! UKS PBE NORI DEFGRID3 VeryTightSCF PAL8`. `NORI` keeps the four-centre
integrals exact, so no auxiliary basis has to be invented for an ECP basis
that ships without one; the grid and convergence keywords push those error
sources below the comparison level.

## Tolerance

Agreement at the 1e-4 Ha level is the expected outcome, because the reference
is quoted to five decimals and the two programs use different DFT integration
grids. That is still a sharp test of the conversion:

* a misread `N_core`, or a channel written to the wrong angular momentum,
  moves the energy by whole Hartrees
* a misread contraction coefficient or exponent moves it by several mHa

so anything that would actually be wrong in the converted file is far outside
the tolerance.

The cases were picked from elements whose ccECP data had already been checked
independently against the same references with PySCF, covering light
main-group, heavy main-group, a 3d transition metal and a 5d transition metal.
