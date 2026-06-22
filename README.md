# mof-starter — seeds `~/software/mof` (the RASPA campaign repo)

Drop these into `~/software/mof/` on admuchem. They are the version-controlled
half of the RASPA engine (the bot/watcher/pp half lives in `~/bots`). `/pull`
updates both repos, so a driver fix ships without a restart.

## What's in here

| File | Goes to | Notes |
|---|---|---|
| `run_batch.py` | `~/software/mof/run_batch.py` | the driver |
| `build_manifest.py` | `~/software/mof/build_manifest.py` | production manifest generator |
| `templates/gcmc_lj.input` | `~/software/mof/templates/` | single-site H2 |
| `templates/gcmc_dl.input` | `~/software/mof/templates/` | Darkrim–Levesque + Ewald |
| `benchmark.csv.example` | copy to `benchmark.csv`, edit | the R/C1/C2/C3 ladder |
| `.gitignore` | `~/software/mof/.gitignore` | ignores jobs/, Output/, *.data, voids.csv, repo/ |

## Already on the box (do NOT overwrite from here)

`templates/void.input` (helium Widom) and `cells.py` (unit-cell sizer) exist
already. This package does not ship them.

## Manifest schema (both manifests)

```
row_id,mof,cif,jobtype,template,ff,cells,void_fraction,temperature,pressure,
ncyc,ninit,status,uptake_abs,uptake_exc,uptake_err,wall_s,note
```

- `jobtype` = `void` | `gcmc`
- `void` row → runs `void.input`, caches the He void fraction in `voids.csv[mof]`
- `gcmc` row → void fraction is taken from the row, else `voids.csv[mof]`, else
  the `HeliumVoidFraction` line is stripped (the config-R "reproduce the bug" case)
- `uptake_abs` / `uptake_exc` / `uptake_err` are **mol/kg** (gravimetric);
  volumetric is re-read from job outputs by `/pp deliverable`
- defaults: `ncyc=10000`, `ninit=5000`

## Production manifest

```
python3 build_manifest.py            # all *.cif in cifs/  -> manifest.csv
python3 build_manifest.py list.txt   # only the CIFs named in list.txt
```

Emits 4 rows per MOF: a void run + three GCMC state points (100 bar/77 K,
5 bar/77 K, 5 bar/160 K). The H2 model is the two constants `MODEL` / `FF` at the
top of `build_manifest.py` — change them once if the benchmark picks single-site.

**`cells.py` contract assumed:** `python3 cells.py <cif>` prints `<nx> <ny> <nz>`
on stdout. If yours differs, fix the one `cells_for()` call.

## Benchmark manifest

`benchmark.csv` is hand-authored (no generator). The example encodes the ladder
as `<mof>__<config>__p<tag>` row_ids — `/pp metrics` groups computed isotherms by
that `<config>` segment. The benchmark has **no void rows**, so a blank
`void_fraction` always means "strip" (config R); C1/C2/C3 carry an explicit void
fraction (replace the placeholder values with real ones).

## Deliverables (`/pp`)

- `/pp deliverable production` → PS = u(100 bar,77 K) − u(5 bar,77 K);
  TPS = u(100 bar,77 K) − u(5 bar,160 K); gravimetric + volumetric, absolute uptake.
- `/pp metrics benchmark ref/hkust1.csv` → computed-vs-reference isotherm overlay,
  per-config mean abs %err, parity PNG.
