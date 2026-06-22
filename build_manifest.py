#!/usr/bin/env python3
"""build_manifest.py — generate the production manifest (4 rows / MOF).

Per MOF: one He void row + three GCMC state points (100 bar/77 K, 5 bar/77 K,
5 bar/160 K). Reads cells from cells.py per CIF. stdlib only.

  python3 build_manifest.py [cif_list.txt] [-o manifest.csv]

cif_list.txt: one CIF filename (or path) per line; defaults to every *.cif in
cifs/. The GCMC model/forcefield are the two constants below — change MODEL/FF
once if the benchmark picks single-site (gcmc_lj/Generic_FH or Generic).
"""
import argparse
import csv
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
CIFS = os.path.join(ROOT, "cifs")
PY = sys.executable
CELLS_PY = os.path.join(ROOT, "cells.py")

# --- model selection (one-line change after the benchmark) -----------------
MODEL = "gcmc_dl"      # gcmc_dl (Darkrim-Levesque, expected winner) or gcmc_lj
FF = "Generic_FH"      # Generic_FH (Feynman-Hibbs) or Generic
VOID_FF = "Generic"    # forcefield for the helium void run

NCYC, NINIT = "10000", "5000"
# (jobtype, template, ff, T, P, tag) per GCMC state point
STATE_POINTS = [
    ("gcmc", MODEL, FF, "77",  "1.0e7", "100bar_77"),   # 100 bar storage (PS+TPS)
    ("gcmc", MODEL, FF, "77",  "5.0e5", "5bar_77"),     # 5 bar PS discharge
    ("gcmc", MODEL, FF, "160", "5.0e5", "5bar_160"),    # 5 bar TPS discharge
]

FIELDS = ["row_id", "mof", "cif", "jobtype", "template", "ff", "cells",
          "void_fraction", "temperature", "pressure", "ncyc", "ninit", "status",
          "uptake_abs", "uptake_exc", "uptake_err", "wall_s", "note"]


def cells_for(cif_path):
    """Call cells.py to size the supercell. ASSUMED CONTRACT: prints three
    ints '<nx> <ny> <nz>' on stdout. If cells.py differs, change this one call."""
    out = subprocess.run([PY, CELLS_PY, cif_path], capture_output=True, text=True)
    toks = out.stdout.split()
    ints = [t for t in toks if t.lstrip("-").isdigit()]
    if len(ints) < 3:
        raise RuntimeError(f"cells.py gave no '<nx> <ny> <nz>' for {cif_path}: "
                           f"{out.stdout!r} {out.stderr!r}")
    return " ".join(ints[:3])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cif_list", nargs="?", help="file with one CIF per line (default: all in cifs/)")
    ap.add_argument("-o", "--out", default=os.path.join(ROOT, "manifest.csv"))
    args = ap.parse_args()

    if args.cif_list:
        cifs = [l.strip() for l in open(args.cif_list) if l.strip()]
    else:
        cifs = sorted(f for f in os.listdir(CIFS) if f.endswith(".cif"))
    if not cifs:
        sys.exit("no CIFs found")

    rows = []
    for cif in cifs:
        cif_path = cif if os.path.isfile(cif) else os.path.join(CIFS, cif)
        mof = os.path.splitext(os.path.basename(cif))[0]
        try:
            cells = cells_for(cif_path)
        except Exception as e:                       # noqa: BLE001
            sys.stderr.write(f"skip {mof}: {e}\n")
            continue
        base = dict.fromkeys(FIELDS, "")
        base.update(cif=os.path.basename(cif), mof=mof, cells=cells,
                    ncyc=NCYC, ninit=NINIT, status="pending")
        # void row first
        v = dict(base); v.update(row_id=f"{mof}__void", jobtype="void",
                                 template="void.input", ff=VOID_FF, temperature="298")
        rows.append(v)
        for jt, tpl, ff, T, P, tag in STATE_POINTS:
            g = dict(base); g.update(row_id=f"{mof}__{tag}", jobtype=jt, template=tpl,
                                     ff=ff, temperature=T, pressure=P)
            rows.append(g)

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"wrote {len(rows)} rows for {len(cifs)} MOFs -> {args.out} (model={MODEL}, ff={FF})")


if __name__ == "__main__":
    main()
