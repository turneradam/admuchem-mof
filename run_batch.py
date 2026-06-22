#!/usr/bin/env python3
"""run_batch.py — RASPA2 driver for the MOF H2 campaign (void + GCMC).

Lives in ~/software/mof/, version-controlled, so a fix deploys via /pull + /run
with no watcher restart. Plain procedural python, stdlib + numpy.

Two jobtypes share one manifest:
  void  -> run templates/void.input (helium Widom), parse the void fraction,
           cache it in voids.csv keyed by `mof`.
  gcmc  -> fill a GCMC template, run, parse absolute+excess loading. The void
           fraction comes from the row, else voids.csv[mof], else the
           HeliumVoidFraction line is stripped (the config-R baseline).

Manifest schema (manifest.csv / benchmark.csv):
  row_id,mof,cif,jobtype,template,ff,cells,void_fraction,temperature,pressure,
  ncyc,ninit,status,uptake_abs,uptake_exc,uptake_err,wall_s,note
  status in {pending,running,done,failed,skipped}

Ordering is two-phase: ALL pending void rows (so voids.csv is fully populated)
then ALL pending gcmc rows, each phase by ascending cell multiplier (cheap MOFs
bank first, the 64x tail runs last). N=12 concurrent serial RASPA jobs. The pool
is threads owning `simulate` subprocesses, so a SIGTERM tears every worker down.

Success = RASPA's final-averages block parsed to finite value(s); the stdout
completion banner is NOT used as a gate (its wording varies by build). A `failed`
row that has in fact completed is re-read and salvaged to `done` on the next run
rather than recomputed. Atomic manifest/voids writes. Exit 0 if the queue ran
(failures included); nonzero only if it cannot run at all.
"""
import argparse
import csv
import math
import os
import re
import shutil
import signal
import statistics
import subprocess
import sys
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor

# --- configuration ---------------------------------------------------------
N = int(os.environ.get("RASPA_N", "12"))            # worker pool; benchmark may retune
ROOT = os.path.dirname(os.path.abspath(__file__))
CIFS = os.path.join(ROOT, "cifs")
BENCH = os.path.join(ROOT, "benchmark")
JOBS = os.path.join(ROOT, "jobs")
TEMPLATES = os.path.join(ROOT, "templates")
VOIDS = os.path.join(ROOT, "voids.csv")
RASPA_DIR = os.environ.get("RASPA_DIR", os.path.expanduser("~/software/raspa"))
RASPA_BIN = os.path.join(RASPA_DIR, "bin", "simulate")

FIELDS = ["row_id", "mof", "cif", "jobtype", "template", "ff", "cells",
          "void_fraction", "temperature", "pressure", "ncyc", "ninit", "status",
          "uptake_abs", "uptake_exc", "uptake_err", "wall_s", "note"]
DEFAULT_NCYC, DEFAULT_NINIT = "10000", "5000"

# RASPA output labels (verify against your build's Output/System_0/*.data).
RE_ABS_GRAV = re.compile(r"Average loading absolute \[mol/kg framework\]\s+([-+0-9.eE]+)\s*(?:\+/-\s*([-+0-9.eE]+))?")
RE_EXC_GRAV = re.compile(r"Average loading excess \[mol/kg framework\]\s+([-+0-9.eE]+)")
RE_ABS_VOL = re.compile(r"Average loading absolute \[cm\^3 \(STP\)/cm\^3 framework\]\s+([-+0-9.eE]+)")
# void.input is upstream-authored — try the common forms; FIRST thing to verify.
RE_VOID = [
    re.compile(r"[Vv]oid[_ ]fraction[^0-9-]*([-+0-9.eE]+)"),
    re.compile(r"[Hh]elium void fraction[^0-9-]*([-+0-9.eE]+)"),
    re.compile(r"[Ww]idom Rosenbluth[- ]weight[^0-9-]*([-+0-9.eE]+)"),
    re.compile(r"Average Widom Rosenbluth factor[^0-9-]*([-+0-9.eE]+)"),
]

# --- shared state ----------------------------------------------------------
_mlock = threading.Lock()
_vlock = threading.Lock()
_live_lock = threading.Lock()
_live = set()
_stop = threading.Event()


def _fail(msg, code=1):
    sys.stderr.write(msg.rstrip() + "\n")
    sys.exit(code)


def _atomic_write_csv(path, fields, rows):
    d = os.path.dirname(path) or "."
    fd, tmp = tempfile.mkstemp(dir=d, prefix=".tmp.", suffix=".csv")
    try:
        with os.fdopen(fd, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def read_manifest(path):
    if not os.path.isfile(path):
        _fail(f"no manifest at {path}")
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        for k in FIELDS:
            r.setdefault(k, "")
        if not r.get("ncyc"):
            r["ncyc"] = DEFAULT_NCYC
        if not r.get("ninit"):
            r["ninit"] = DEFAULT_NINIT
    return rows


def write_manifest(path, rows):
    _atomic_write_csv(path, FIELDS, rows)


def update_row(path, row_id, **changes):
    with _mlock:
        rows = read_manifest(path)
        for r in rows:
            if r["row_id"] == row_id:
                r.update({k: ("" if v is None else str(v)) for k, v in changes.items()})
                break
        write_manifest(path, rows)


def read_voids():
    d = {}
    if os.path.isfile(VOIDS):
        with open(VOIDS, newline="") as f:
            for r in csv.DictReader(f):
                d[r["mof"]] = r.get("void_fraction", "")
    return d


def write_void(mof, vf):
    with _vlock:
        d = read_voids()
        d[mof] = f"{vf:.6g}"
        rows = [{"mof": k, "void_fraction": v} for k, v in sorted(d.items())]
        _atomic_write_csv(VOIDS, ["mof", "void_fraction"], rows)


def on_signal(signum, _frame):
    _stop.set()
    with _live_lock:
        for p in list(_live):
            try:
                p.terminate()
            except Exception:
                pass


def cell_mult(row):
    try:
        return math.prod(int(x) for x in row.get("cells", "").split())
    except (ValueError, TypeError):
        return 1


def resolve_cif(cif):
    cands = [cif, os.path.join(ROOT, cif), os.path.join(CIFS, cif)]
    for c in cands:
        if os.path.isfile(c):
            return c
    if os.path.isdir(BENCH):
        for sub in os.listdir(BENCH):
            p = os.path.join(BENCH, sub, cif)
            if os.path.isfile(p):
                return p
    return None


def template_path(name):
    fn = name if name.endswith(".input") else name + ".input"
    return os.path.join(TEMPLATES, fn)


def render_job(row, vf):
    """Build jobs/<row_id>/ from the row's template + CIF. vf is the void
    fraction to inject, or None to strip the HeliumVoidFraction line."""
    rid = row["row_id"]
    jobdir = os.path.join(JOBS, rid)
    shutil.rmtree(jobdir, ignore_errors=True)
    os.makedirs(jobdir)

    cif = resolve_cif(row["cif"])
    if cif is None:
        raise FileNotFoundError(f"cif not found: {row['cif']}")
    framework = os.path.splitext(os.path.basename(cif))[0]
    shutil.copy(cif, os.path.join(jobdir, framework + ".cif"))

    tpl = template_path(row["template"])
    if not os.path.isfile(tpl):
        raise FileNotFoundError(f"template not found: {tpl}")
    with open(tpl) as fh:
        text = fh.read()

    subs = {
        "__NCYC__": row.get("ncyc") or DEFAULT_NCYC,
        "__NINIT__": row.get("ninit") or DEFAULT_NINIT,
        "__FF__": row.get("ff", "Generic"),
        "__MOF__": framework,
        "__CELLS__": row.get("cells", "1 1 1"),
        "__T__": row.get("temperature", ""),
        "__P__": row.get("pressure", ""),
        "__VF__": "" if vf is None else f"{float(vf):.6g}",
    }
    for k, v in subs.items():
        text = text.replace(k, str(v))
    if vf is None:
        text = "\n".join(l for l in text.splitlines() if "HeliumVoidFraction" not in l) + "\n"
    with open(os.path.join(jobdir, "simulation.input"), "w") as fh:
        fh.write(text)
    return jobdir


def _run_simulate(jobdir):
    p = None
    try:
        with open(os.path.join(jobdir, "raspa.log"), "w") as logf:
            p = subprocess.Popen([RASPA_BIN, "-i", "simulation.input"],
                                 cwd=jobdir, stdout=logf, stderr=subprocess.STDOUT)
            with _live_lock:
                _live.add(p)
            p.wait()
    finally:
        if p is not None:
            with _live_lock:
                _live.discard(p)


def _read_output(jobdir):
    """RASPA's stdout (raspa.log) plus every Output/System_0/*.data file, as one
    blob. Completion is judged downstream by whether the final-averages block
    actually parsed -- not by a stdout 'Simulation finished' banner. That banner's
    exact wording varies between RASPA builds, so gating on it silently failed
    runs that had in fact completed and written their results to disk."""
    blob = ""
    try:
        blob += open(os.path.join(jobdir, "raspa.log"), errors="ignore").read()
    except OSError:
        pass
    outdir = os.path.join(jobdir, "Output", "System_0")
    if os.path.isdir(outdir):
        for fn in sorted(os.listdir(outdir)):
            if fn.endswith(".data"):
                try:
                    blob += "\n" + open(os.path.join(outdir, fn), errors="ignore").read()
                except OSError:
                    pass
    return blob


def parse_gcmc(jobdir):
    blob = _read_output(jobdir)
    if not blob.strip():
        return None, "no RASPA output (raspa.log + Output/System_0 both empty)"
    m_abs = RE_ABS_GRAV.search(blob)
    m_exc = RE_EXC_GRAV.search(blob)
    m_vol = RE_ABS_VOL.search(blob)
    if not (m_abs and m_exc and m_vol):
        tail = " ".join(blob.split())[-180:]
        return None, f"no final-loading block -- RASPA did not complete (tail: ...{tail})"
    try:
        abs_g = float(m_abs.group(1))
        err = float(m_abs.group(2)) if m_abs.group(2) else 0.0
        exc_g = float(m_exc.group(1))
        vol = float(m_vol.group(1))
    except ValueError:
        return None, "loading value not numeric"
    for v in (abs_g, exc_g, vol):
        if not math.isfinite(v):
            return None, "non-finite loading"
    return {"abs": abs_g, "exc": exc_g, "err": err}, ""


def parse_void(jobdir):
    blob = _read_output(jobdir)
    if not blob.strip():
        return None, "no RASPA output (raspa.log + Output/System_0 both empty)"
    for rx in RE_VOID:
        m = rx.search(blob)
        if m:
            try:
                vf = float(m.group(1))
            except ValueError:
                continue
            if math.isfinite(vf):
                return vf, ""
    return None, "void fraction not parsed (check void.input output labels)"


def worker_void(row, manifest):
    if _stop.is_set():
        return
    rid = row["row_id"]
    try:
        jobdir = render_job(row, None)
    except Exception as e:                          # noqa: BLE001
        update_row(manifest, rid, status="failed", note=f"setup: {e}")
        return
    update_row(manifest, rid, status="running", note="")
    t0 = time.time()
    _run_simulate(jobdir)
    wall = round(time.time() - t0, 1)
    if _stop.is_set():
        return
    vf, note = parse_void(jobdir)
    if vf is None:
        update_row(manifest, rid, status="failed", wall_s=wall, note=note)
        return
    write_void(row["mof"], vf)
    update_row(manifest, rid, status="done", void_fraction=f"{vf:.6g}", wall_s=wall, note="")


def worker_gcmc(row, manifest, voids):
    if _stop.is_set():
        return
    rid = row["row_id"]
    vf_raw = (row.get("void_fraction") or "").strip()
    if not vf_raw:
        vf_raw = (voids.get(row["mof"]) or "").strip()
    vf = vf_raw if vf_raw else None                 # None -> strip He line (config R)
    try:
        jobdir = render_job(row, vf)
    except Exception as e:                          # noqa: BLE001
        update_row(manifest, rid, status="failed", note=f"setup: {e}")
        return
    update_row(manifest, rid, status="running", note="")
    t0 = time.time()
    _run_simulate(jobdir)
    wall = round(time.time() - t0, 1)
    if _stop.is_set():
        return
    res, note = parse_gcmc(jobdir)
    if res is None:
        update_row(manifest, rid, status="failed", wall_s=wall, note=note)
        return
    update_row(manifest, rid, status="done",
               uptake_abs=f"{res['abs']:.6g}", uptake_exc=f"{res['exc']:.6g}",
               uptake_err=f"{res['err']:.6g}", wall_s=wall, note="")


def run_phase(rows, fn):
    if not rows:
        return
    with ThreadPoolExecutor(max_workers=N) as ex:
        futs = [ex.submit(fn, r) for r in rows]
        for fut in futs:
            try:
                fut.result()
            except Exception as e:                  # noqa: BLE001
                sys.stderr.write(f"worker error (non-fatal): {e}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="manifest.csv", help="manifest filename under the mof root")
    ap.add_argument("--only", metavar="ROW_ID", help="reset and run just this row")
    args = ap.parse_args()

    manifest = args.manifest if os.path.isabs(args.manifest) else os.path.join(ROOT, args.manifest)
    if not os.path.isfile(RASPA_BIN) or not os.access(RASPA_BIN, os.X_OK):
        _fail(f"raspa simulate not found/executable at {RASPA_BIN}")
    os.makedirs(JOBS, exist_ok=True)

    rows = read_manifest(manifest)
    ids = {r["row_id"] for r in rows}
    if args.only is not None and args.only not in ids:
        _fail(f"--only {args.only}: not in {os.path.basename(manifest)}")

    # salvage: rows an earlier build marked `failed` purely on the banner gate may
    # actually have completed and written their results. Re-read their existing
    # output (no recompute) and flip any that now parse to `done`, before we work
    # out what is still pending. Genuine failures (no job dir, or no final block)
    # are left as-is.
    salvaged = 0
    with _mlock:
        rows = read_manifest(manifest)
        for r in rows:
            if r["status"] != "failed":
                continue
            if args.only and r["row_id"] != args.only:
                continue
            jd = os.path.join(JOBS, r["row_id"])
            if not os.path.isdir(jd):
                continue
            if r["jobtype"] == "gcmc":
                res, _ = parse_gcmc(jd)
                if res is not None:
                    r.update(status="done", note="salvaged",
                             uptake_abs=f"{res['abs']:.6g}",
                             uptake_exc=f"{res['exc']:.6g}",
                             uptake_err=f"{res['err']:.6g}")
                    salvaged += 1
            elif r["jobtype"] == "void":
                vf, _ = parse_void(jd)
                if vf is not None:
                    r.update(status="done", note="salvaged",
                             void_fraction=f"{vf:.6g}")
                    salvaged += 1
        if salvaged:
            write_manifest(manifest, rows)
    if salvaged:
        print(f"salvaged {salvaged} previously-failed row(s) by re-reading existing output")

    # resume: stale `running` -> pending + wipe job dir; --only forces a *re-run*
    # of its row -- but not one salvage just marked done (don't bin good output).
    with _mlock:
        rows = read_manifest(manifest)
        for r in rows:
            forced = bool(args.only) and r["row_id"] == args.only and r["status"] != "done"
            if r["status"] == "running" or forced:
                r["status"] = "pending"
                r["note"] = ""
                shutil.rmtree(os.path.join(JOBS, r["row_id"]), ignore_errors=True)
        write_manifest(manifest, rows)

    def selectable(r):
        return r["status"] == "pending" and (args.only is None or r["row_id"] == args.only)

    voids_rows = sorted([r for r in rows if r["jobtype"] == "void" and selectable(r)], key=cell_mult)
    gcmc_rows = sorted([r for r in rows if r["jobtype"] == "gcmc" and selectable(r)], key=cell_mult)
    if not voids_rows and not gcmc_rows:
        print("nothing pending — queue already drained")
        return 0

    signal.signal(signal.SIGTERM, on_signal)
    signal.signal(signal.SIGINT, on_signal)

    print(f"{os.path.basename(manifest)}: {len(voids_rows)} void + {len(gcmc_rows)} gcmc pending, "
          f"N={N}, raspa={RASPA_BIN}")
    t0 = time.time()

    run_phase(voids_rows, lambda r: worker_void(r, manifest))
    voids = read_voids()
    run_phase(gcmc_rows, lambda r: worker_gcmc(r, manifest, voids))

    final = read_manifest(manifest)
    done = sum(1 for r in final if r["status"] == "done")
    failed = sum(1 for r in final if r["status"] == "failed")
    pend = sum(1 for r in final if r["status"] == "pending")
    print(f"finished in {round(time.time() - t0, 1)}s — done={done} failed={failed} pending={pend} "
          f"({'interrupted' if _stop.is_set() else 'queue drained'})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
