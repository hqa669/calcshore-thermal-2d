#!/usr/bin/env python3
"""Pre-flight smoke test for the Sprint 9 sweep wrapper.

Runs ONE combination (mix01, Hu_factor=0.95, c_multiplier=1.00) and:
  1. Compares scalars against PILOT_REPORT.md baselines
     (Hu_factor=0.951403 vs 0.95 → expect ≤ 0.05°F drift)
  2. Reports wall-clock time of the solve_hydration_2d call
"""
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from run_engine_sweep import run_one, load_cw_once, PILOT, C_ANCHOR

# Pilot baselines from PILOT_REPORT.md §1.3 / §1.6
PILOT_REF = {
    "max_R_F":   1.0508,
    "di_at_max": 47,
    "R_di47_F": -1.0508,   # max|R| sits at di=47 with negative sign per §1.6 dip
    "R_di48_F":  0.2055,   # §1.6
}
TOL_F = 0.05

print("=" * 70)
print("Pre-flight smoke test — MIX-01, Hu_factor=0.95, c_multiplier=1.00")
print(f"(Pilot used Hu_factor=0.951403; this grid point is 0.95 → expect ≤ {TOL_F}°F drift)")
print(f"c_eff = 1.00 × {C_ANCHOR} = {1.00 * C_ANCHOR:.4f}  (matches pilot)")
print("=" * 70)

cw_folder = PILOT / "cw_data/mix01"
cw = load_cw_once(cw_folder)

t0 = time.perf_counter()
scalars, log = run_one("mix01", cw_folder, Hu_fac=0.95, c_mult=1.00, cw=cw)
t1 = time.perf_counter()
elapsed = t1 - t0

print("\nResult scalars:")
for k, v in scalars.items():
    print(f"  {k:12s} = {v}")

print("\nWrapper log:")
for k, v in log.items():
    print(f"  {k:14s} = {v}")

print("\n--- Comparison against PILOT_REPORT.md ---")
print(f"{'metric':12s} | {'sweep':>10s} | {'pilot':>10s} | {'|delta|':>10s} | status")
print("-" * 60)
all_ok = True
for key, ref in PILOT_REF.items():
    sweep_val = scalars[key]
    if key == "di_at_max":
        ok = sweep_val == ref
        delta_str = f"{abs(sweep_val - ref):>10d}"
    else:
        delta = abs(sweep_val - ref)
        ok = delta <= TOL_F
        delta_str = f"{delta:>10.4f}"
    status = "OK" if ok else "FAIL"
    if not ok:
        all_ok = False
    print(f"{key:12s} | {sweep_val:>10} | {ref:>10} | {delta_str} | {status}")

print()
print(f"Wall-clock for solve_hydration_2d + extraction: {elapsed:.1f} s")
print(f"Estimated full-sweep runtime (80 × {elapsed:.1f} s): "
      f"{80 * elapsed / 60.0:.1f} min  ({80 * elapsed / 3600.0:.2f} hr)")

print()
if not all_ok:
    print("HALT: smoke test FAILED — scalars differ from pilot by more than tolerance.")
    print("Do NOT proceed to the full sweep until investigated.")
    sys.exit(1)

if elapsed > 90.0:
    print(f"WARN: per-run time {elapsed:.1f}s exceeds 90s threshold — full sweep will exceed plan estimate.")
    print(f"      Projected runtime: {80 * elapsed / 3600.0:.2f} hr")
    sys.exit(2)

print(f"PASS: all scalars within {TOL_F}°F of pilot; per-run time {elapsed:.1f}s ≤ 90s.")
print("OK to proceed to full 80-run sweep.")
sys.exit(0)
