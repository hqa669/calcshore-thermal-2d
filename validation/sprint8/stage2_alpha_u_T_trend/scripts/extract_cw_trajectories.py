#!/usr/bin/env python3
"""§3.3 — Extract CW T_core(t) trajectories for all 32 Stage 2-alpha-u-T-trend datasets.

Reads output.txt from each of the 32 cw_data/ folders.
Extracts T_core_CW at (wi=0, di=24) = centerline midpoint (12.19 m depth, 6.1 m width).

Saves data/cw_trajectories.npz with:
  t_hrs          (N_t,)      — common time grid in hours
  T_pc_F         (32,)       — placement temperature in °F
  alpha_u        (32,)       — nominal α_u (0.20/0.40/0.60/0.80)
  folder_names   (32,)       — string labels
  T_core_CW_C    (32, N_t)   — T_core in °C, row i matches T_pc_F[i]/alpha_u[i]
  T_core_CW_F    (32, N_t)   — T_core in °F (for easy residual comparison)

Sorted by (alpha_u, T_pc_F) — ascending in both.
"""
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[1]
REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))

from cw_scenario_loader import parse_cw_temp_output

CW_DATA = HERE / "cw_data"
DATA    = HERE / "data"

ALPHA_TAG_MAP = {"02": 0.20, "04": 0.40, "06": 0.60, "08": 0.80}
WI_CENTER = 0   # wi=0 = 6.1 m (centerline for 40 ft = 12.2 m half-mat)
DI_CENTER = 24  # di=24 = 12.19 m (mid-depth for 80 ft slab)


def main():
    DATA.mkdir(parents=True, exist_ok=True)

    folders = sorted(CW_DATA.glob("thermal_temp_alpha*"))
    assert len(folders) == 32, f"Expected 32 folders, found {len(folders)}"

    # Parse folder metadata
    records = []
    for folder in folders:
        m = re.match(r"thermal_temp_alpha(\d{2})_(\d+)$", folder.name)
        assert m, f"Unexpected folder name: {folder.name}"
        alpha_nom = ALPHA_TAG_MAP[m.group(1)]
        T_pc      = float(m.group(2))
        records.append((alpha_nom, T_pc, folder))

    # Sort by (alpha_u, T_pc)
    records.sort(key=lambda r: (r[0], r[1]))

    print(f"§3.3 CW trajectory extraction — {len(records)} datasets")
    print(f"  Anchor: wi={WI_CENTER} ({6.1:.1f} m from corner = centerline), "
          f"di={DI_CENTER} (12.19 m depth = mid-slab)")
    print()

    ref_t_hrs = None
    T_core_rows_F = []

    for alpha_nom, T_pc, folder in records:
        out_path = folder / "output.txt"
        assert out_path.exists(), f"output.txt missing: {folder.name}"
        v = parse_cw_temp_output(str(out_path))

        nD, nW = v.T_field_F.shape[1], v.T_field_F.shape[2]
        assert nD > DI_CENTER, f"{folder.name}: nD={nD} < {DI_CENTER+1}"
        assert nW > WI_CENTER, f"{folder.name}: nW={nW} < {WI_CENTER+1}"

        # Verify / set common time grid
        if ref_t_hrs is None:
            ref_t_hrs = v.time_hrs
        else:
            if not np.allclose(v.time_hrs, ref_t_hrs, atol=1e-6):
                print(f"  WARNING: {folder.name} has different time grid "
                      f"(n_t={len(v.time_hrs)} vs ref {len(ref_t_hrs)})")
                # Interpolate to common grid
                T_core_F_interp = np.interp(ref_t_hrs, v.time_hrs,
                                             v.T_field_F[:, DI_CENTER, WI_CENTER])
                T_core_rows_F.append(T_core_F_interp)
                print(f"  {folder.name}  α_u={alpha_nom:.2f}  T_pc={T_pc:.0f}°F  "
                      f"T_core[t=168]={T_core_F_interp[-1]:.2f}°F  [interpolated]")
                continue

        T_core_F = v.T_field_F[:, DI_CENTER, WI_CENTER]
        T_core_rows_F.append(T_core_F.copy())
        print(f"  {folder.name:<35}  α_u={alpha_nom:.2f}  T_pc={T_pc:>5.0f}°F  "
              f"T_core_t0={T_core_F[0]:.2f}°F  T_core_t168={T_core_F[-1]:.2f}°F")

    T_core_CW_F = np.array(T_core_rows_F)          # (32, N_t)
    T_core_CW_C = (T_core_CW_F - 32.0) * 5.0 / 9.0

    alpha_u_arr  = np.array([r[0] for r in records])
    T_pc_F_arr   = np.array([r[1] for r in records])
    folder_names = np.array([r[2].name for r in records])

    out_npz = DATA / "cw_trajectories.npz"
    np.savez(out_npz,
             t_hrs=ref_t_hrs,
             T_pc_F=T_pc_F_arr,
             alpha_u=alpha_u_arr,
             folder_names=folder_names,
             T_core_CW_C=T_core_CW_C,
             T_core_CW_F=T_core_CW_F)
    print(f"\nWrote {out_npz}")
    print(f"  t_hrs shape: {ref_t_hrs.shape}")
    print(f"  T_core_CW_C shape: {T_core_CW_C.shape}")
    print(f"  T_pc_F range: {T_pc_F_arr.min():.0f}–{T_pc_F_arr.max():.0f}°F")
    print(f"  alpha_u values: {sorted(set(alpha_u_arr))}")
    print("Done.")


if __name__ == "__main__":
    main()
