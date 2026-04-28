"""
Generate a CW adiabatic reference CSV from a CW output.txt file.

Usage:
    python generate_adiabatic_reference.py <path_to_output.txt> [mix_label]

Example:
    python generate_adiabatic_reference.py ~/Downloads/HydrationCenter_mix02/output.txt mix02
    python generate_adiabatic_reference.py "C:\\Users\\Qinang\\Downloads\\HydrationCenter_mix04\\output.txt" mix04

If mix_label is not provided, it's inferred from the parent directory name.

Outputs (in the same directory as output.txt):
    cw_adiabatic_reference_<label>.csv
    adiabatic_diagnostic_<label>.csv  (debug — full spatial probes)

The CSV column 'T_center_F_adiabatic' is the kinetics ground truth — use
this for engine validation. See cw_adiabatic_reference_mix01_README.md
for the full methodology and CW configuration recipe.

REQUIREMENTS:
    - This script must be run from a directory where cw_scenario_loader.py
      is importable, OR the loader's directory must be on PYTHONPATH.
    - The CW output.txt must come from an adiabatic-isolation scenario:
      large geometry (>=200x200 ft), constant ambient = soil = placement
      temperature, sides shaded, 100% cloud, low wind, full curing.
"""
import argparse
import os
import sys
import numpy as np


def find_loader():
    """Try to import cw_scenario_loader from common locations."""
    # Try current directory first
    sys.path.insert(0, os.getcwd())
    try:
        from cw_scenario_loader import parse_cw_temp_output
        return parse_cw_temp_output
    except ImportError:
        pass

    # Try the script's own directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, script_dir)
    try:
        from cw_scenario_loader import parse_cw_temp_output
        return parse_cw_temp_output
    except ImportError:
        pass

    raise ImportError(
        "Cannot find cw_scenario_loader.py. Place this script in the same "
        "directory as cw_scenario_loader.py, or add that directory to "
        "PYTHONPATH."
    )


def find_centerline_axis(widths_m, depths_m, T_field_F):
    """
    CW exports widths_m sorted descending: widths_m[0] = form face,
    widths_m[-1] = centerline (0 m). But this convention has been
    observed to vary; let's empirically find which width index is
    actually the centerline by checking which column has the highest
    final temperature (= farthest from any boundary = centerline).
    """
    nt, nd, nw = T_field_F.shape

    # The centerline column should have the highest temperature at
    # the final timestep (deepest interior, no boundary cooling).
    # Check the mid-depth row.
    mid_d = nd // 2
    final_row = T_field_F[-1, mid_d, :]  # T at mid-depth, all widths, final time

    # Centerline is wherever the temperature is highest at the end
    centerline_w_idx = int(np.argmax(final_row))

    # Also figure out which depth row is the geometric center
    # (not necessarily the hottest — bottom can sometimes be hotter
    # depending on T_gw). Use the widthwise-hottest column and find
    # the hottest depth in it.
    final_col = T_field_F[-1, :, centerline_w_idx]
    mid_d_idx = int(np.argmax(final_col))

    return centerline_w_idx, mid_d_idx


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("output_txt",
                        help="Path to CW output.txt file")
    parser.add_argument("mix_label", nargs="?", default=None,
                        help="Mix label (e.g. 'mix02'). "
                             "Defaults to parent directory name.")
    parser.add_argument("--save-dir", default=None,
                        help="Where to save outputs (default: same dir as output.txt)")
    args = parser.parse_args()

    output_path = os.path.expanduser(args.output_txt)
    if not os.path.isfile(output_path):
        print(f"ERROR: file not found: {output_path}")
        sys.exit(1)

    scenario_dir = os.path.dirname(os.path.abspath(output_path))

    # Infer label
    if args.mix_label:
        label = args.mix_label
    else:
        # Use parent dir name, e.g. 'HydrationCenter_mix02' -> 'mix02'
        parent = os.path.basename(scenario_dir)
        if "_" in parent:
            label = parent.split("_")[-1]
        else:
            label = parent
    print(f"Mix label: {label}")

    save_dir = os.path.expanduser(args.save_dir) if args.save_dir else scenario_dir

    # Load loader and parse
    parse_cw_temp_output = find_loader()
    print(f"\nLoading {output_path}...")
    print("(this may take 30-60 seconds for large files)")
    v = parse_cw_temp_output(output_path)
    nt, nd, nw = v.T_field_F.shape
    print(f"  T_field shape: {v.T_field_F.shape}  (n_time, n_depth, n_width)")
    print(f"  time range:    {v.time_hrs[0]:.2f} to {v.time_hrs[-1]:.2f} hrs")
    print(f"  widths_m:      {v.widths_m[0]:.2f} to {v.widths_m[-1]:.2f} m  ({nw} nodes)")
    print(f"  depths_m:      {v.depths_m[0]:.2f} to {v.depths_m[-1]:.2f} m  ({nd} nodes)")
    print(f"  T range:       {v.T_field_F.min():.1f} to {v.T_field_F.max():.1f} F")
    print(f"  T_ambient:     {v.T_ambient_F.min():.2f} to {v.T_ambient_F.max():.2f} F")

    # Sanity check: ambient should be near-constant for a proper adiabatic run
    amb_range = v.T_ambient_F.max() - v.T_ambient_F.min()
    if amb_range > 1.0:
        print(f"  WARNING: T_ambient varies by {amb_range:.2f} F — was the CW "
              "scenario set up for adiabatic isolation? Expected constant ambient.")

    # Find the centerline
    cl_w, cl_d = find_centerline_axis(v.widths_m, v.depths_m, v.T_field_F)
    print(f"\nDetected centerline:")
    print(f"  width index  = {cl_w} (out of {nw-1})  -> width = {v.widths_m[cl_w]:.2f} m")
    print(f"  depth index  = {cl_d} (out of {nd-1})  -> depth = {v.depths_m[cl_d]:.2f} m")
    print(f"  Final T at this point: {v.T_field_F[-1, cl_d, cl_w]:.2f} F")

    # Spatial uniformity check at the centerline column
    cl_column = v.T_field_F[-1, :, cl_w]  # full depth column at the centerline
    interior_band = cl_column[nd//4 : 3*nd//4]
    interior_uniformity = interior_band.max() - interior_band.min()
    print(f"  Centerline column interior uniformity: {interior_uniformity:.2f} F"
          f"  ({'GOOD' if interior_uniformity < 0.5 else 'CHECK BOUNDARIES'})")

    # Extract trajectories
    T_center = v.T_field_F[:, cl_d, cl_w]            # adiabatic core
    T_max_xs = v.T_max_xs_F                          # max anywhere in xs
    T_min_xs = v.T_min_xs_F                          # min anywhere in xs

    # ---- WRITE PRIMARY REFERENCE CSV ----
    ref_path = os.path.join(save_dir, f"cw_adiabatic_reference_{label}.csv")
    np.savetxt(
        ref_path,
        np.column_stack([v.time_hrs, T_center, T_max_xs, v.T_ambient_F]),
        header="time_hrs,T_center_F_adiabatic,T_max_xs_F,T_ambient_F",
        delimiter=",", comments="", fmt="%.4f",
    )
    print(f"\nSaved {ref_path}")

    # ---- WRITE DIAGNOSTIC CSV (extra spatial probes) ----
    # mid-depth temps at several width positions, for spatial uniformity check
    diag_path = os.path.join(save_dir, f"adiabatic_diagnostic_{label}.csv")
    np.savetxt(
        diag_path,
        np.column_stack([
            v.time_hrs,
            v.T_field_F[:, 0, cl_w],          # centerline TOP surface
            v.T_field_F[:, cl_d, cl_w],       # centerline mid-depth (= ref)
            v.T_field_F[:, -1, cl_w],         # centerline bottom
            v.T_field_F[:, cl_d, nw//2],      # mid-depth, mid-width
            v.T_field_F[:, cl_d, 0 if cl_w != 0 else nw - 1],  # mid-depth at form face
        ]),
        header="time_hrs,c_top,c_mid,c_bot,mid_quarter,mid_corner",
        delimiter=",", comments="", fmt="%.4f",
    )
    print(f"Saved {diag_path}")

    # ---- SUMMARY ----
    print(f"\n=== Summary for {label} ===")
    peak_idx = int(np.argmax(T_center))
    print(f"  Peak adiabatic T:     {T_center.max():.2f} F at t = {v.time_hrs[peak_idx]:.2f} hr")
    print(f"  Initial T:            {T_center[0]:.2f} F")
    print(f"  Adiabatic rise:       {T_center.max() - T_center[0]:.2f} F")
    print(f"  Ambient T (constant): {v.T_ambient_F.mean():.2f} F")

    # Time-to-temperature milestones
    t0 = T_center[0]
    for target_above_initial in [10.0, 25.0, 50.0]:
        target = t0 + target_above_initial
        if T_center.max() > target:
            t_reach = float(np.interp(target, T_center, v.time_hrs))
            print(f"  Time to T = {target:.0f} F (T0+{target_above_initial:.0f}): {t_reach:.2f} hr")

    print(f"\nUse `cw_adiabatic_reference_{label}.csv` for engine kinetics validation.")
    print(f"Column 'T_center_F_adiabatic' is the kinetics ground truth.")


if __name__ == "__main__":
    main()
