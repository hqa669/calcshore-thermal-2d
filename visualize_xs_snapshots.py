"""
Visualize CW 2D temperature distribution at fixed time intervals.

Generates a grid of subplots showing the half-mat cross-section temperature
field at every N hours. Each subplot is a heatmap with depth on the y-axis
and width on the x-axis (width=0 is the centerline / symmetry plane).

Run from a directory where cw_scenario_loader.py is importable:
    python visualize_xs_snapshots.py
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# ---- CONFIG ----
SCENARIO_DIR = os.path.expanduser("~/Downloads/HydrationCenter_mix01")
OUTPUT_FILE = os.path.join(SCENARIO_DIR, "output.txt")

# Time interval for snapshots (hours)
SNAPSHOT_INTERVAL_HR = 5.0

# Output figure
FIG_PATH = os.path.join(SCENARIO_DIR, "xs_temperature_snapshots.png")

# Optional: zoom in on the corner region (where boundary effects live)
# Set to None to show the full half-mat. For 200x200 ft, the full mat is mostly
# uniform — set to e.g. 0.1 to show only the corner-most 10% of the width.
ZOOM_FRACTION_FROM_FORM_FACE = None  # e.g., 0.05 zooms to the corner-most 5%

# ---- LOAD DATA ----
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cw_scenario_loader import parse_cw_temp_output

print(f"Loading {OUTPUT_FILE}...")
v = parse_cw_temp_output(OUTPUT_FILE)
nt, nd, nw = v.T_field_F.shape
print(f"  Shape: {v.T_field_F.shape}  (n_time, n_depth, n_width)")
print(f"  Time:  {v.time_hrs[0]:.2f} to {v.time_hrs[-1]:.2f} hrs")
print(f"  Depth: {v.depths_m[0]:.2f} to {v.depths_m[-1]:.2f} m  ({nd} nodes)")
print(f"  Width: {v.widths_m[0]:.2f} to {v.widths_m[-1]:.2f} m  ({nw} nodes)")
print(f"  T:     {v.T_field_F.min():.1f} to {v.T_field_F.max():.1f} F")

# ---- PICK SNAPSHOT TIMES ----
t_targets = np.arange(0.0, v.time_hrs[-1] + 0.01, SNAPSHOT_INTERVAL_HR)
# Find the time index closest to each target
snapshot_indices = [int(np.abs(v.time_hrs - tt).argmin()) for tt in t_targets]
n_snap = len(snapshot_indices)
print(f"\nGenerating {n_snap} snapshots at {SNAPSHOT_INTERVAL_HR} hr intervals")

# ---- ZOOM (OPTIONAL) ----
if ZOOM_FRACTION_FROM_FORM_FACE is not None:
    n_zoom_cols = max(2, int(nw * ZOOM_FRACTION_FROM_FORM_FACE))
    w_slice = slice(nw - n_zoom_cols, nw)
    print(f"  Zoom: showing last {n_zoom_cols} of {nw} width nodes (corner)")
else:
    w_slice = slice(0, nw)
field_slice = v.T_field_F[:, :, w_slice]
widths_plot = v.widths_m[w_slice]

# ---- PLOT ----
# Use a common color scale across all snapshots for fair comparison
T_min_global = field_slice.min()
T_max_global = field_slice.max()
print(f"  Color scale: {T_min_global:.1f} to {T_max_global:.1f} F\n")

# Grid of subplots: aim for roughly square layout
n_cols = int(np.ceil(np.sqrt(n_snap)))
n_rows = int(np.ceil(n_snap / n_cols))

fig, axes = plt.subplots(n_rows, n_cols, figsize=(3.5 * n_cols, 2.8 * n_rows),
                         squeeze=False)

# CW convention: width axis 0 = centerline, last index = form face
# Depth axis 0 = top, last = bottom
# Use widths in feet for readability
widths_ft = widths_plot * 3.28084
depths_ft = v.depths_m * 3.28084

extent = [widths_ft[0], widths_ft[-1], depths_ft[-1], depths_ft[0]]
# extent: [left, right, bottom, top] in data coords
# We want depth=0 (top) at the TOP of the figure → y-axis inverted via extent

for i_ax, (t_idx, ax) in enumerate(zip(snapshot_indices, axes.flat)):
    field = field_slice[t_idx]  # shape (nd, nw_slice)
    im = ax.imshow(field, aspect='auto', extent=extent,
                   vmin=T_min_global, vmax=T_max_global,
                   cmap='inferno', origin='upper')
    ax.set_title(f"t = {v.time_hrs[t_idx]:.1f} hr   "
                 f"(min={field.min():.1f}, max={field.max():.1f} F)",
                 fontsize=9)
    if i_ax % n_cols == 0:
        ax.set_ylabel("Depth (ft)\n(top → bottom)", fontsize=8)
    if i_ax >= (n_rows - 1) * n_cols:
        ax.set_xlabel("Width (ft)\n(0 = centerline, → form face)", fontsize=8)
    ax.tick_params(labelsize=7)

# Hide any unused subplots
for j in range(n_snap, n_rows * n_cols):
    axes.flat[j].set_visible(False)

# Shared colorbar
fig.subplots_adjust(right=0.90, hspace=0.45, wspace=0.3, top=0.94, bottom=0.08)
cbar_ax = fig.add_axes([0.92, 0.10, 0.015, 0.80])
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.set_label('Temperature (°F)', fontsize=9)
cbar.ax.tick_params(labelsize=8)

zoom_label = f" (zoomed to corner, last {len(widths_plot)}/{nw} cols)" if \
    ZOOM_FRACTION_FROM_FORM_FACE else ""
fig.suptitle(f"CW 2D Temperature Distribution — every {SNAPSHOT_INTERVAL_HR} hr"
             + zoom_label,
             fontsize=11, fontweight='bold')

print(f"Saving figure to {FIG_PATH}...")
plt.savefig(FIG_PATH, dpi=120, bbox_inches='tight')
print("Done.")

# ---- ALSO: 1D CENTERLINE COLUMN PROFILE AT EACH SNAPSHOT ----
# This shows how temperature varies with depth at the centerline (width=0).
# In a perfectly adiabatic interior, this should be flat.
fig2, ax2 = plt.subplots(1, 1, figsize=(8, 6))
cmap = plt.cm.viridis
for k, t_idx in enumerate(snapshot_indices):
    color = cmap(k / max(1, n_snap - 1))
    centerline_col = v.T_field_F[t_idx, :, 0]  # all depths, width=0
    ax2.plot(centerline_col, depths_ft, color=color,
             label=f"t={v.time_hrs[t_idx]:.0f} hr", linewidth=1.5)
ax2.invert_yaxis()
ax2.set_xlabel("Temperature (°F)")
ax2.set_ylabel("Depth (ft)  (top → bottom)")
ax2.set_title("Centerline column temperature profile vs depth\n"
              "(width=0, every 5 hr)")
ax2.grid(True, alpha=0.3)
ax2.legend(loc='center right', fontsize=8, ncol=2)
fig2.tight_layout()
fig2_path = os.path.join(SCENARIO_DIR, "centerline_column_profiles.png")
plt.savefig(fig2_path, dpi=120, bbox_inches='tight')
print(f"Saved {fig2_path}")

print("\nDone. Open the PNG files to inspect.")
