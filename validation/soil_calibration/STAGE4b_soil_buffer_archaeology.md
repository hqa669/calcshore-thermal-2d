# STAGE4b — Soil Buffer Archaeology

## Question

Who added the 3 m soil buffer below the concrete, in what commit, with what
stated intent?  Is the 3 m / 15-cell specification deliberate or an
undocumented incidental default?

## Findings

### Soil row allocation (grid builder)

`git blame -L 976,981 thermal_engine_2d.py` attributes lines 976–981 to
commit **`fd19c9a`** (2026-04-22):

```
M0: 2D grid builder + domain composition
```

This is the very first commit in the project — the soil buffer was present
from the initial grid design.  The commit message says:

> - Grid2D dataclass with material_id, masks, and index metadata
> - build_grid_half_mat() with primitive args, uniform dx, non-uniform dy
> - 11 pytest unit tests covering shape, extents, topology, CW width alignment
> - Interface node ownership convention: interface rows belong to concrete

No explicit rationale is given for the 3 m / 15-cell defaults (`n_soil_y=15`,
`soil_depth_below_m=3.0`).  The parameters appear as keyword-argument defaults
in `build_grid_half_mat` without any docstring or comment explaining the
physical justification.

The M0 code comment at lines 976–979 reads:

```python
# Soil below concrete: n_soil_y nodes, uniform dy_soil, starting one dy_soil below concrete bottom.
dy_soil = soil_depth_below_m / n_soil_y
y_concrete_bottom = blanket_thickness_m + depth_m
y_soil = y_concrete_bottom + dy_soil * np.arange(1, n_soil_y + 1)
```

This describes the mechanics of the allocation (uniform spacing, starting one
`dy_soil` below the concrete face) but does not state the physical reason for
choosing 3 m or 15 cells.

### Deep BC (bottom-row Dirichlet)

`git blame -L 2183,2184 thermal_engine_2d.py` attributes those lines to commit
**`f7b1826`** (2026-04-22):

```
M4: full 2D BCs, Sprint 0 complete
```

M4 was the same day as M0 (the initial "Sprint 0 complete" milestone).  The
commit message includes:

> - Far-left soil Dirichlet T_gw

The code comment added in M4 reads:

> Design doc D5: "assume T = T_gw at x=0 in soil (soil extends to deep ground
> temperature far from concrete)." Dirichlet rather than adiabatic prevents the
> cold soil-extension from acting as a lateral heat sink for the form-edge
> concrete bottom — matching CW's domain which has no soil to the left of the
> form face.

This confirms that the original design intent was `T = T_gw` at the far
boundary of the soil domain.  The "far" boundary was assumed to be 3 m below
the concrete face (the M0 default), but the M4 comment does not explain why
3 m was chosen or how it relates to CW.

### Summary: deliberate or undocumented?

**Deliberate, but not justified against CW.**  The soil buffer and its 3 m
depth were present from M0 as a design choice (not an accidental omission), and
M4 added an explicit Dirichlet BC at that depth with a reference to "Design doc
D5."  However, the buffer depth was never calibrated against CW — CW does not
model a soil domain at all and applies `T_soil` directly at the concrete face
(as Stage 4a confirmed).  The M4 comment explicitly notes that the left-side
Dirichlet "matching CW's domain which has no soil to the left of the form face"
— ironically, the analogous bottom-face match was never made.

The 3 m / 15-cell specification appears to be a reasonable engineering estimate
for far-field ground temperature (thermal diffusion length for Δt=168 hr ≈
√(α·t) ≈ √(2.86×10⁻³ × 168) ≈ 0.69 m, so 3 m ≈ 4.3 penetration depths) but
was not documented as such.

## Implication for `model_soil` docstring

The `model_soil=True` path represents the **original engine behavior (M0+M4)**:
a transient soil mesh with a 3 m buffer below (and, when `is_submerged=True`,
alongside) the concrete.  The soil temperature at depth is held at `T_soil`
via a Dirichlet BC.

`model_soil=False` (new default) is the **CW-matching path**: no soil domain
at all; `T_soil` applied as a Dirichlet BC directly at the concrete face.
