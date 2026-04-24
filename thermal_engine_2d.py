"""
CalcShore Thermal Engine v3 — 2D Finite-Difference Solver
==========================================================
Milestone M0: Grid builder and domain composition.
Milestone M1: Explicit 2D conduction solver, constant properties,
              Dirichlet BC, analytical Fourier-series validation helper.
Milestone M2: Hydration heat + variable properties (Schindler-Folliard).
Milestone M3: Top boundary condition — convection + blanket R + Menzel evaporation.
Milestone M4: Geometry rework (air/inactive region), top BC fix, side BC
              (form face), centerline symmetry BC, far-left soil Dirichlet,
              boundary_mode='full_2d'.  Sprint-0 complete.

M4 architecture decisions:
  - Blanket is modelled as pure thermal resistance (not a physical node) in
    full_2d mode.  The blanket row (j=0) exists in the grid for index
    alignment but carries no thermal mass; the top BC is applied directly at
    j=iy_concrete_start using h_top_series(wind, blanket_R).
  - Side BC at x=0 (form face): h_side = 1/(1/h_conv + R_FORM_EFFECTIVE_SI)
    where R_FORM_EFFECTIVE_SI = 0.0862 m²·K/W (empirically calibrated against
    CW MIX-01; see Sprint 2 for ACI Eq 27 orientation model revision).
  - Corner cell (iy_concrete_start, ix_concrete_start) receives a quarter-cell
    energy balance that accounts simultaneously for the top and side BCs.

Known limitations (deferred to later sprints):
  - Corner RMS ~4.0°F on MIX-01 (tolerance 3.0°F).  Root cause: no solar
    radiation on the form face.  CW's afternoon corner temperature exceeds
    ambient by several °F due to solar gain that our model does not include.
    Fix: Sprint 1 (solar + longwave top BC).
  - Peak Max T runs ~1.0°F below CW on MIX-01 (at the ±1.0°F boundary).
    Same root cause: no solar gain during warm hours depresses the peak
    slightly.  Fix: Sprint 1.
  - Peak Max T timing ~6 hr early vs CW on MIX-01.  Same root cause: without
    solar augmentation the hydration peak occurs before CW's solar-driven
    profile does.  Fix: Sprint 1.

Coordinate convention (half-mat, vertex-centered):
  x = 0      : concrete/form edge (soil-extension at x < 0)
  x = W/2    : centerline (symmetry BC)
  y = 0      : top of blanket row (blanket is pure-R in full_2d mode)
  y > 0      : downward into the structure

Material IDs: 0 = blanket, 1 = concrete, 2 = soil, 3 = air/inactive (not solved)

============================================================
BUGS FOUND IN reference/thermal_engine_v2.py DURING M3 PORT
============================================================
Three bugs were discovered while porting v2's top BC to this 2D engine.
Flag these upstream to v2 before its next release.

(a) ambient_temp_F SIGN FLIP
    v2's formula: avg - amp * cos(2π*(h-15)/24)
    This gives MINIMUM at h=15 (3 PM) and MAXIMUM at h=3 (3 AM) —
    the opposite of physical reality and CW's own T_ambient_F output.
    Fix applied here: avg + amp * cos(2π*(h-15)/24).
    Impact in v2: ambient temperature at peak-concrete time (~7 AM) is
    ~12°F too warm, reducing the computed surface heat loss and inflating
    the simulated peak temperature.

(b) Menzel evaporation mmHg / kPa UNIT MISMATCH
    v2 calls saturated_vapor_pressure_mmHg() (returns mmHg, ~17 mmHg at
    20°C) but passes the result into Ew = 0.315 * (e0 - RH*ea) * f(V),
    where the coefficient 0.315 was calibrated for vapor pressures in kPa
    (~2.3 kPa at 20°C). Using mmHg inflates the evaporation rate by
    ~7.5× (mmHg/kPa ratio), producing ~230 W/m² instead of ~31 W/m².
    Fix applied here: saturated_vapor_pressure_kPa() wrapper; Menzel
    uses kPa throughout.
    Impact in v2: enormous first-hour evaporative flux — partially masked
    in v2 by the blanket's explicit thermal mass absorbing the spike, but
    physically incorrect.

(c) Ghost-node BC NUMERICAL INSTABILITY with nonlinear Menzel
    The ghost-node formula T[0] = T[1] - q_top·dy/k is conditionally
    unstable when Menzel activates: the effective surface conductance
    (dq_evap/dT ≈ 15 W/(m²·K) even with the kPa fix) exceeds the
    stability limit k_face/dy_top ≈ 2 W/(m²·K) for the thin blanket
    node. This causes the surface temperature to oscillate and blow up
    within a few time steps.
    Fix applied here: half-cell energy balance for row 0, which carries
    the blanket's thermal inertia (ρCp·dy/2) and is stable at the CFL
    dt. This matches v2's own approach (v2 uses explicit blanket nodes
    with thermal mass, not a ghost-node).
    Impact in v2: v2 avoids this by using explicit blanket nodes — the
    instability would surface if v2's top BC were ever rewritten as a
    ghost-node.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np

FT_TO_M = 0.3048

# Physical constants
# Copied from reference/thermal_engine_v2.py:33-34, verbatim.
R_GAS = 8.314       # J/(mol·K)
T_REF_K = 296.15    # 23°C — Arrhenius reference temperature

# Imperial-to-SI unit conversions (single source of truth)
LB_FT3_TO_KG_M3 = 16.0185
BTU_LB_F_TO_J_KG_K = 4186.8
BTU_HR_FT_F_TO_W_M_K = 1.7307
LB_YD3_TO_KG_M3 = 0.593276
R_IMP_TO_R_SI = 0.1761   # (hr·ft²·°F/BTU) → (m²·K/W)

# Steel form thermal resistance (negligible for bare steel).
FORM_R_SI = 0.0   # m²·K/W

# Effective side-face resistance calibrated against CW MIX-01 (steel form +
# wet-concrete contact film, Sprint 0).  Sprint 2 will revisit with the ACI
# Eq 27 orientation-dependent model; removing it caused a 10°F gradient
# regression (full LW + conv overcooled the side face at night).
R_FORM_EFFECTIVE_SI = 0.0862   # m²·K/W  (= 0.490 hr·ft²·°F/BTU × 0.1761)

# Vertical-surface solar projection factor: fraction of horizontal irradiance
# incident on a vertical form face. F_vert = 0.5 is the correct flat-projection
# value for a randomly oriented vertical surface at mid-latitude (geometry only).
#
# WHY CWConstruction.vertical_solar_factor DEFAULTS TO 0.0 (DARK INFRASTRUCTURE):
#
# An empirical F_vert sweep on MIX-01 (168-hr solve, R_FORM_EFFECTIVE_SI in place)
# showed Corner RMS monotonically worsening with increasing F_vert — the opposite
# of what the retrospective "missing solar" hypothesis predicted:
#
#   F_vert | Corner RMS | Peak Grad Δ  | Peak Max T Δ
#   -------|------------|-------------|-------------
#     0.0  |  4.10°F   |  -0.95°F ✓  |  -0.33°F ✓
#     0.3  |  9.23°F   |  -5.67°F ✗  |  -0.33°F ✓
#     0.5  | 14.58°F   |  -4.78°F ✗  |  +0.02°F ✓
#     0.7  | 20.13°F   |  +4.35°F ✗  |  +7.76°F ✗
#     1.0  | 28.58°F   |  +17.94°F ✗ | +22.28°F ✗
#
# Diagnosis: CW's corner diurnal mismatch is NOT a missing-solar problem; it is a
# missing-LW-at-night problem on the vertical form face. Engine corner is ~3°F too
# warm vs CW at night (not at solar noon — they match at noon), meaning CW applies
# nighttime LW cooling at the form face that the Sprint 1 model lacks. Increasing
# F_vert injects daytime solar heat that pushes corner further from CW by day while
# doing nothing to close the nighttime gap. F_vert=0.0 is the best-performing value
# across all metrics and the only configuration where all other S0 gates pass.
#
# SPRINT 2 SOLUTION PATH: T_outer solve on the vertical form face (analogous to
# Option-B Newton solve on the horizontal blanket outer surface from PR 3). Side-face
# LW q = ε·σ·(T_outer_form^4 − T_sky^4) applied at the OUTER form surface, conducted
# through R_FORM_EFFECTIVE_SI to the concrete, will provide the missing nighttime
# radiative cooling without bloating the daytime temperature. Once that is in place,
# F_vert can be re-swept against the corrected thermal baseline to calibrate the
# solar contribution correctly. ACI 207.2R Eq 27 provides the orientation model.
VERTICAL_SOLAR_FACTOR_DEFAULT = 0.5   # dimensionless — correct geometry; see above

# Solar absorptivity and emissivity defaults for the top concrete surface.
# absorptivity: ASHRAE default for light-gray concrete/steel form.
SOLAR_ABSORPTIVITY_DEFAULT = 0.65
EMISSIVITY_DEFAULT = 0.88
STEFAN_BOLTZMANN = 5.670374419e-8   # W/(m²·K⁴)

# View factor from a vertical form face to the sky hemisphere.
# For a vertical flat surface in contact with horizontal ground, half
# the hemispheric solid angle sees sky and half sees ground.
# Used in PR 6 side-face LW balance:
#   q_lw_form = ε·σ·[F_SKY_VERT·(T_form⁴ − T_sky⁴)
#                  + (1−F_SKY_VERT)·(T_form⁴ − T_ground⁴)]
# Ref: Modest "Radiative Heat Transfer" 3rd ed., view-factor algebra
# for perpendicular planes in infinite geometry.
F_SKY_VERT: float = 0.5

# Emissivity of bare soil / aggregate ground surface.
# ASHRAE Fundamentals 2021 Ch. 14 Table 1: "Dry Ground" ε ≈ 0.92-0.96,
# typical near-blackbody value for natural terrain.
# Hardcoded (not parameterized) in Sprint 2; parameterizes in Sprint 3
# when Barber soil model lands and surface type becomes an input.
EMIS_GROUND: float = 0.95

# Default material properties for non-concrete zones
SOIL_PROPERTIES_2D: dict = {'k': 1.50, 'rho': 2100.0, 'Cp': 900.0}
BLANKET_PROPERTIES_2D: dict = {'rho': 100.0, 'Cp': 1500.0}
# k_blanket is derived from R-value and thickness at solver call-site.


# ============================================================
# Hydration physics (Schindler-Folliard model)
# ============================================================

def arrhenius_vec(T_K: np.ndarray, Ea: float) -> np.ndarray:
    """Arrhenius acceleration factor (dimensionless, vectorized).

    Parameters
    ----------
    T_K : array-like
        Temperature in Kelvin.
    Ea : float
        Activation energy J/mol.

    Returns
    -------
    np.ndarray  dimensionless acceleration factor relative to T_REF_K (23°C).

    # Copied from reference/thermal_engine_v2.py:209-211, verbatim.
    # DO NOT edit in place — changes must be made to both files and documented.
    """
    return np.exp(-Ea / R_GAS * (1.0 / T_K - 1.0 / T_REF_K))


def hydration_rate_vec(te: np.ndarray, tau: float, beta: float, au: float) -> np.ndarray:
    """dα/d(te) — hydration rate per equivalent hour (vectorized).

    Parameters
    ----------
    te : array-like
        Equivalent age in hours.
    tau : float
        Hydration time parameter (hours).
    beta : float
        Hydration shape parameter (dimensionless).
    au : float
        Ultimate degree of hydration (dimensionless).

    Returns
    -------
    np.ndarray  dα/dte in units of [1/hour].

    # Copied from reference/thermal_engine_v2.py:213-217, verbatim.
    # DO NOT edit in place — changes must be made to both files and documented.
    """
    te_safe = np.maximum(te, 0.01)
    r = tau / te_safe
    return au * beta * r**beta / te_safe * np.exp(-r**beta)


def hydration_alpha_vec(te: np.ndarray, tau: float, beta: float, au: float) -> np.ndarray:
    """α(te) — degree of hydration as a function of equivalent age (vectorized).

    Parameters
    ----------
    te : array-like
        Equivalent age in hours.
    tau, beta, au : float
        Schindler-Folliard model parameters.

    Returns
    -------
    np.ndarray  α dimensionless, range [0, au].

    # Copied from reference/thermal_engine_v2.py:219-222, verbatim.
    # DO NOT edit in place — changes must be made to both files and documented.
    """
    te_safe = np.maximum(te, 0.01)
    return au * np.exp(-(tau / te_safe)**beta)


def thermal_conductivity_variable(k_uc: float, alpha_node: np.ndarray) -> np.ndarray:
    """CW Eq 23: thermal conductivity as a function of hydration degree.

    Parameters
    ----------
    k_uc : float
        Uncured concrete thermal conductivity W/(m·K).
    alpha_node : array-like
        Degree of hydration per node (dimensionless).

    Returns
    -------
    np.ndarray  k in W/(m·K). Ranges from 1.33·k_uc at α=0 to k_uc at α=α_u≈1.

    # Copied from reference/thermal_engine_v2.py:281-283, verbatim.
    # DO NOT edit in place — changes must be made to both files and documented.
    """
    return k_uc * (1.33 - 0.33 * alpha_node)


def specific_heat_variable(
    Wc: float,
    Wa: float,
    Ww: float,
    Ca: float,
    alpha_node: np.ndarray,
    T_C_node: np.ndarray,
    rho: float,
) -> np.ndarray:
    """CW Eq 24-25: specific heat as a function of mix, hydration, and temperature.

    Van Breugel model. All weights are per unit volume of concrete (kg/m³).

    Parameters
    ----------
    Wc : float
        Cement weight kg/m³.
    Wa : float
        Aggregate weight kg/m³.
    Ww : float
        Water weight kg/m³.
    Ca : float
        Aggregate specific heat J/(kg·K).
    alpha_node : array-like
        Degree of hydration (dimensionless).
    T_C_node : array-like
        Concrete temperature in °C per node.
    rho : float
        Concrete density kg/m³.

    Returns
    -------
    np.ndarray  Cp in J/(kg·K).

    # Copied from reference/thermal_engine_v2.py:285-298, verbatim.
    # DO NOT edit in place — changes must be made to both files and documented.
    """
    Cw = 4186.0                             # water specific heat J/(kg·K)
    c_cef = 8.4 * T_C_node + 339.0         # fictitious specific heat of hydrating cement
    Cp_node = (1.0 / rho) * (
        Wc * alpha_node * c_cef
        + Wc * (1.0 - alpha_node) * Ca
        + Wa * Ca
        + Ww * Cw
    )
    return Cp_node


# ============================================================
# M3 physics helpers — top boundary condition
# ============================================================

# Copied from reference/thermal_engine_v2.py:241-251, verbatim.
def saturated_vapor_pressure_mmHg(T_C: float) -> float:
    """ASHRAE saturation vapor pressure.

    Parameters
    ----------
    T_C : float  Temperature in °C.

    Returns
    -------
    float  Saturation vapor pressure in mmHg.
    """
    T_K = T_C + 273.15
    if T_C >= 0:
        C8, C9, C10, C11, C12, C13 = -5.8002206e3, -5.516256, -4.8640239e-2, 4.1764768e-5, -1.4452093e-8, 6.5459673
        ln_Pws = C8/T_K + C9 + C10*T_K + C11*T_K**2 + C12*T_K**3 + C13*np.log(T_K)
    else:
        C1, C2, C3, C4, C5, C6, C7 = -5.6745359e3, -5.1523058e-1, -9.677843e-3, 6.2215701e-7, 2.0747825e-9, -9.484024e-13, 4.1635019
        ln_Pws = C1/T_K + C2 + C3*T_K + C4*T_K**2 + C5*T_K**3 + C6*T_K**4 + C7*np.log(T_K)
    Pws_kPa = np.exp(ln_Pws)
    return Pws_kPa * 7.50062  # kPa to mmHg


def saturated_vapor_pressure_kPa(T_C: float) -> float:
    """ASHRAE saturation vapor pressure in kPa. Wrapper around mmHg version.

    Parameters
    ----------
    T_C : float  Temperature in °C.

    Returns
    -------
    float  Saturation vapor pressure in kPa.
    """
    return saturated_vapor_pressure_mmHg(T_C) / 7.50062


# Copied from reference/thermal_engine_v2.py:253-275, adapted.
# Change: uses saturated_vapor_pressure_kPa instead of mmHg, because the
# Menzel coefficient 0.315 was calibrated for vapor pressures in kPa.
# Using mmHg (7.5× larger) overestimates evaporation 7.5× and causes
# numerical instability in the ghost-node BC. Confirmed by comparing against
# CW ambient T_ambient_F series for MIX-01.
# Cure-time suppression is applied by the caller, not here.
def menzel_evaporation(t_hrs: float, T_surface_C: float, T_air_C: float,
                       RH: float, wind_m_s: float) -> float:
    """Menzel evaporative cooling rate. CW Eq 38-39.

    Parameters
    ----------
    t_hrs : float     Hours since concrete placement.
    T_surface_C : float  Concrete surface temperature °C.
    T_air_C : float   Ambient air temperature °C.
    RH : float        Relative humidity as fraction 0–1.
    wind_m_s : float  Wind speed m/s (already derated by caller if needed).

    Returns
    -------
    float  Evaporative heat flux W/m² (positive = heat loss from surface).
    """
    e0 = saturated_vapor_pressure_kPa(T_surface_C)
    ea = saturated_vapor_pressure_kPa(T_air_C)

    Ew = 0.315 * (e0 - RH * ea) * (0.253 + 0.060 * wind_m_s)
    Ew = max(Ew, 0.0)

    # Bleed water depletion decay
    a_evap = 3.75  # hours
    if t_hrs < 0.01:
        decay = 1.0
    else:
        decay = np.exp(-(t_hrs / a_evap)**1.5)

    Ec = Ew * decay

    hfg = 2_260_000  # J/kg latent heat of vaporization
    q_evap = Ec * hfg / 3600.0  # W/m²
    return q_evap


def h_forced_convection(wind_m_s_max: float) -> float:
    """Pure forced-convection coefficient for the top surface.

    Does NOT include a longwave-equivalent term; longwave is handled
    explicitly in PR 3. For legacy regression/ablation modes that expect
    the combined coefficient, use _h_convective_legacy instead.

    Applies 0.4× wind derating internally to convert avg-max wind to
    effective surface wind, matching reference/thermal_engine_v2.py:591.

    Parameters
    ----------
    wind_m_s_max : float  Average-max wind speed m/s (env.cw_ave_max_wind_m_s).

    Returns
    -------
    float  Forced-convection heat transfer coefficient W/(m²·K).
    """
    return 5.6 + 3.5 * (0.4 * wind_m_s_max)


def _h_convective_legacy(wind_m_s_max: float) -> float:
    """LEGACY: convection with baked-in longwave equivalent (+5.5 W/(m²·K)).

    DO NOT USE for new validation. Preserved for regression/ablation
    testing of top_bc_only and v2_equivalent boundary modes only.
    Production path is h_forced_convection + explicit LW in full_2d mode.

    Formula: h = 5.6 + 3.5*(0.4*wind) + 5.5 (the trailing 5.5 is an
    empirical longwave-equivalent term from v2).

    Parameters
    ----------
    wind_m_s_max : float  Average-max wind speed m/s.

    Returns
    -------
    float  Combined convection + LW-equivalent coefficient W/(m²·K).
    """
    return 5.6 + 3.5 * (0.4 * wind_m_s_max) + 5.5


def h_top_series(wind_m_s_max: float, blanket_R_imp: float) -> float:
    """Series-resistance combination of pure forced-convection + cure blanket.

    NOTE (M4+): This function is no longer called by the solver. In M4 the
    blanket is modelled as pure-R resistance in full_2d mode; the solver
    precomputes _h_top_combined directly. This function is retained as a
    unit-tested mathematical reference (pure forced-convection + blanket R,
    no LW-equivalent term).

    NOTE (PR 2+): Uses h_forced_convection (no +5.5 LW term). The legacy
    formula (with +5.5) is h_top_series with _h_convective_legacy, which the
    solver precomputes as _h_top_legacy for top_bc_only / v2_equivalent modes.

    Parameters
    ----------
    wind_m_s_max : float  Average-max wind speed m/s.
    blanket_R_imp : float  Blanket R-value in hr·ft²·°F/BTU.

    Returns
    -------
    float  Effective top-surface heat transfer coefficient W/(m²·K).
    """
    h_conv = h_forced_convection(wind_m_s_max)
    R_blanket_SI = blanket_R_imp * R_IMP_TO_R_SI  # m²·K/W
    return 1.0 / (1.0 / h_conv + R_blanket_SI)


def h_side_convective(wind_m_s_max: float, blanket_R_imp: float = 0.0) -> float:
    """Effective heat transfer coefficient at the concrete form face.

    Physical path: ambient → convection → form face → concrete.
    R_FORM_EFFECTIVE_SI captures the empirically calibrated combined resistance
    of the steel form and wet-concrete contact film (≈ 0.49 hr·ft²·°F/BTU).
    The cure blanket is a top-surface treatment only; blanket_R_imp is accepted
    for call-site compatibility but does not affect the side BC result.

    PR 4 keeps R_FORM_EFFECTIVE_SI for convection but adds solar (F_vert
    projection with daytime gate). Removing R_FORM entirely caused a 10°F
    gradient regression due to nighttime overcooling from unattenuated LW + h_forced.
    Sprint 2 will replace R_FORM with the ACI 207.2R Eq 27 orientation model.

    Parameters
    ----------
    wind_m_s_max : float  Average-max wind speed m/s.
    blanket_R_imp : float  Accepted but ignored; blanket is on top only.

    Returns
    -------
    float  Effective side-surface heat transfer coefficient W/(m²·K).
    """
    return 1.0 / (1.0 / h_forced_convection(wind_m_s_max) + R_FORM_EFFECTIVE_SI)


# Adapted from reference/thermal_engine_v2.py:197-203.
# Change: placement_hour passed as explicit argument (lives on CWConstruction
# in this project, not CWEnvironment).
# Change (PR 2): when env.T_air_F is populated, use hourly linear interpolation
# instead of the sinusoidal formula. The sinusoid remains as a fallback for
# test fixtures that construct minimal CWEnvironment objects without weather data.
# Original sign correction note: `+amp * cos(...)` so maximum is at h=15 (3 PM)
# and minimum at h=3 (3 AM), matching physical reality and CW ambient output.
def ambient_temp_F(t_hrs: float, env, placement_hour: int) -> float:
    """Ambient temperature at time t_hrs since placement, in °F.

    When env.T_air_F is populated (hourly weather data loaded), returns
    linearly interpolated value from the hourly array. Otherwise falls
    back to sinusoidal daily min/max formula.

    The hourly path is the production path for CW-loaded scenarios; the
    sinusoidal fallback exists only for synthetic test fixtures.

    Parameters
    ----------
    t_hrs : float        Hours since concrete placement.
    env : duck-typed     For hourly path: must have T_air_F, hours (aligned).
                         For sinusoidal fallback: must have daily_max_F,
                         daily_min_F (per-day lists).
    placement_hour : int Hour of day (0–23) when placement started.
                         Unused on the hourly path (t_hrs=0 maps to placement).

    Returns
    -------
    float  Ambient temperature in °F.
    """
    if getattr(env, 'T_air_F', None) is not None and np.asarray(env.T_air_F).size > 0:
        return float(np.interp(t_hrs, env.hours, env.T_air_F))
    # Sinusoidal fallback for test fixtures without weather data:
    h = (placement_hour + t_hrs) % 24.0
    d = min(int((placement_hour + t_hrs) / 24.0), len(env.daily_max_F) - 1)
    avg = (env.daily_max_F[d] + env.daily_min_F[d]) / 2
    amp = (env.daily_max_F[d] - env.daily_min_F[d]) / 2
    return avg + amp * np.cos(2 * np.pi * (h - 15.0) / 24.0)


def is_daytime(t_hrs: float, placement_hour: int) -> float:
    """Binary daytime indicator for side-face solar gating.

    Returns 1.0 between 6 AM and 6 PM local clock time, 0.0 otherwise.

    Rationale: the top BC uses env.solar_W_m2(t) which naturally zeros at night
    (the CW weather file contains measured irradiance). The side BC instead applies
    a flat F_vert = 0.5 projection to approximate the fraction of horizontal solar
    that reaches a vertical form face averaged over a day. Without a daytime gate
    this flat factor would smear daily solar energy into nighttime hours through an
    averaging bias. The 6-18 window is a rough daylight proxy for mid-latitude
    locations; the residual error from astronomical dawn/dusk (solar zenith > 75°)
    is small and is deferred to Sprint 2.

    Parameters
    ----------
    t_hrs : float        Hours since concrete placement.
    placement_hour : int Hour of day (0–23) when placement started.

    Returns
    -------
    float  1.0 during 6 AM ≤ local hour < 6 PM, else 0.0.
    """
    hour_of_day = (placement_hour + t_hrs) % 24.0
    return 1.0 if 6.0 <= hour_of_day < 18.0 else 0.0


def compute_T_gw_C(env) -> float:
    """Deep ground temperature from mean of daily max/min over the run window.

    CW Eq 44: T_gw_C = 0.83 * T_aat_C + 3.7.
    For Austin Jul: T_aat ≈ 83.5°F ≈ 28.6°C → T_gw ≈ 27.4°C.

    Parameters
    ----------
    env : duck-typed  Must have daily_max_F, daily_min_F (per-day lists).

    Returns
    -------
    float  Deep ground (water table) temperature in °C.
    """
    T_aat_F = (np.mean(env.daily_max_F) + np.mean(env.daily_min_F)) / 2.0
    T_aat_C = (T_aat_F - 32.0) * 5.0 / 9.0
    return 0.83 * T_aat_C + 3.7


@dataclass
class ConductionResult:
    """Output from solve_conduction_2d.

    All temperatures in °C, times in seconds.
    """

    t_s: np.ndarray           # shape (n_out,) — sampled time points
    T_field_C: np.ndarray     # shape (n_out, ny, nx)
    dt_inner_s: float         # actual inner time step used
    n_inner_steps: int        # total inner time steps executed
    n_output_samples: int     # len(t_s)
    peak_T_C: float           # max over all time and space
    min_T_C: float            # min over all time and space


@dataclass
class HydrationResult:
    """Output from solve_hydration_2d.

    All temperatures in °C, times in seconds, equivalent age in hours.
    """

    t_s: np.ndarray              # shape (n_out,) — sampled output times
    T_field_C: np.ndarray        # shape (n_out, ny, nx)
    alpha_field: np.ndarray      # shape (n_out, ny, nx) — degree of hydration
    t_e_field_hrs: np.ndarray    # shape (n_out, ny, nx) — equivalent age (hours)
    dt_inner_s: float            # nominal CFL-limited inner time step (seconds)
    n_inner_steps: int           # total inner steps executed
    n_output_samples: int        # len(t_s)
    peak_T_C: float              # max temperature over all time and space
    peak_T_location: tuple       # (j, i) grid indices of peak temperature
    peak_alpha: float            # max degree of hydration over all time and space
    centerline_T_C: np.ndarray   # shape (n_out, ny) — T at x = x_max (centerline column)
    # M3 diagnostic fields — None when boundary_mode is 'adiabatic' or 'dirichlet'
    top_flux_W_m2_history: np.ndarray | None = None  # shape (n_out, nx) surface heat loss
    T_amb_C_history: np.ndarray | None = None        # shape (n_out,) ambient air temperature
    # PR 2 solar diagnostics — both None when diagnostic_outputs=False or non-full_2d mode.
    # q_solar_history         : flux that ENTERS the concrete after blanket attenuation
    #                           = −α·G·f  where f = h_top_combined/h_forced_convection ≈ 0.047
    #                           Use this for the PR 4 total-top-flux panel; it is additive with
    #                           q_conv and q_lw (PR 3) which are all expressed at the concrete face.
    # q_solar_incident_history: raw absorbed flux at the blanket outer surface = −α·G
    #                           Use this for validation against measured/CW solar inputs.
    q_solar_history: np.ndarray | None = None           # (n_out, nx) W/m²; negative = heat in
    q_solar_incident_history: np.ndarray | None = None  # (n_out, nx) W/m²; negative = heat in
    # PR 3 LW diagnostics — all None when diagnostic_outputs=False or non-full_2d mode.
    q_LW_history: np.ndarray | None = None              # (n_out, nx) W/m², effective (positive=heat out)
    q_LW_incident_history: np.ndarray | None = None     # (n_out, nx) W/m², raw at blanket outer surface
    T_outer_C_history: np.ndarray | None = None         # (n_out, nx) °C, blanket outer-surface temperature
    # PR 4 top-face completion — None when diagnostic_outputs=False or non-full_2d mode.
    # q_top_total_history == top_flux_W_m2_history to numerical precision; both kept for compat.
    q_conv_history: np.ndarray | None = None            # (n_out, nx) W/m², convective flux at concrete top
    q_evap_history: np.ndarray | None = None            # (n_out, nx) W/m², evaporative flux at concrete top
    q_top_total_history: np.ndarray | None = None       # (n_out, nx) W/m², sum of all 4 top components
    # PR 4 side-face diagnostics — None when diagnostic_outputs=False or non-full_2d mode.
    # Shape (n_out, n_side_rows) where n_side_rows = iy_concrete_end - 1.
    q_side_solar_history: np.ndarray | None = None      # (n_out, n_side_rows) W/m², solar on form face
    q_side_LW_history: np.ndarray | None = None         # (n_out, n_side_rows) W/m², LW on form face
    q_side_conv_history: np.ndarray | None = None       # (n_out, n_side_rows) W/m², convection on form face
    q_side_total_history: np.ndarray | None = None      # (n_out, n_side_rows) W/m², total side flux
    # PR 5 (M6a) Sprint 2 plumbing — populated by PR 6's vertical-form
    # T_outer solve. Declared here but never written in PR 5 (remains None
    # in all code paths). Shape (n_out, n_side_rows) matching the sibling
    # q_side_*_history fields. n_side_rows = iy_concrete_end - 1.
    T_outer_form_C_history: np.ndarray | None = None    # (n_out, n_side_rows) °C


@dataclass
class Grid2D:
    """Half-mat 2D grid for the thermal engine.

    Vertex-centered: nodes sit at material interfaces.
    Node [j, i] is at (x[i], y[j]).  Shape convention: (ny, nx).
    """

    # Coordinate arrays (meters)
    x: np.ndarray  # shape (nx,), monotonically increasing
    y: np.ndarray  # shape (ny,), monotonically increasing (downward)

    # Material ID per cell. Shape (ny, nx). Values: 0=blanket, 1=concrete, 2=soil, 3=air/inactive
    material_id: np.ndarray

    # Boolean convenience masks, all shape (ny, nx)
    is_blanket: np.ndarray
    is_concrete: np.ndarray
    is_soil: np.ndarray
    is_air: np.ndarray    # material_id == 3; cells not solved by stencil

    # Index metadata for later milestones
    ix_concrete_start: int  # first x-index where material is concrete/blanket (= n_soil_x_ext)
    ix_concrete_end: int    # last x-index of concrete/blanket (= nx - 1, inclusive)
    iy_blanket_end: int     # y-index of blanket/concrete interface (= n_blanket; node belongs to concrete)
    iy_concrete_start: int  # first y-index of concrete (= n_blanket)
    iy_concrete_end: int    # last y-index of concrete (= n_blanket + n_concrete_y - 1, inclusive)
    iy_centerline: int      # last x-index (= nx - 1), where symmetry BC is applied in M4

    # Raw build params (for debugging / regeneration)
    nx: int
    ny: int
    dx: float   # uniform across x
    n_blanket: int
    n_concrete_x: int
    n_concrete_y: int
    n_soil_x_ext: int
    n_soil_y: int

    def concrete_slice(self) -> tuple[slice, slice]:
        """Returns (y_slice, x_slice) to extract concrete subdomain from a full field."""
        return (
            slice(self.iy_concrete_start, self.iy_concrete_end + 1),
            slice(self.ix_concrete_start, self.ix_concrete_end + 1),
        )


def build_grid_half_mat(
    width_ft: float,
    depth_ft: float,
    n_concrete_x: int = 21,
    n_concrete_y: int = 13,
    blanket_thickness_m: float = 0.02,
    n_blanket: int = 1,
    soil_depth_below_m: float = 3.0,
    n_soil_y: int = 15,
    n_soil_x_ext: int = 12,
    soil_ext_lateral_m: float | None = None,
) -> Grid2D:
    """Build a half-mat 2D grid for the thermal solver.

    The grid is vertex-centered.  Material boundaries align with node rows /
    columns so interface nodes are unambiguously on the boundary.

    Parameters
    ----------
    width_ft : float
        Full mat width in feet.  Half-width (centerline) is modelled.
    depth_ft : float
        Concrete depth in feet.
    n_concrete_x : int
        Number of nodes across the concrete half-width (vertex-centered,
        so n_concrete_x - 1 intervals span exactly W/2).
    n_concrete_y : int
        Number of nodes through the concrete depth.
    blanket_thickness_m : float
        Physical blanket thickness in metres.  Blanket gets n_blanket nodes.
    n_blanket : int
        Number of y-nodes allocated to the blanket layer.
    soil_depth_below_m : float
        Depth of soil domain below the concrete bottom.
    n_soil_y : int
        Number of y-nodes in the soil below concrete.
    n_soil_x_ext : int
        Number of soil-extension cells to the left of x=0.  Cell edges align
        with the uniform dx so the full x-axis has a single dx.
    soil_ext_lateral_m : float | None
        Optional target lateral extent.  If it disagrees with
        n_soil_x_ext * dx by more than 1 cm, a UserWarning is emitted and the
        exact-alignment value is used instead (dx stays uniform).

    Returns
    -------
    Grid2D
    """
    half_width_m = width_ft * FT_TO_M / 2.0
    depth_m = depth_ft * FT_TO_M

    # Uniform dx derived from concrete node spacing (vertex-centered: n-1 intervals)
    dx = half_width_m / (n_concrete_x - 1)

    # Soil-extension lateral extent aligned to exact integer multiples of dx
    L_ext_actual = n_soil_x_ext * dx

    if soil_ext_lateral_m is not None and abs(soil_ext_lateral_m - L_ext_actual) > 0.01:
        warnings.warn(
            f"soil_ext_lateral_m={soil_ext_lateral_m:.4f} m disagrees with "
            f"n_soil_x_ext * dx = {L_ext_actual:.4f} m by more than 1 cm. "
            f"Using exact-alignment value {L_ext_actual:.4f} m to keep dx uniform.",
            UserWarning,
            stacklevel=2,
        )

    nx = n_soil_x_ext + n_concrete_x
    ny = n_blanket + n_concrete_y + n_soil_y

    # --- x coordinates (uniform, monotonically increasing) ---
    x = np.linspace(-L_ext_actual, half_width_m, nx)

    # --- y coordinates (non-uniform: blanket << concrete/soil dy) ---
    # Blanket: n_blanket nodes from y=0 to y=blanket_thickness_m.
    y_blanket = np.linspace(0.0, blanket_thickness_m, n_blanket + 1)[:n_blanket]

    # Concrete: n_concrete_y nodes from blanket_thickness_m to blanket_thickness_m + depth_m.
    y_concrete = np.linspace(blanket_thickness_m, blanket_thickness_m + depth_m, n_concrete_y)

    # Soil below concrete: n_soil_y nodes, uniform dy_soil, starting one dy_soil below concrete bottom.
    dy_soil = soil_depth_below_m / n_soil_y
    y_concrete_bottom = blanket_thickness_m + depth_m
    y_soil = y_concrete_bottom + dy_soil * np.arange(1, n_soil_y + 1)

    y = np.concatenate([y_blanket, y_concrete, y_soil])

    # --- Index metadata ---
    ix_concrete_start = n_soil_x_ext
    ix_concrete_end = nx - 1
    iy_blanket_end = n_blanket  # interface row; y[n_blanket] = blanket_thickness_m
    iy_concrete_start = n_blanket
    iy_concrete_end = n_blanket + n_concrete_y - 1
    iy_centerline = nx - 1  # last x-index is the symmetry centerline

    # --- Material ID array: default fill = 2 (soil) ---
    material_id = np.full((ny, nx), 2, dtype=np.int8)

    # Blanket: rows 0..n_blanket-1, cols ix_concrete_start..nx-1
    material_id[:n_blanket, ix_concrete_start:] = 0

    # Concrete: rows iy_concrete_start..iy_concrete_end, cols ix_concrete_start..nx-1
    material_id[iy_concrete_start : iy_concrete_end + 1, ix_concrete_start:] = 1

    # Air/inactive: soil-extension strip (x < 0) for rows above the concrete-soil
    # interface. These cells are NOT solved by the stencil; they exist only so that
    # indices align. CW models air on the sides of the concrete slab, not soil.
    material_id[: iy_concrete_end + 1, :ix_concrete_start] = 3

    # --- Boolean masks ---
    is_blanket = material_id == 0
    is_concrete = material_id == 1
    is_soil = material_id == 2
    is_air = material_id == 3

    return Grid2D(
        x=x,
        y=y,
        material_id=material_id,
        is_blanket=is_blanket,
        is_concrete=is_concrete,
        is_soil=is_soil,
        is_air=is_air,
        ix_concrete_start=ix_concrete_start,
        ix_concrete_end=ix_concrete_end,
        iy_blanket_end=iy_blanket_end,
        iy_concrete_start=iy_concrete_start,
        iy_concrete_end=iy_concrete_end,
        iy_centerline=iy_centerline,
        nx=nx,
        ny=ny,
        dx=dx,
        n_blanket=n_blanket,
        n_concrete_x=n_concrete_x,
        n_concrete_y=n_concrete_y,
        n_soil_x_ext=n_soil_x_ext,
        n_soil_y=n_soil_y,
    )


def build_grid_rectangular(
    Lx_m: float,
    Ly_m: float,
    nx: int,
    ny: int,
) -> Grid2D:
    """Build a uniform rectangular grid with all nodes labelled as concrete.

    Used for M1 analytical validation and future pavement scenarios.
    x in [0, Lx_m], y in [0, Ly_m], uniform spacing on both axes.

    Parameters
    ----------
    Lx_m : float
        Domain width in metres (x-extent).
    Ly_m : float
        Domain depth in metres (y-extent).
    nx : int
        Number of nodes in x (>= 3).
    ny : int
        Number of nodes in y (>= 3).

    Returns
    -------
    Grid2D
    """
    if Lx_m <= 0.0 or Ly_m <= 0.0:
        raise ValueError(f"Lx_m and Ly_m must be positive, got {Lx_m}, {Ly_m}")
    if nx < 3 or ny < 3:
        raise ValueError(f"nx and ny must be >= 3, got {nx}, {ny}")

    x = np.linspace(0.0, Lx_m, nx)
    y = np.linspace(0.0, Ly_m, ny)
    dx = Lx_m / (nx - 1)

    material_id = np.ones((ny, nx), dtype=np.int8)
    is_concrete = np.ones((ny, nx), dtype=bool)
    is_blanket = np.zeros((ny, nx), dtype=bool)
    is_soil = np.zeros((ny, nx), dtype=bool)
    is_air = np.zeros((ny, nx), dtype=bool)

    return Grid2D(
        x=x,
        y=y,
        material_id=material_id,
        is_blanket=is_blanket,
        is_concrete=is_concrete,
        is_soil=is_soil,
        is_air=is_air,
        ix_concrete_start=0,
        ix_concrete_end=nx - 1,
        iy_blanket_end=-1,
        iy_concrete_start=0,
        iy_concrete_end=ny - 1,
        iy_centerline=ny // 2,
        nx=nx,
        ny=ny,
        dx=dx,
        n_blanket=0,
        n_concrete_x=nx,
        n_concrete_y=ny,
        n_soil_x_ext=0,
        n_soil_y=0,
    )


def _stencil_step(
    T: np.ndarray,
    T_new: np.ndarray,
    k_x_face: np.ndarray,
    k_y_face: np.ndarray,
    rho_cp: np.ndarray,
    Q: np.ndarray,
    inv_dxsq: float,
    dy_plus: np.ndarray,
    dy_minus: np.ndarray,
    y_coef: np.ndarray,
    dt: float,
) -> None:
    """Vectorized 5-point FD interior update, shared by both solvers.

    Updates T_new[1:-1, 1:-1] in place. Boundary rows/columns are NOT touched;
    the caller must apply BCs after this call.

    Parameters
    ----------
    T : (ny, nx)  current temperature field °C.
    T_new : (ny, nx)  output buffer — interior updated, edges unchanged.
    k_x_face : (ny, nx-1)  harmonic-mean k on x-directed faces W/(m·K).
    k_y_face : (ny-1, nx)  harmonic-mean k on y-directed faces W/(m·K).
    rho_cp : (ny, nx)  product ρ·Cp J/(m³·K).
    Q : (ny, nx)  volumetric heat source W/m³.
    inv_dxsq : float  1/dx² m⁻².
    dy_plus : (ny-2, 1)  y[j+1]-y[j] for interior rows j=1..ny-2.
    dy_minus : (ny-2, 1)  y[j]-y[j-1].
    y_coef : (ny-2, 1)  2/(dy_plus + dy_minus).
    dt : float  time step seconds.
    """
    Tc  = T[1:-1, 1:-1]
    Tip = T[1:-1, 2:  ]
    Tim = T[1:-1,  :-2]
    Tjp = T[2:  , 1:-1]
    Tjm = T[ :-2, 1:-1]

    # Face-k arrays sliced to interior cell faces
    kxp = k_x_face[1:-1, 1:  ]   # x+ face of interior cells
    kxm = k_x_face[1:-1,  :-1]   # x- face
    kyp = k_y_face[1:,    1:-1]  # y+ face
    kym = k_y_face[ :-1,  1:-1]  # y- face

    flux_x = (kxp * (Tip - Tc) - kxm * (Tc - Tim)) * inv_dxsq
    flux_y = y_coef * (kyp * (Tjp - Tc) / dy_plus - kym * (Tc - Tjm) / dy_minus)

    T_new[1:-1, 1:-1] = (
        Tc + dt * (flux_x + flux_y + Q[1:-1, 1:-1]) / rho_cp[1:-1, 1:-1]
    )


def solve_conduction_2d(
    grid: Grid2D,
    k_W_m_K: float,
    rho_kg_m3: float,
    Cp_J_kg_K: float,
    T_initial_C: np.ndarray,
    T_boundary_C: float,
    duration_s: float,
    output_interval_s: float = 300.0,
    cfl_safety: float = 0.4,
) -> ConductionResult:
    """Explicit 5-point 2D conduction solver with constant properties.

    Solves  ρ Cp ∂T/∂t = k (∂²T/∂x² + ∂²T/∂y²)  on the Grid2D domain
    with fixed Dirichlet BCs on all four sides.

    Stencil (D4, non-uniform dy):
        T_new[j,i] = T[j,i] + dt·α·(
            (T[j,i+1] - 2T[j,i] + T[j,i-1]) / dx²
          + (2/(dy+ + dy-)) · ((T[j+1,i]-T[j,i])/dy+ - (T[j,i]-T[j-1,i])/dy-)
        )

    Parameters
    ----------
    grid : Grid2D
    k_W_m_K : float
        Thermal conductivity W/(m·K).
    rho_kg_m3 : float
        Density kg/m³.
    Cp_J_kg_K : float
        Specific heat J/(kg·K).
    T_initial_C : np.ndarray
        Shape (ny, nx), initial temperature field in °C.
    T_boundary_C : float
        Dirichlet BC value applied to all four edges every step.
    duration_s : float
        Total simulation duration in seconds.
    output_interval_s : float
        Time between output samples in seconds.
    cfl_safety : float
        Multiplier on CFL-limited dt; must be <= 1.0.

    Returns
    -------
    ConductionResult
    """
    if cfl_safety > 1.0:
        raise ValueError(
            f"cfl_safety={cfl_safety} > 1.0 violates explicit stability; "
            "use cfl_safety <= 1.0"
        )

    alpha = k_W_m_K / (rho_kg_m3 * Cp_J_kg_K)
    dx = grid.dx
    dy_all = np.diff(grid.y)          # (ny-1,)
    dy_min = float(dy_all.min())

    dt_cfl = 0.5 / (alpha * (1.0 / dx**2 + 1.0 / dy_min**2))
    dt = cfl_safety * dt_cfl

    # per-interior-row spacings broadcast over x automatically
    dy_plus  = dy_all[1:].reshape(-1, 1)    # (ny-2, 1): y[j+1]-y[j] for j=1..ny-2
    dy_minus = dy_all[:-1].reshape(-1, 1)   # (ny-2, 1): y[j]-y[j-1]
    inv_dxsq = 1.0 / dx**2
    y_coef   = 2.0 / (dy_plus + dy_minus)   # (ny-2, 1)

    # Build uniform face-k and source arrays for the constant-property case.
    # Harmonic mean of equal values equals the value, so this is algebraically
    # identical to the M1 scalar k_face stencil.
    ny, nx = grid.ny, grid.nx
    k_cell   = np.full((ny, nx), k_W_m_K)
    rho_cp   = np.full((ny, nx), rho_kg_m3 * Cp_J_kg_K)
    Q_zero   = np.zeros((ny, nx))
    k_x_face = 2.0 * k_cell[:, :-1] * k_cell[:, 1:] / (k_cell[:, :-1] + k_cell[:, 1:])
    k_y_face = 2.0 * k_cell[:-1, :] * k_cell[1:, :] / (k_cell[:-1, :] + k_cell[1:, :])

    # Determine output sample times (include t=0)
    n_samples = int(round(duration_s / output_interval_s)) + 1
    t_samples = np.linspace(0.0, duration_s, n_samples)

    T_out = np.empty((n_samples, grid.ny, grid.nx), dtype=np.float64)
    t_out = np.empty(n_samples, dtype=np.float64)

    T = T_initial_C.astype(np.float64, copy=True)
    T_new = np.empty_like(T)

    # Apply BC to initial field so boundary is consistent from the start
    T[0, :]  = T_boundary_C
    T[-1, :] = T_boundary_C
    T[:, 0]  = T_boundary_C
    T[:, -1] = T_boundary_C

    # Store t=0 sample
    T_out[0] = T
    t_out[0] = 0.0
    next_sample = 1

    n_steps = 0
    sim_time = 0.0

    while sim_time < duration_s - 1e-12:
        # Don't overshoot duration
        dt_step = min(dt, duration_s - sim_time)

        _stencil_step(T, T_new, k_x_face, k_y_face, rho_cp, Q_zero,
                      inv_dxsq, dy_plus, dy_minus, y_coef, dt_step)

        T_new[0, :]  = T_boundary_C
        T_new[-1, :] = T_boundary_C
        T_new[:, 0]  = T_boundary_C
        T_new[:, -1] = T_boundary_C

        T, T_new = T_new, T
        sim_time += dt_step
        n_steps += 1

        # Capture output samples whose scheduled time has been reached
        while next_sample < n_samples and sim_time >= t_samples[next_sample] - dt * 0.5:
            T_out[next_sample] = T
            t_out[next_sample] = sim_time
            next_sample += 1

    # Ensure all samples filled (handles edge case where duration is exact)
    while next_sample < n_samples:
        T_out[next_sample] = T
        t_out[next_sample] = sim_time
        next_sample += 1

    return ConductionResult(
        t_s=t_out,
        T_field_C=T_out,
        dt_inner_s=dt,
        n_inner_steps=n_steps,
        n_output_samples=n_samples,
        peak_T_C=float(T_out.max()),
        min_T_C=float(T_out.min()),
    )


def analytical_square_slab(
    x: np.ndarray,
    y: np.ndarray,
    t: float,
    Lx: float,
    Ly: float,
    alpha: float,
    T_init: float,
    T_bc: float,
    n_terms: int = 11,
) -> np.ndarray:
    """Analytical solution for 2D diffusion in a rectangle with Dirichlet BCs.

    Solves ∂T/∂t = α(∂²T/∂x² + ∂²T/∂y²) on [0,Lx]×[0,Ly] with
    T=T_bc on all four sides and T(x,y,0)=T_init (uniform).

    θ(x,y,t) = θ₀·(16/π²)·Σ_m Σ_n sin(mπx/Lx)·sin(nπy/Ly)/(m·n)
                           ·exp(-α·π²·t·(m²/Lx²+n²/Ly²))
    where m, n = 1, 3, 5, … (odd integers only), θ = T - T_bc.

    Parameters
    ----------
    x : np.ndarray
        Shape (nx,), node x-coordinates (domain should start at 0).
    y : np.ndarray
        Shape (ny,), node y-coordinates (domain should start at 0).
    t : float
        Evaluation time in seconds.
    Lx, Ly : float
        Domain extents in metres.
    alpha : float
        Thermal diffusivity m²/s.
    T_init : float
        Uniform initial temperature °C.
    T_bc : float
        Dirichlet boundary temperature °C.
    n_terms : int
        Number of odd Fourier modes (1, 3, …, 2*n_terms-1).

    Returns
    -------
    np.ndarray
        Shape (ny, nx), temperature field in °C.
    """
    odd = 2 * np.arange(n_terms) + 1   # [1, 3, 5, ...]
    theta0 = T_init - T_bc
    prefactor = theta0 * (16.0 / np.pi**2)

    X = x[np.newaxis, :]    # (1, nx)
    Y = y[:, np.newaxis]    # (ny, 1)

    theta = np.zeros((y.size, x.size), dtype=np.float64)
    for m in odd:
        sx = np.sin(m * np.pi * X / Lx)         # (1, nx)
        for n in odd:
            sy = np.sin(n * np.pi * Y / Ly)     # (ny, 1)
            decay = np.exp(-alpha * np.pi**2 * t * (m**2 / Lx**2 + n**2 / Ly**2))
            theta += sx * sy * (decay / (m * n))

    return T_bc + prefactor * theta


def solve_hydration_2d(
    grid: Grid2D,
    mix,
    T_initial_C: np.ndarray,
    T_boundary_C: float = 20.0,
    duration_s: float = 0.0,
    output_interval_s: float = 300.0,
    cfl_safety: float = 0.4,
    boundary_mode: str = "adiabatic",
    soil_props: dict | None = None,
    blanket_R_value: float = 5.67,
    blanket_thickness_m: float = 0.02,
    # M3 parameters — required for 'top_bc_only' and 'v2_equivalent' modes:
    environment=None,
    construction=None,
    T_ground_deep_C: float | None = None,
    # Debug flag: force blanket-as-pure-R in non-full_2d modes for ablation
    # tests.  For full_2d this is always True (permanent M4 architecture).
    skip_blanket_node: bool = False,
    # Gate detailed per-sample flux histories.  Set False for memory-sensitive
    # long runs.  PR 2 gates q_solar_history only; PR 4 will extend this gate
    # to top_flux_W_m2_history and T_amb_C_history.
    diagnostic_outputs: bool = True,
) -> HydrationResult:
    """Explicit 2D conduction solver with hydration heat and variable properties.

    Solves  ρ Cp ∂T/∂t = ∇·(k∇T) + Q_hyd  with the Schindler-Folliard
    hydration model (equivalent age, degree of hydration α) and CW Eqs 23-25
    for variable k(α) and Cp(α, T).

    Boundary modes
    --------------
    'adiabatic'    — zero-flux ghost-node on all four edges (M2 default).
    'dirichlet'    — fixed T = T_boundary_C on all four edges.
    'top_bc_only'  — convection+blanket+Menzel on top; adiabatic on other three.
    'v2_equivalent'— convection+blanket+Menzel on top; Dirichlet T_gw on bottom;
                     adiabatic sides.  Matches v2 1D BC topology.

    Parameters
    ----------
    grid : Grid2D
    mix : duck-typed CWMixDesign
    T_initial_C : np.ndarray  shape (ny, nx), initial temperature °C.
    T_boundary_C : float      Dirichlet value for 'dirichlet' mode (°C).
    duration_s : float        Total simulation duration in seconds.
    output_interval_s : float Target output sample interval in seconds.
    cfl_safety : float        CFL multiplier, must be <= 1.0.
    boundary_mode : str       See above.
    soil_props : dict | None  Override SOIL_PROPERTIES_2D entries.
    blanket_R_value : float   Blanket R-value hr·ft²·°F/BTU (used as fallback
                               when construction is None).
    blanket_thickness_m : float  Blanket thickness metres.
    environment : duck-typed  Required for top_bc modes.  Must have:
                               daily_max_F, daily_min_F, cw_ave_max_wind_m_s,
                               cw_ave_max_RH_pct, cw_ave_min_RH_pct.
    construction : duck-typed Required for top_bc modes.  Must have:
                               blanket_R_value, top_cure_blanket_time_hrs,
                               placement_hour.
    T_ground_deep_C : float | None  Deep ground Dirichlet for 'v2_equivalent'.
                               If None, computed from env via CW Eq 44.

    Returns
    -------
    HydrationResult
    """
    _TOP_BC_MODES = ("top_bc_only", "v2_equivalent", "full_2d")
    if cfl_safety > 1.0:
        raise ValueError(
            f"cfl_safety={cfl_safety} > 1.0 violates explicit stability; "
            "use cfl_safety <= 1.0"
        )
    if boundary_mode not in {"adiabatic", "dirichlet"} | set(_TOP_BC_MODES):
        raise ValueError(
            f"boundary_mode must be one of 'adiabatic', 'dirichlet', "
            f"'top_bc_only', 'v2_equivalent', 'full_2d'; got {boundary_mode!r}"
        )
    if boundary_mode in _TOP_BC_MODES and (environment is None or construction is None):
        raise ValueError(
            f"boundary_mode={boundary_mode!r} requires environment and construction"
        )

    # ------------------------------------------------------------------ #
    # A. Unit conversions (imperial → SI, from mix)                        #
    # ------------------------------------------------------------------ #
    Hu     = mix.Hu_J_kg                                          # J/kg_cement
    tau    = mix.tau_hrs                                          # hrs
    beta_h = mix.beta                                             # dimensionless
    au     = mix.alpha_u                                          # dimensionless
    Ea     = mix.activation_energy_J_mol                         # J/mol
    rho_c  = mix.concrete_density_lb_ft3 * LB_FT3_TO_KG_M3      # kg/m³
    Cc     = mix.total_cementitious_lb_yd3 * LB_YD3_TO_KG_M3    # kg_cement/m³_concrete
    k_uc   = mix.thermal_conductivity_BTU_hr_ft_F * BTU_HR_FT_F_TO_W_M_K  # W/(m·K)
    Ca     = mix.aggregate_Cp_BTU_lb_F * BTU_LB_F_TO_J_KG_K     # J/(kg·K)
    Wc     = mix.cement_type_I_II_lb_yd3 * LB_YD3_TO_KG_M3      # kg/m³
    Ww     = mix.water_lb_yd3 * LB_YD3_TO_KG_M3                 # kg/m³
    Wa     = (mix.coarse_agg_lb_yd3 + mix.fine_agg_lb_yd3) * LB_YD3_TO_KG_M3  # kg/m³

    # ------------------------------------------------------------------ #
    # A2. M3/M4 top-BC pre-computation                                   #
    # ------------------------------------------------------------------ #
    _use_top_bc = boundary_mode in _TOP_BC_MODES
    if _use_top_bc:
        _r_blanket = getattr(construction, "blanket_R_value", blanket_R_value)
        _R_blanket_SI = _r_blanket * R_IMP_TO_R_SI  # m²·K/W
        # PR 2: h_forced_convection (no +5.5 LW term) is the production coefficient.
        # _h_top_combined = forced_conv + blanket R in series → used by full_2d.
        # _h_top_legacy   = legacy conv (+5.5 LW) + blanket R → used by top_bc_only /
        #                    v2_equivalent to preserve pre-PR-2 numeric behaviour.
        _h_conv_top     = h_forced_convection(environment.cw_ave_max_wind_m_s)  # kept for non-blanket branch
        _h_top_combined = 1.0 / (1.0 / h_forced_convection(environment.cw_ave_max_wind_m_s) + _R_blanket_SI)
        _h_top_legacy   = 1.0 / (1.0 / _h_convective_legacy(environment.cw_ave_max_wind_m_s) + _R_blanket_SI)
        # Solar absorptivity and LW emissivity for the top blanket outer surface (PR 2/3).
        _alpha_sol_top  = getattr(construction, "solar_absorptivity_top", SOLAR_ABSORPTIVITY_DEFAULT)
        _emis_top       = getattr(construction, "emissivity_top", EMISSIVITY_DEFAULT)
        _wind_eff = environment.cw_ave_max_wind_m_s * 0.4  # v2 derating
        _RH_frac = 0.005 * (environment.cw_ave_max_RH_pct + environment.cw_ave_min_RH_pct)
        _cure_time = getattr(construction, "top_cure_blanket_time_hrs", 2.0)
        _placement_hour = construction.placement_hour
        _T_gw_C = (
            T_ground_deep_C if T_ground_deep_C is not None
            else compute_T_gw_C(environment)
        )
        # Columns at top row that are air/inactive (not solved by stencil)
        _top_row_is_air = grid.is_air[0, :]
        # Legacy mask for top_bc_only / v2_equivalent (soil beside mat)
        _top_row_is_soil = (grid.material_id[0, :] == 2)
        # M4 side BC: h_side for the form face (PR 4: pure forced convection)
        _h_side = h_side_convective(environment.cw_ave_max_wind_m_s, _r_blanket)
        # PR 6: pure convective h for the T_outer_form Newton solve (no R_form baked in)
        _h_conv_vert = h_forced_convection(environment.cw_ave_max_wind_m_s)
        _ix_cs = grid.ix_concrete_start     # column index of concrete form edge
        # PR 4 side-face radiation properties
        _alpha_sol_side = getattr(construction, "solar_absorptivity_side", SOLAR_ABSORPTIVITY_DEFAULT)
        _emis_side      = getattr(construction, "emissivity_side", EMISSIVITY_DEFAULT)
        _F_vert         = getattr(construction, "vertical_solar_factor", VERTICAL_SOLAR_FACTOR_DEFAULT)
        # _F_vert = 0.0 by default (CWConstruction.vertical_solar_factor=0.0); Sprint 2 will
        # calibrate using ACI 207.2R Eq 27 orientation model. F_vert=0.5 (physics) causes the
        # side-face minimum cell to warm, reducing the gradient metric below the ±2.0°F gate.

    # ------------------------------------------------------------------ #
    # B. Soil and blanket material properties                              #
    # ------------------------------------------------------------------ #
    soil = {**SOIL_PROPERTIES_2D, **(soil_props or {})}
    R_SI = blanket_R_value * R_IMP_TO_R_SI                       # m²·K/W
    k_blanket = blanket_thickness_m / R_SI                       # W/(m·K)

    # ------------------------------------------------------------------ #
    # C. Per-cell constant property arrays (shape ny, nx)                  #
    # ------------------------------------------------------------------ #
    ny, nx = grid.ny, grid.nx
    rho_cell = np.empty((ny, nx), dtype=np.float64)
    rho_cell[grid.is_concrete] = rho_c
    rho_cell[grid.is_soil]     = soil['rho']
    rho_cell[grid.is_blanket]  = BLANKET_PROPERTIES_2D['rho']
    rho_cell[grid.is_air]      = 1.0   # placeholder; air cells are never updated

    # Cp for soil/blanket constant; concrete recomputed each step
    Cp_cell = np.empty((ny, nx), dtype=np.float64)
    Cp_cell[grid.is_soil]    = soil['Cp']
    Cp_cell[grid.is_blanket] = BLANKET_PROPERTIES_2D['Cp']
    Cp_cell[grid.is_air]     = 1.0   # placeholder
    # Concrete Cp initialised at α=0:
    Cp_cell[grid.is_concrete] = specific_heat_variable(
        Wc, Wa, Ww, Ca,
        np.zeros(grid.is_concrete.sum()),
        T_initial_C[grid.is_concrete],
        rho_c,
    )

    # k_cell: concrete initialised at α=0 (maximum k = 1.33·k_uc)
    k_cell = np.empty((ny, nx), dtype=np.float64)
    k_cell[grid.is_concrete] = thermal_conductivity_variable(k_uc, 0.0)
    k_cell[grid.is_soil]     = soil['k']
    k_cell[grid.is_blanket]  = k_blanket
    k_cell[grid.is_air]      = 0.0   # zero k so face conductances toward air are also zero

    # ------------------------------------------------------------------ #
    # M4 architecture: blanket as pure-R in full_2d mode                   #
    # In full_2d the blanket row carries no thermal mass — it is modelled   #
    # only as a thermal resistance (R_blanket) in the top BC coefficient    #
    # applied at j=iy_concrete_start.  This avoids the blanket-as-node      #
    # timing error (blanket absorbs ambient heat during hot days and drives  #
    # it into the 60°F concrete, shifting peak T ~20 hr early vs CW).       #
    # The skip_blanket_node flag activates the same behaviour for ablation   #
    # tests in non-full_2d modes.                                            #
    # ------------------------------------------------------------------ #
    _use_pure_r_blanket = (boundary_mode == "full_2d") or skip_blanket_node
    if _use_pure_r_blanket:
        k_cell[grid.is_blanket]    = 0.0
        rho_cell[grid.is_blanket]  = 1.0
        Cp_cell[grid.is_blanket]   = 1.0

    # ------------------------------------------------------------------ #
    # D. Hydration state                                                   #
    # ------------------------------------------------------------------ #
    te    = np.full((ny, nx), 0.01, dtype=np.float64)  # equivalent age (hrs)
    alpha = np.zeros((ny, nx), dtype=np.float64)
    # Non-concrete cells stay at te=0.01, alpha=0 permanently.

    # ------------------------------------------------------------------ #
    # E. Grid geometry for stencil                                         #
    # ------------------------------------------------------------------ #
    dx     = grid.dx
    dy_all = np.diff(grid.y)
    dy_min = float(dy_all.min())
    dy_plus  = dy_all[1:].reshape(-1, 1)
    dy_minus = dy_all[:-1].reshape(-1, 1)
    inv_dxsq = 1.0 / dx**2
    y_coef   = 2.0 / (dy_plus + dy_minus)

    # ------------------------------------------------------------------ #
    # F. CFL bound (conservative: max k, min Cp across all materials)     #
    # ------------------------------------------------------------------ #
    alpha_diff_c   = 1.33 * k_uc / (rho_c * 800.0)            # concrete worst-case
    alpha_diff_soil = soil['k'] / (soil['rho'] * soil['Cp'])
    alpha_diff_blk  = k_blanket / (BLANKET_PROPERTIES_2D['rho'] * BLANKET_PROPERTIES_2D['Cp'])
    alpha_diff_max  = max(alpha_diff_c, alpha_diff_soil, alpha_diff_blk)
    dt_cfl   = 0.5 / (alpha_diff_max * (1.0 / dx**2 + 1.0 / dy_min**2))
    dt_inner = cfl_safety * dt_cfl

    # ------------------------------------------------------------------ #
    # G. Output schedule (exact sample times)                              #
    # ------------------------------------------------------------------ #
    n_samples = int(round(duration_s / output_interval_s)) + 1
    t_samples = np.linspace(0.0, duration_s, n_samples)
    # Clamp last sample to exactly duration_s
    t_samples[-1] = duration_s

    T_out     = np.empty((n_samples, ny, nx), dtype=np.float64)
    alpha_out = np.empty((n_samples, ny, nx), dtype=np.float64)
    te_out    = np.empty((n_samples, ny, nx), dtype=np.float64)
    t_out     = np.empty(n_samples, dtype=np.float64)

    if _use_top_bc:
        top_flux_out = np.empty((n_samples, nx), dtype=np.float64)
        T_amb_out    = np.empty(n_samples, dtype=np.float64)
        if diagnostic_outputs and boundary_mode == "full_2d":
            _n_side_rows = grid.iy_concrete_end - 1   # j=1..iy_concrete_end-1
            q_solar_out          = np.empty((n_samples, nx), dtype=np.float64)
            q_solar_incident_out = np.empty((n_samples, nx), dtype=np.float64)
            q_lw_out             = np.empty((n_samples, nx), dtype=np.float64)
            q_lw_incident_out    = np.empty((n_samples, nx), dtype=np.float64)
            T_outer_out          = np.empty((n_samples, nx), dtype=np.float64)
            # PR 4: top-face completion
            q_conv_out           = np.empty((n_samples, nx), dtype=np.float64)
            q_evap_out           = np.empty((n_samples, nx), dtype=np.float64)
            q_top_total_out      = np.empty((n_samples, nx), dtype=np.float64)
            # PR 4: side-face diagnostics
            q_side_solar_out     = np.empty((n_samples, _n_side_rows), dtype=np.float64)
            q_side_lw_out        = np.empty((n_samples, _n_side_rows), dtype=np.float64)
            q_side_conv_out      = np.empty((n_samples, _n_side_rows), dtype=np.float64)
            q_side_total_out     = np.empty((n_samples, _n_side_rows), dtype=np.float64)
            # PR 6: vertical-form T_outer history
            T_outer_form_out     = np.empty((n_samples, _n_side_rows), dtype=np.float64)
        else:
            q_solar_out          = None
            q_solar_incident_out = None
            q_lw_out             = None
            q_lw_incident_out    = None
            T_outer_out          = None
            q_conv_out           = None
            q_evap_out           = None
            q_top_total_out      = None
            q_side_solar_out     = None
            q_side_lw_out        = None
            q_side_conv_out      = None
            q_side_total_out     = None
            T_outer_form_out     = None
    else:
        top_flux_out         = None
        T_amb_out            = None
        q_solar_out          = None
        q_solar_incident_out = None
        q_lw_out             = None
        q_lw_incident_out    = None
        T_outer_out          = None
        q_conv_out           = None
        q_evap_out           = None
        q_top_total_out      = None
        q_side_solar_out     = None
        q_side_lw_out        = None
        q_side_conv_out      = None
        q_side_total_out     = None
        T_outer_form_out     = None

    # ------------------------------------------------------------------ #
    # H. Initial conditions                                                #
    # ------------------------------------------------------------------ #
    T     = T_initial_C.astype(np.float64, copy=True)
    T_new = np.empty_like(T)

    T_out[0]     = T
    alpha_out[0] = alpha
    te_out[0]    = te
    t_out[0]     = 0.0
    if _use_top_bc:
        _T_amb_F_0 = ambient_temp_F(0.0, environment, _placement_hour)
        _T_amb_C_0 = (_T_amb_F_0 - 32.0) * 5.0 / 9.0
        T_amb_out[0]    = _T_amb_C_0
        top_flux_out[0] = np.zeros(nx)   # no flux at t=0 (no stencil yet)
        if q_solar_out is not None:
            q_solar_out[0]          = np.zeros(nx)  # no solar at t=0 (no stencil yet)
            q_solar_incident_out[0] = np.zeros(nx)
            q_lw_out[0]             = np.zeros(nx)
            q_lw_incident_out[0]    = np.zeros(nx)
            T_outer_out[0]          = np.zeros(nx)
            q_conv_out[0]           = np.zeros(nx)
            q_evap_out[0]           = np.zeros(nx)
            q_top_total_out[0]      = np.zeros(nx)
            q_side_solar_out[0]     = np.zeros(_n_side_rows)
            q_side_lw_out[0]        = np.zeros(_n_side_rows)
            q_side_conv_out[0]      = np.zeros(_n_side_rows)
            q_side_total_out[0]     = np.zeros(_n_side_rows)
            T_outer_form_out[0]     = T[1:grid.iy_concrete_end, _ix_cs].copy()  # initial T_outer ≈ T_conc
    next_idx = 1

    n_steps = 0
    t = 0.0

    # ------------------------------------------------------------------ #
    # I. Main time loop                                                    #
    # ------------------------------------------------------------------ #
    while next_idx < n_samples:
        dt_step = dt_inner
        if t + dt_step > t_samples[next_idx] + 1e-10:
            dt_step = t_samples[next_idx] - t

        # -- Hydration update on concrete cells (forward-Euler in time) --
        T_K = T + 273.15
        af = arrhenius_vec(T_K[grid.is_concrete], Ea)
        te[grid.is_concrete] += (dt_step / 3600.0) * af

        alpha[grid.is_concrete] = hydration_alpha_vec(
            te[grid.is_concrete], tau, beta_h, au
        )
        da_dte = hydration_rate_vec(te[grid.is_concrete], tau, beta_h, au)
        da_dt_per_s = da_dte * af / 3600.0   # [1/s] (chain rule + time unit)

        Q = np.zeros((ny, nx), dtype=np.float64)
        Q[grid.is_concrete] = Hu * Cc * da_dt_per_s   # W/m³

        # -- Refresh variable k and Cp on concrete cells --
        k_cell[grid.is_concrete] = thermal_conductivity_variable(
            k_uc, alpha[grid.is_concrete]
        )
        Cp_cell[grid.is_concrete] = specific_heat_variable(
            Wc, Wa, Ww, Ca,
            alpha[grid.is_concrete],
            T[grid.is_concrete],   # °C
            rho_c,
        )

        # -- Harmonic-mean interface conductivities --
        # Safe division: when both cells have k=0 (air-air face), denom=0 → 0/0.
        # np.where prevents the NaN; harmonic mean of (0,0) is correctly 0.
        _denom_x = k_cell[:, :-1] + k_cell[:, 1:]
        k_x_face = np.where(
            _denom_x == 0.0, 0.0,
            2.0 * k_cell[:, :-1] * k_cell[:, 1:] / _denom_x,
        )
        _denom_y = k_cell[:-1, :] + k_cell[1:, :]
        k_y_face = np.where(
            _denom_y == 0.0, 0.0,
            2.0 * k_cell[:-1, :] * k_cell[1:, :] / _denom_y,
        )
        # Air cell k=0 ensures all faces touching air are already 0 from the
        # formula above (harmonic mean of (0, k) = 0 for any k>=0). No
        # additional zeroing needed; this comment kept for clarity.

        rho_cp = rho_cell * Cp_cell

        # -- Interior FD update --
        _stencil_step(
            T, T_new, k_x_face, k_y_face, rho_cp, Q,
            inv_dxsq, dy_plus, dy_minus, y_coef, dt_step,
        )

        # Diagnostic placeholders (defined for all modes; populated by top-BC paths)
        _T_amb_C       = 0.0
        _q_top         = np.zeros(nx)
        _q_solar_array = np.zeros(nx)  # PR 2 solar flux; zero unless full_2d with solar

        # -- Boundary conditions --
        if boundary_mode == "adiabatic":
            T_new[0, :]  = T_new[1, :]
            T_new[-1, :] = T_new[-2, :]
            T_new[:, 0]  = T_new[:, 1]
            T_new[:, -1] = T_new[:, -2]

        elif boundary_mode == "dirichlet":
            T_new[0, :]  = T_boundary_C
            T_new[-1, :] = T_boundary_C
            T_new[:, 0]  = T_boundary_C
            T_new[:, -1] = T_boundary_C

        elif boundary_mode == "full_2d":
            # ---- M4 production mode ----
            t_hrs = t / 3600.0
            _T_amb_F = ambient_temp_F(t_hrs, environment, _placement_hour)
            _T_amb_C = (_T_amb_F - 32.0) * 5.0 / 9.0

            # (a) Top surface BC
            if _use_pure_r_blanket:
                # M4 full_2d: blanket is pure-R.  Apply combined BC at concrete top (j=1)
                # using _h_top_combined (h_forced_convection + blanket R in series).
                # PR 2 also adds solar absorption q_sw = -alpha_sol * G(t).
                _j_top = grid.iy_concrete_start   # j=1
                _dy_conc_top = float(grid.y[_j_top + 1] - grid.y[_j_top])
                _dy_half_c   = _dy_conc_top / 2.0
                _T_surf_c    = T[_j_top, :]
                _q_conv_c    = _h_top_combined * (_T_surf_c - _T_amb_C)
                _q_evap_c    = np.zeros(nx)
                if t_hrs < _cure_time:
                    for _i in range(grid.ix_concrete_start, nx):
                        _q_evap_c[_i] = menzel_evaporation(
                            t_hrs, float(_T_surf_c[_i]), _T_amb_C,
                            _RH_frac, _wind_eff
                        )
                # ── Blanket-surface quasi-steady energy balance ─────────────────────
                # The cure blanket is pure thermal resistance R_blanket [m²·K/W].
                # Solar radiation (and PR 3 longwave) is absorbed at the OUTER surface of
                # the blanket, not at the concrete.  The outer surface has no thermal mass
                # (pure-R model), so it is in quasi-steady state every time step.
                #
                # Setting net flux = 0 on the outer blanket surface (sign convention:
                # positive = leaving the surface):
                #
                #   h_conv·(T_outer − T_amb)        [convection to air]
                # + (T_outer − T_concrete)/R_blanket [conduction to concrete]
                # − α·G                              [absorbed solar, heat in]
                # = 0
                #
                # Solving for the downward flux into concrete:
                #
                #   q_down = (T_outer − T_concrete)/R_blanket
                #          = [α·G − h_conv·(T_concrete − T_amb)] / (h_conv·R_blanket + 1)
                #
                # Since h_top_combined = h_conv/(h_conv·R_blanket + 1), this simplifies to:
                #
                #   q_down = (h_top_combined/h_conv)·α·G
                #            − h_top_combined·(T_concrete − T_amb)
                #
                # The second term is already _q_conv_c.  The first term is the solar
                # contribution reaching the concrete, attenuated by factor
                #
                #   f = h_top_combined / h_forced_convection
                #
                # MIX-01 Austin: h_conv ≈ 20.3, R_blanket ≈ 0.998, h_top ≈ 0.954 W/(m²·K)
                #   f ≈ 0.047  →  ~4.7% of α·G reaches the concrete.
                #   Peak solar 837 W/m² → α·G ≈ 544 W/m² incident, ~25 W/m² to concrete.
                #
                # ============================================================
                # LONGWAVE RADIATION (PR 3) — IMPLEMENTATION NOTES
                # ============================================================
                # LW emission happens at the blanket outer surface, not at the
                # concrete. T_outer is solved from the quasi-steady blanket
                # energy balance (see derivation above) extended to include LW.
                #
                # The full balance has a T_outer^4 term (nonlinear). The
                # linearized Option A (T_ref = (T_conc+T_sky)/2) was tested
                # and gave >0.5°F error for cold-sky/high-solar combinations
                # (e.g., T_sky=-10°C, T_conc=55°C, G=900 W/m²) — exceeding
                # the 0.5°F threshold — so Option B is used: the linearized
                # estimate is refined by 2 Newton steps on the full T^4 balance.
                # 2 Newton steps from the linearized start converge to <0.01°C
                # of the fully converged solution (see test_lw_linearization_error_bounded).
                #
                # Linearized initial guess (T_ref = midpoint(T_conc, T_sky)):
                # h_rad = 4·ε·σ·T_ref^3
                # T_outer_0 = [α·G + T_conc/R + h_conv·T_amb + h_rad·T_sky]
                #             / (h_conv + h_rad + 1/R)
                #
                # Newton refinement (2 steps; F = residual of full T^4 balance):
                # F(T_o) = h_conv·(T_o−T_amb) + ε·σ·(T_o_K^4−T_sky_K^4) − α·G − (T_conc−T_o)/R
                # dF/dT_o = h_conv + 4·ε·σ·T_o_K^3 + 1/R
                # T_o ← T_o − F/dF
                #
                # LW flux into concrete is attenuated by the same factor f as
                # solar (derivation above), since both fluxes travel the same
                # path from blanket outer surface through R_blanket to concrete:
                # q_lw_effective = ε·σ·(T_outer^4 − T_sky^4) · f
                # where f = h_top_combined / h_forced_convection.
                #
                # DIAGNOSTIC SEMANTICS (consistent with PR 2 solar):
                # q_LW_incident_history: raw LW exchange at blanket outer surface
                # q_LW_history:          effective LW flux into concrete (= incident · f)
                # T_outer_C_history:     blanket outer-surface temperature
                # All three populated when diagnostic_outputs=True.
                # ─────────────────────────────────────────────────────────────────────
                _solar_W_m2 = getattr(environment, 'solar_W_m2', None)
                _env_hours  = getattr(environment, 'hours', None)
                if _solar_W_m2 is not None and len(_solar_W_m2) > 0:
                    _G_t = float(np.interp(t_hrs, _env_hours, _solar_W_m2))
                else:
                    _G_t = 0.0
                # Incident flux at blanket outer surface (negative = heat into system)
                _q_solar_scalar   = -_alpha_sol_top * _G_t
                _q_solar_incident = np.where(grid.is_air[_j_top, :], 0.0, _q_solar_scalar)
                # Attenuated flux reaching concrete surface: multiply by f = h_top/h_conv
                _q_solar_eff      = _q_solar_incident * (_h_top_combined / _h_conv_top)
                # PR 3: LW balance at blanket outer surface — Newton iteration (Option B)
                # Starts from linearized estimate, then refines 2 steps on full T^4 balance.
                # Option A (pure linearization) was tested and gave >0.5°F error for
                # cold-sky / high-solar combinations → Option B required per spec.
                _T_sky_arr = getattr(environment, 'T_sky_C', None)
                if _T_sky_arr is not None and len(_T_sky_arr) > 0 and _env_hours is not None:
                    _T_sky_C_t = float(np.interp(t_hrs, _env_hours, _T_sky_arr))
                    _T_sky_K   = _T_sky_C_t + 273.15
                    _T_conc_C  = _T_surf_c   # (nx,) at j=iy_concrete_start
                    # Linearized initial guess: T_ref = midpoint(T_conc, T_sky)
                    _T_ref_K   = 0.5 * (_T_conc_C + _T_sky_C_t) + 273.15
                    _h_rad0    = 4.0 * _emis_top * STEFAN_BOLTZMANN * _T_ref_K ** 3
                    _denom_o   = _h_conv_top + _h_rad0 + 1.0 / _R_blanket_SI
                    _num_o     = (_alpha_sol_top * _G_t
                                  + _T_conc_C / _R_blanket_SI
                                  + _h_conv_top * _T_amb_C
                                  + _h_rad0 * _T_sky_C_t)
                    _T_outer_C = _num_o / _denom_o   # initial (linearized) estimate (nx,) °C
                    # 2 Newton steps on full T^4 balance (vectorised over nx columns)
                    for _lw_iter in range(2):
                        _T_o_K = _T_outer_C + 273.15
                        _F_lw  = (_h_conv_top * (_T_outer_C - _T_amb_C)
                                  + _emis_top * STEFAN_BOLTZMANN * (_T_o_K ** 4 - _T_sky_K ** 4)
                                  - _alpha_sol_top * _G_t
                                  - (_T_conc_C - _T_outer_C) / _R_blanket_SI)
                        _dF_lw = (_h_conv_top
                                  + 4.0 * _emis_top * STEFAN_BOLTZMANN * _T_o_K ** 3
                                  + 1.0 / _R_blanket_SI)
                        _T_outer_C = _T_outer_C - _F_lw / _dF_lw
                    _T_outer_K = _T_outer_C + 273.15
                    _q_lw_incident = _emis_top * STEFAN_BOLTZMANN * (_T_outer_K ** 4 - _T_sky_K ** 4)
                    _q_lw_incident = np.where(grid.is_air[_j_top, :], 0.0, _q_lw_incident)
                    _q_lw_eff  = _q_lw_incident * (_h_top_combined / _h_conv_top)
                else:
                    _T_outer_C     = _T_surf_c.copy()    # no sky data: T_outer ≈ T_concrete
                    _T_sky_C_t     = float(_T_surf_c.mean())  # fallback: LW term zeros when T_conc ≈ T_sky
                    _T_sky_K       = _T_sky_C_t + 273.15
                    _q_lw_incident = np.zeros(nx)
                    _q_lw_eff      = np.zeros(nx)
                _q_top_c   = _q_conv_c + _q_evap_c + _q_solar_eff + _q_lw_eff
                _q_cond_c  = k_y_face[_j_top, :] * (T[_j_top + 1, :] - _T_surf_c) / _dy_conc_top
                T_new[_j_top, :] = np.where(
                    grid.is_air[_j_top, :],
                    T_initial_C[_j_top, :],
                    _T_surf_c + dt_step * (_q_cond_c - _q_top_c) / (rho_cp[_j_top, :] * _dy_half_c),
                )
                _q_top = _q_top_c   # for diagnostic output
                # blanket row j=0 is inactive — pin to initial value
                T_new[0, :] = T_initial_C[0, :]
            else:
                # Normal: h_conv only at blanket surface (j=0); blanket node is active
                _dy_top  = float(grid.y[1] - grid.y[0])
                _dy_half = _dy_top / 2.0
                _T_surf  = T[0, :]
                _q_conv_top = _h_conv_top * (_T_surf - _T_amb_C)
                _q_evap  = np.zeros(nx)
                if t_hrs < _cure_time:
                    for _i in range(nx):
                        if not _top_row_is_air[_i]:
                            _q_evap[_i] = menzel_evaporation(
                                t_hrs, float(_T_surf[_i]), _T_amb_C,
                                _RH_frac, _wind_eff
                            )
                _q_top = _q_conv_top + _q_evap
                _q_cond_top = k_y_face[0, :] * (T[1, :] - _T_surf) / _dy_top
                T_new[0, :] = np.where(
                    _top_row_is_air,
                    T_initial_C[0, :],
                    _T_surf + dt_step * (_q_cond_top - _q_top) / (rho_cp[0, :] * _dy_half),
                )

            # (b) Side column (i = ix_cs, j = 1..iy_concrete_end) — form face BC
            #     PR 4: adds side-face solar (F_vert projection + daytime gate) to
            #     the Sprint-0 form-R convective path. Direct LW deferred to Sprint 2:
            #     removing R_FORM_EFFECTIVE_SI caused a 10°F gradient regression from
            #     nighttime overcooling (unattenuated h_forced + LW). The solar term
            #     is the physically dominant contribution to the corner RMS error.
            #
            #     The interior stencil ran with full-dx cell mass for column ix_cs.
            #     For the actual half-cell (dx/2 wide), lateral conductance per unit
            #     volume doubles. Post-correct: add missing extra lateral term and
            #     subtract the side-face flux.
            _ics = _ix_cs
            # Side BC applies to the FORM FACE: from iy_concrete_start to
            # iy_concrete_end-1. Bottom row (iy_concrete_end) is dominated by
            # the soil BC; blanket row (j=0) at i=ix_cs is not a form face.
            _js_side  = slice(1, grid.iy_concrete_end)        # j=1..iy_concrete_end-1
            _T_side_c = T[_js_side, _ics]                     # concrete surface temps at form face

            # Solar: flat F_vert projection gated to daytime (6 AM – 6 PM).
            # Default F_vert = 0.0 (construction.vertical_solar_factor defaults to 0.0).
            # α_form·F_vert·G term is zero in PR 6; PR 7 activates it.
            _daytime = is_daytime(t_hrs, _placement_hour)

            # PR 6: quasi-steady T_outer_form solve on steel-form outer surface.
            # Mirrors the top-BC T_outer pattern (PR 3, lines ~1641-1691).
            # h_conv_vert: pure forced convection (no R_form baked in).
            # R_FORM_EFFECTIVE_SI: form + contact-film resistance (calibrated 0.0862 m²·K/W).
            # F_SKY_VERT=0.5: view factor vertical face → sky.
            # T_ground = T_amb (Sprint 3 Barber model will upgrade).
            _T_sky_arr = getattr(environment, "T_sky_C", None)
            if _T_sky_arr is not None and len(_T_sky_arr) > 0 and _env_hours is not None:
                _T_sky_C_t = float(np.interp(t_hrs, _env_hours, _T_sky_arr))
                _T_sky_K   = _T_sky_C_t + 273.15

                # Linearized initial guess using a sky/ground blended reference temp
                _T_eff_sky_C = F_SKY_VERT * _T_sky_C_t + (1.0 - F_SKY_VERT) * _T_amb_C
                _T_ref_K     = 0.5 * (_T_side_c + _T_eff_sky_C) + 273.15
                _h_rad0      = 4.0 * _emis_side * STEFAN_BOLTZMANN * _T_ref_K ** 3
                _denom_f     = _h_conv_vert + _h_rad0 + 1.0 / R_FORM_EFFECTIVE_SI
                _num_f       = (_T_side_c / R_FORM_EFFECTIVE_SI
                                + _h_conv_vert * _T_amb_C
                                + _h_rad0 * _T_eff_sky_C)
                _T_outer_form_C = _num_f / _denom_f

                # 2 Newton steps on full T⁴ residual (same convergence guarantee as PR 3)
                for _f_iter in range(2):
                    _T_o_K   = _T_outer_form_C + 273.15
                    _T_gnd_K = _T_amb_C + 273.15   # T_ground = T_amb (PR 6 baseline)
                    _F_lw = (
                        _h_conv_vert * (_T_outer_form_C - _T_amb_C)
                        + _emis_side * STEFAN_BOLTZMANN * F_SKY_VERT
                          * (_T_o_K ** 4 - _T_sky_K ** 4)
                        + _emis_side * EMIS_GROUND * STEFAN_BOLTZMANN * (1.0 - F_SKY_VERT)
                          * (_T_o_K ** 4 - _T_gnd_K ** 4)
                        - (_T_side_c - _T_outer_form_C) / R_FORM_EFFECTIVE_SI
                        # α_form · F_vert · G · daytime — zero in PR 6 (F_vert=0.0); PR 7 adds.
                    )
                    _dF_lw = (
                        _h_conv_vert
                        + 4.0 * _emis_side * STEFAN_BOLTZMANN * F_SKY_VERT
                          * _T_o_K ** 3
                        + 4.0 * _emis_side * EMIS_GROUND * STEFAN_BOLTZMANN
                          * (1.0 - F_SKY_VERT) * _T_o_K ** 3
                        + 1.0 / R_FORM_EFFECTIVE_SI
                    )
                    _T_outer_form_C = _T_outer_form_C - _F_lw / _dF_lw

                # Heat flux positive = OUT of concrete (T_conc > T_outer means heat leaves)
                _q_side_total = -(_T_outer_form_C - _T_side_c) / R_FORM_EFFECTIVE_SI
            else:
                # No sky data: T_outer ≈ T_conc, LW path off
                _T_outer_form_C = _T_side_c.copy()
                _q_side_total   = _h_conv_vert * (_T_side_c - _T_amb_C)

            _kxp_side  = k_x_face[_js_side, _ics]
            _lat_extra = _kxp_side * (T[_js_side, _ics + 1] - T[_js_side, _ics]) * inv_dxsq
            T_new[_js_side, _ics] += (
                dt_step * (_lat_extra - 2.0 * _q_side_total / dx) / rho_cp[_js_side, _ics]
            )

            # (c) Centerline (i = nx-1): zero-flux symmetry
            #     The stencil does NOT update i=nx-1 (edge column). We compute it
            #     from scratch as a half-cell with mirror neighbour at i=nx-2.
            #     All rows j=1..ny-2 (interior) computed vectorially.
            _kxl_cl  = k_x_face[1:-1, -1]   # face between i=nx-2 and i=nx-1
            _lat_cl  = 2.0 * _kxl_cl * (T[1:-1, -2] - T[1:-1, -1]) * inv_dxsq
            _vert_cl = y_coef[:, 0] * (
                k_y_face[1:, -1] * (T[2:, -1] - T[1:-1, -1]) / dy_plus[:, 0]
                - k_y_face[:-1, -1] * (T[1:-1, -1] - T[:-2, -1]) / dy_minus[:, 0]
            )
            T_new[1:-1, -1] = (
                T[1:-1, -1]
                + dt_step * (_lat_cl + _vert_cl + Q[1:-1, -1]) / rho_cp[1:-1, -1]
            )
            # j=0 centerline: top BC already set T_new[0, -1]; add lateral correction.
            _kxl_cl0 = k_x_face[0, -1]
            T_new[0, -1] += dt_step * (
                2.0 * _kxl_cl0 * (T[0, -2] - T[0, -1]) * inv_dxsq
            ) / rho_cp[0, -1]

            # (b2) Top-form corner: quarter-cell energy balance
            #
            # The corner cell (j=iy_concrete_start, i=ix_concrete_start) lies on
            # BOTH the top surface (half-cell-y) and the form face (half-cell-x)
            # simultaneously — it is a QUARTER cell.  Stages (a) and (b) handle
            # each half-cell independently and their combined result uses
            # inconsistent thermal mass assumptions.  This overwrite replaces
            # both with a single consistent quarter-cell energy balance where
            # each face area is correctly halved.
            _jc = grid.iy_concrete_start    # j = 1
            _ic = grid.ix_concrete_start    # i = ix_cs
            _dy_q = (grid.y[_jc + 1] - grid.y[_jc]) / 2.0   # dy_quarter = dy_concrete_top / 2
            _dx_q = dx / 2.0                                   # dx_quarter

            _k_right = k_x_face[_jc, _ic]     # face between ix_cs and ix_cs+1
            _k_below = k_y_face[_jc, _ic]     # face between j=1 and j=2
            _dy_full = grid.y[_jc + 1] - grid.y[_jc]

            _T_c   = T[_jc, _ic]
            _q_r   = _k_right * (T[_jc, _ic + 1] - _T_c) / dx          # W/m² in from right
            _q_b   = _k_below * (T[_jc + 1, _ic] - _T_c) / _dy_full    # W/m² in from below

            # Top face: full radiation balance through blanket (same physics as main top BC)
            _f_atten         = _h_top_combined / _h_conv_top
            _q_top_conv      = _h_top_combined * (_T_c - _T_amb_C)
            _q_top_solar     = -_alpha_sol_top * _G_t * _f_atten
            _T_outer_corner_K = float(_T_outer_C[_ic]) + 273.15
            _q_top_lw        = _emis_top * STEFAN_BOLTZMANN * (_T_outer_corner_K ** 4 - _T_sky_K ** 4) * _f_atten
            _q_top_c         = _q_top_conv + _q_top_solar + _q_top_lw
            if t_hrs < _cure_time:
                _q_top_c += menzel_evaporation(
                    t_hrs, float(_T_c), _T_amb_C, _RH_frac, _wind_eff
                )

            # Corner side face: reuse T_outer_form_C[0] from the main side-BC Newton solve
            # (row 0 of _js_side = j=iy_concrete_start = corner cell). Mirrors how the top
            # corner reuses _T_outer_C[_ic] at line 1806. PR 7 activates solar term.
            if _T_sky_arr is not None and len(_T_sky_arr) > 0 and _env_hours is not None:
                _T_outer_form_corner_C = float(_T_outer_form_C[0])
                _q_s = -(_T_outer_form_corner_C - _T_c) / R_FORM_EFFECTIVE_SI
            else:
                _q_s = _h_conv_vert * (_T_c - _T_amb_C)

            # Quarter-cell area (per unit depth) for mass and volumetric source
            _cell_A = _dx_q * _dy_q
            T_new[_jc, _ic] = _T_c + dt_step * (
                  _q_r   * _dy_q      # flux in × right-face area
                + _q_b   * _dx_q      # flux in × bottom-face area
                - _q_top_c * _dx_q    # flux out × top-face area
                - _q_s   * _dy_q      # flux out × side-face area
                + Q[_jc, _ic] * _cell_A
            ) / (rho_cp[_jc, _ic] * _cell_A)

            # (d) Far-left soil (i = 0, all soil rows): Dirichlet T_gw.
            #     Design doc D5: "assume T = T_gw at x=0 in soil (soil extends to
            #     deep ground temperature far from concrete)." Dirichlet rather than
            #     adiabatic prevents the cold soil-extension from acting as a lateral
            #     heat sink for the form-edge concrete bottom — matching CW's domain
            #     which has no soil to the left of the form face.
            T_new[grid.iy_concrete_end + 1 :, 0] = _T_gw_C

            # (e) Bottom row: Dirichlet T_gw
            T_new[-1, :] = _T_gw_C

            # (f) Air cells: hold at initial value (belt-and-suspenders against drift)
            T_new[grid.is_air] = T_initial_C[grid.is_air]

        else:
            # LEGACY MODE — do not use for new validation.  Preserves pre-PR-2
            # convection (+5.5 LW-equivalent baked in) for regression / ablation runs.
            # Uses _h_top_legacy (= _h_convective_legacy + blanket R in series).
            t_hrs = t / 3600.0
            _T_amb_F = ambient_temp_F(t_hrs, environment, _placement_hour)
            _T_amb_C = (_T_amb_F - 32.0) * 5.0 / 9.0

            _dy_top = float(grid.y[1] - grid.y[0])
            _T_surf = T[0, :]
            _q_conv = _h_top_legacy * (_T_surf - _T_amb_C)
            _q_evap = np.zeros(nx)
            if t_hrs < _cure_time:
                for _i in range(nx):
                    if not _top_row_is_soil[_i]:
                        _q_evap[_i] = menzel_evaporation(
                            t_hrs, float(_T_surf[_i]), _T_amb_C,
                            _RH_frac, _wind_eff
                        )
            _q_top = _q_conv + _q_evap
            _q_cond_in = k_y_face[0, :] * (T[1, :] - _T_surf) / _dy_top
            _dy_half = _dy_top / 2.0
            T_new[0, :] = _T_surf + dt_step * (_q_cond_in - _q_top) / (rho_cp[0, :] * _dy_half)

            # Bottom BC
            if boundary_mode == "v2_equivalent":
                T_new[-1, :] = _T_gw_C
            else:
                T_new[-1, :] = T_new[-2, :]

            # Sides: adiabatic for both legacy top_bc modes
            T_new[:, 0]  = T_new[:, 1]
            T_new[:, -1] = T_new[:, -2]

        # Pin air cells to initial value for ALL modes (air has rho_cp=1 placeholder;
        # leaving them at the BC-updated value would corrupt next-step flux calcs).
        if grid.is_air.any():
            T_new[grid.is_air] = T_initial_C[grid.is_air]
        # When blanket is pure-R (full_2d or skip_blanket_node), pin it to initial value.
        if _use_pure_r_blanket:
            T_new[grid.is_blanket] = T_initial_C[grid.is_blanket]

        T, T_new = T_new, T
        t += dt_step
        n_steps += 1

        if abs(t - t_samples[next_idx]) < 1e-6:
            T_out[next_idx]     = T
            alpha_out[next_idx] = alpha
            te_out[next_idx]    = te
            t_out[next_idx]     = t
            if _use_top_bc:
                top_flux_out[next_idx] = _q_top
                T_amb_out[next_idx]    = _T_amb_C
                if q_solar_out is not None:
                    # Recompute at exact sample time (BC was computed at t-dt).
                    _jt  = grid.iy_concrete_start
                    _sol = getattr(environment, 'solar_W_m2', None)
                    _hrs = getattr(environment, 'hours', None)
                    _t_s_hr = t / 3600.0
                    if _sol is not None and len(_sol) > 0:
                        _G_s      = float(np.interp(_t_s_hr, _hrs, _sol))
                        _inc_s    = np.where(grid.is_air[_jt, :], 0.0, -_alpha_sol_top * _G_s)
                        _atten_s  = _inc_s * (_h_top_combined / _h_conv_top)
                    else:
                        _G_s     = 0.0
                        _inc_s   = np.zeros(nx)
                        _atten_s = np.zeros(nx)
                    q_solar_out[next_idx]          = _atten_s   # flux entering concrete
                    q_solar_incident_out[next_idx] = _inc_s     # raw at blanket outer surface

                    # PR 4: top-face conv and evap (always available, no sky data needed)
                    _T_amb_s_F = ambient_temp_F(_t_s_hr, environment, _placement_hour)
                    _T_amb_s_C = (_T_amb_s_F - 32.0) * 5.0 / 9.0
                    _T_c_s     = T[_jt, :]               # post-swap, shape (nx,)
                    _q_conv_s  = _h_top_combined * (_T_c_s - _T_amb_s_C)
                    _q_evap_s  = np.zeros(nx)
                    if _t_s_hr < _cure_time:
                        for _i_s in range(grid.ix_concrete_start, nx):
                            _q_evap_s[_i_s] = menzel_evaporation(
                                _t_s_hr, float(_T_c_s[_i_s]), _T_amb_s_C, _RH_frac, _wind_eff
                            )
                    q_conv_out[next_idx] = _q_conv_s
                    q_evap_out[next_idx] = _q_evap_s

                    # PR 3: recompute LW at exact sample time
                    _T_sky_arr_s = getattr(environment, 'T_sky_C', None)
                    if _T_sky_arr_s is not None and len(_T_sky_arr_s) > 0 and _hrs is not None:
                        _T_sky_s   = float(np.interp(_t_s_hr, _hrs, _T_sky_arr_s))
                        _T_sky_s_K = _T_sky_s + 273.15
                        _T_ref_s_K = 0.5 * (_T_c_s + _T_sky_s) + 273.15
                        _h_rad_s0  = 4.0 * _emis_top * STEFAN_BOLTZMANN * _T_ref_s_K ** 3
                        _denom_s   = _h_conv_top + _h_rad_s0 + 1.0 / _R_blanket_SI
                        _num_s     = (_alpha_sol_top * _G_s
                                      + _T_c_s / _R_blanket_SI
                                      + _h_conv_top * _T_amb_s_C
                                      + _h_rad_s0 * _T_sky_s)
                        _T_outer_s = _num_s / _denom_s   # linearized initial estimate
                        for _lw_iter_s in range(2):
                            _T_o_s_K = _T_outer_s + 273.15
                            _F_s  = (_h_conv_top * (_T_outer_s - _T_amb_s_C)
                                     + _emis_top * STEFAN_BOLTZMANN * (_T_o_s_K ** 4 - _T_sky_s_K ** 4)
                                     - _alpha_sol_top * _G_s
                                     - (_T_c_s - _T_outer_s) / _R_blanket_SI)
                            _dF_s = (_h_conv_top
                                     + 4.0 * _emis_top * STEFAN_BOLTZMANN * _T_o_s_K ** 3
                                     + 1.0 / _R_blanket_SI)
                            _T_outer_s = _T_outer_s - _F_s / _dF_s
                        _T_outer_s_K = _T_outer_s + 273.15
                        _q_lw_inc_s  = _emis_top * STEFAN_BOLTZMANN * (_T_outer_s_K ** 4 - _T_sky_s_K ** 4)
                        _q_lw_inc_s  = np.where(grid.is_air[_jt, :], 0.0, _q_lw_inc_s)
                        q_lw_out[next_idx]          = _q_lw_inc_s * (_h_top_combined / _h_conv_top)
                        q_lw_incident_out[next_idx] = _q_lw_inc_s
                        T_outer_out[next_idx]       = np.where(grid.is_air[_jt, :], np.nan, _T_outer_s)
                    else:
                        q_lw_out[next_idx]          = np.zeros(nx)
                        q_lw_incident_out[next_idx] = np.zeros(nx)
                        T_outer_out[next_idx]       = np.zeros(nx)
                    # PR 4: top-face total (sum of all 4 components)
                    q_top_total_out[next_idx] = (
                        q_conv_out[next_idx] + q_evap_out[next_idx]
                        + q_solar_out[next_idx] + q_lw_out[next_idx]
                    )

                    # PR 6: side-face Newton re-solve at sample time (mirrors top-BC
                    # sample-time recompute at lines ~1929-1961 for T_outer_C_history).
                    _T_side_c_s  = T[1:grid.iy_concrete_end, _ix_cs]
                    _daytime_s   = is_daytime(_t_s_hr, _placement_hour)
                    _T_sky_arr_s = getattr(environment, "T_sky_C", None)
                    if _T_sky_arr_s is not None and len(_T_sky_arr_s) > 0 and _hrs is not None:
                        _T_sky_C_s     = float(np.interp(_t_s_hr, _hrs, _T_sky_arr_s))
                        _T_sky_K_s     = _T_sky_C_s + 273.15
                        _T_eff_sky_C_s = F_SKY_VERT * _T_sky_C_s + (1.0 - F_SKY_VERT) * _T_amb_s_C
                        _T_ref_K_s     = 0.5 * (_T_side_c_s + _T_eff_sky_C_s) + 273.15
                        _h_rad0_s      = 4.0 * _emis_side * STEFAN_BOLTZMANN * _T_ref_K_s ** 3
                        _denom_s       = _h_conv_vert + _h_rad0_s + 1.0 / R_FORM_EFFECTIVE_SI
                        _num_s         = (_T_side_c_s / R_FORM_EFFECTIVE_SI
                                          + _h_conv_vert * _T_amb_s_C
                                          + _h_rad0_s * _T_eff_sky_C_s)
                        _T_outer_form_s = _num_s / _denom_s
                        _T_gnd_K_s = _T_amb_s_C + 273.15
                        for _fs_iter in range(2):
                            _T_o_K_s = _T_outer_form_s + 273.15
                            _F_s = (
                                _h_conv_vert * (_T_outer_form_s - _T_amb_s_C)
                                + _emis_side * STEFAN_BOLTZMANN * F_SKY_VERT
                                  * (_T_o_K_s ** 4 - _T_sky_K_s ** 4)
                                + _emis_side * EMIS_GROUND * STEFAN_BOLTZMANN * (1.0 - F_SKY_VERT)
                                  * (_T_o_K_s ** 4 - _T_gnd_K_s ** 4)
                                - (_T_side_c_s - _T_outer_form_s) / R_FORM_EFFECTIVE_SI
                            )
                            _dF_s = (
                                _h_conv_vert
                                + 4.0 * _emis_side * STEFAN_BOLTZMANN * F_SKY_VERT
                                  * _T_o_K_s ** 3
                                + 4.0 * _emis_side * EMIS_GROUND * STEFAN_BOLTZMANN
                                  * (1.0 - F_SKY_VERT) * _T_o_K_s ** 3
                                + 1.0 / R_FORM_EFFECTIVE_SI
                            )
                            _T_outer_form_s = _T_outer_form_s - _F_s / _dF_s
                        _T_o_K_s = _T_outer_form_s + 273.15
                        q_side_lw_out[next_idx] = (
                            _emis_side * STEFAN_BOLTZMANN * F_SKY_VERT
                            * (_T_o_K_s ** 4 - _T_sky_K_s ** 4)
                            + _emis_side * EMIS_GROUND * STEFAN_BOLTZMANN * (1.0 - F_SKY_VERT)
                            * (_T_o_K_s ** 4 - _T_gnd_K_s ** 4)
                        )
                        T_outer_form_out[next_idx] = _T_outer_form_s
                        q_side_conv_out[next_idx]  = _h_conv_vert * (_T_outer_form_s - _T_amb_s_C)
                    else:
                        _T_outer_form_s = _T_side_c_s.copy()
                        q_side_lw_out[next_idx]    = np.zeros(_n_side_rows)
                        T_outer_form_out[next_idx] = _T_outer_form_s
                        q_side_conv_out[next_idx]  = _h_conv_vert * (_T_side_c_s - _T_amb_s_C)
                    q_side_solar_out[next_idx] = np.full(
                        _n_side_rows, -_alpha_sol_side * _F_vert * _G_s * _daytime_s
                    )
                    q_side_total_out[next_idx] = (
                        q_side_conv_out[next_idx]
                        + q_side_lw_out[next_idx]
                        + q_side_solar_out[next_idx]
                    )
            next_idx += 1

    # ------------------------------------------------------------------ #
    # J. Assemble result                                                   #
    # ------------------------------------------------------------------ #
    peak_flat = int(np.argmax(T_out))
    peak_loc_3d = np.unravel_index(peak_flat, T_out.shape)   # (t_idx, j, i)
    peak_T_loc: tuple[int, int] = (int(peak_loc_3d[1]), int(peak_loc_3d[2]))

    return HydrationResult(
        t_s=t_out,
        T_field_C=T_out,
        alpha_field=alpha_out,
        t_e_field_hrs=te_out,
        dt_inner_s=dt_inner,
        n_inner_steps=n_steps,
        n_output_samples=n_samples,
        peak_T_C=float(T_out.max()),
        peak_T_location=peak_T_loc,
        peak_alpha=float(alpha_out.max()),
        centerline_T_C=T_out[:, :, -1],
        top_flux_W_m2_history=top_flux_out,
        T_amb_C_history=T_amb_out,
        q_solar_history=q_solar_out,
        q_solar_incident_history=q_solar_incident_out,
        q_LW_history=q_lw_out,
        q_LW_incident_history=q_lw_incident_out,
        T_outer_C_history=T_outer_out,
        q_conv_history=q_conv_out,
        q_evap_history=q_evap_out,
        q_top_total_history=q_top_total_out,
        q_side_solar_history=q_side_solar_out,
        q_side_LW_history=q_side_lw_out,
        q_side_conv_history=q_side_conv_out,
        q_side_total_history=q_side_total_out,
        T_outer_form_C_history=T_outer_form_out,
    )
