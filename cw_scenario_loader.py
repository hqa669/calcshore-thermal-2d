"""
cw_scenario_loader.py
=====================
Load a ConcreteWorks scenario from its native file formats and return a
fully-specified scenario object ready for the CalcShore thermal engine.

Three inputs:
  1. TEST.dat           — CW input file (positional, one value per line)
  2. TX__Austin.dat     — CW weather file (17-column hourly annual data)
  3. temp.txt           — CW 2D temperature export (for validation)

Usage
-----
    from cw_scenario_loader import load_cw_scenario

    scn = load_cw_scenario(
        input_dat="TEST.dat",
        weather_dat="TX__Austin.dat",
        cw_output_txt="temp.txt",              # optional, for validation
    )
    # scn.mix, scn.geometry, scn.construction, scn.environment  (CalcShore dataclasses)
    # scn.cw_validation      → Max/Min/Gradient/Ambient time-series

Author: CalcShore / Qinang Li
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from datetime import datetime
from typing import Dict, List, Optional, Tuple

import numpy as np


# ============================================================
# Hard-coded line indices for CW v2.1.3 .dat file
# Derived empirically from TEST.dat + UI screenshot cross-check.
# If CW changes file format, only this block needs updating.
# ============================================================
CW_DAT_INDEX = {
    # Header/meta
    'length_units':         1,
    'mass_units':           2,
    'temp_units':           4,
    'analysis_duration':    8,
    'placement_time_str':  10,
    'project_location':    11,
    'placement_date':      12,
    # Geometry
    'member_width_ft':     21,
    'member_depth_ft':     22,
    'member_length_ft':    23,
    # Mix design
    'cement_lb_yd3':       54,
    'water_lb_yd3':        55,
    'coarse_agg_lb_yd3':   56,
    'fine_agg_lb_yd3':     57,
    'air_content_pct':     58,
    'silica_fume_lb_yd3':  59,
    'class_c_fly_ash_lb_yd3': 60,
    'class_f_fly_ash_lb_yd3': 61,
    'fly_ash_CaO_pct':     62,
    'ggbfs_lb_yd3':        63,
    # Cement & hydration (Schindler-Folliard regression output)
    'cement_type':        361,
    'activation_energy_J_mol': 385,
    'tau_hrs':            386,
    'beta':               387,
    'alpha_u':            388,
    'Hu_J_kg':            389,
    # Aggregate types
    'coarse_agg_type':    390,
    'fine_agg_type':      393,
    # Derived thermal
    'CTE_microstrain_F':                403,
    'thermal_conductivity_BTU_hr_ft_F': 404,
    'aggregate_Cp_BTU_lb_F':            405,
    # Maturity (informational)
    'maturity_a':         407,
    'maturity_b':         408,
    # Forms / curing
    'form_color':          431,
    'form_type':           432,
    'form_removal_hrs':    433,
    'blanket_R_value':     434,
    'placement_temp_F':    438,
    'soil_temp_F':         439,
    'delay_strip_to_cure_hrs': 445,
    'footing_subbase':     446,
    'steel_cover_in':      447,
    'top_cure_blanket_time_hrs': 478,
    # Service life (informational)
    'Dref_x1e13_m2_s':     491,
    'chloride_aging_m':    492,
    'steel_type':          495,
    'exposure_class':      508,
}


# ============================================================
# Scenario dataclasses (mirrors thermal_engine_v2.py structure)
# ============================================================

@dataclass
class CWMixDesign:
    """Mix design and hydration parameters as CW extracts them from a TCP.

    Customer-set composition fields: cement_type_I_II_lb_yd3,
    fly_ash_F_lb_yd3, fly_ash_C_lb_yd3, ggbfs_lb_yd3, silica_fume_lb_yd3,
    water_lb_yd3, coarse_agg_lb_yd3, fine_agg_lb_yd3, air_content_pct,
    fly_ash_CaO_pct.

    Hydration kinetics (regressed by CW from the TCP, or directly entered):
    activation_energy_J_mol (Ea), tau_hrs (τ), beta (β), alpha_u (αu),
    Hu_J_kg. These drive the Arrhenius hydration model; see
    docs/engine_v3_release_notes.md for the validated kinetics range.

    Derived thermal properties: thermal_conductivity_BTU_hr_ft_F,
    aggregate_Cp_BTU_lb_F, CTE_microstrain_F. Default values reflect the
    39.1% SCM Reference mix (limestone coarse aggregate, siliceous fine sand).
    """
    cement_type_I_II_lb_yd3: float = 350.0
    fly_ash_F_lb_yd3: float = 0.0
    fly_ash_C_lb_yd3: float = 0.0
    ggbfs_lb_yd3: float = 0.0
    silica_fume_lb_yd3: float = 0.0
    water_lb_yd3: float = 253.0
    coarse_agg_lb_yd3: float = 1800.0
    fine_agg_lb_yd3: float = 1100.0
    air_content_pct: float = 5.0
    fly_ash_CaO_pct: float = 0.0

    cement_type: str = "Type I/II"
    coarse_agg_type: str = "Limestone"
    fine_agg_type: str = "Siliceous River Sand"

    # Hydration regression output (directly from CW)
    activation_energy_J_mol: float = 26457.9
    tau_hrs: float = 29.401
    beta: float = 0.895
    alpha_u: float = 0.75852
    Hu_J_kg: float = 424143.0

    # Derived thermal properties
    thermal_conductivity_BTU_hr_ft_F: float = 1.56
    aggregate_Cp_BTU_lb_F: float = 0.20
    CTE_microstrain_F: float = 4.25

    @property
    def total_cementitious_lb_yd3(self) -> float:
        return (self.cement_type_I_II_lb_yd3 + self.fly_ash_F_lb_yd3
                + self.fly_ash_C_lb_yd3 + self.ggbfs_lb_yd3 + self.silica_fume_lb_yd3)

    @property
    def wcm(self) -> float:
        return self.water_lb_yd3 / self.total_cementitious_lb_yd3

    @property
    def concrete_density_lb_ft3(self) -> float:
        """Computed from mix total, corrected for air content (CW approach)."""
        total_solids = (self.cement_type_I_II_lb_yd3 + self.fly_ash_F_lb_yd3
                        + self.fly_ash_C_lb_yd3 + self.ggbfs_lb_yd3
                        + self.silica_fume_lb_yd3 + self.water_lb_yd3
                        + self.coarse_agg_lb_yd3 + self.fine_agg_lb_yd3)
        return total_solids / 27.0 * (1.0 - self.air_content_pct / 100.0)


@dataclass
class CWGeometry:
    """Placement geometry for a half-mat footing.

    Customer-set dimensions: depth_ft (mat thickness), width_ft, length_ft.
    The engine v3 validated envelope is the 40×60×8 ft half-mat footing
    (width=40, length=60, depth=8); geometry outside this envelope is
    out-of-envelope (R6; routes to the geometry coverage series).

    analysis_type and shape default to "2-D" and "Rect. footing" — these
    match the engine v3 solver's 2D half-mat formulation and should not be
    changed unless the solver is extended.
    """
    depth_ft: float = 8.0
    width_ft: float = 40.0
    length_ft: float = 60.0
    analysis_type: str = "2-D"      # CW allows 1-D / 2-D / 3-D
    shape: str = "Rect. footing"


@dataclass
class CWConstruction:
    """Construction and curing parameters for the placement.

    Customer-set fields: placement_temp_F, placement_date, placement_hour,
    blanket_R_value, form_removal_hrs, side_cure_method, top_cure_method_1,
    top_cure_method_2, top_cure_blanket_time_hrs, delay_strip_to_cure_hrs,
    soil_temp_F, footing_subbase, form_color, form_orientation.

    Radiative tunables (customer-adjustable): solar_absorptivity_top,
    solar_absorptivity_side, emissivity_top, emissivity_side,
    vertical_solar_factor (None = use form_orientation lookup; override for
    sweeps or ablation only).

    form_type accepts {"steel", "plywood", "plastic_liner"} (case-insensitive;
    the loader normalizes to lowercase at parse time). Engine v3 is validated
    only for steel forms (ADR-04, reinforced PR 20). Plywood is
    out-of-envelope (R_form = 0.17 m²·K/W from ACI 306R-88, not
    CW-validated); using it emits a UserWarning at runtime via
    resolve_r_form(). plastic_liner is not implemented and raises
    NotImplementedError. See docs/engine_v3_release_notes.md for the
    full envelope.

    soil_lag_hrs and soil_damping are deprecated no-ops (defaulted to no-op
    values in PR 17, ADR-08); retained for future cold-climate work where
    ground-air differentials become identifiable.
    """
    placement_temp_F: float = 60.0
    placement_hour: int = 5
    placement_date: str = "2026/7/15"

    blanket_R_value: float = 5.67
    solar_absorptivity_top: float = 0.65   # ASHRAE default for light-gray concrete/steel
    emissivity_top: float = 0.88           # blanket outer surface emissivity
    solar_absorptivity_side: float = 0.65  # steel form (similar to light concrete)
    emissivity_side: float = 0.88          # steel form emissivity
    # F_vert is normally resolved via form_orientation → F_VERT_BY_ORIENTATION
    # lookup in the engine (PR 8, M6d). This field is an optional override used
    # for sweeps, calibration, and ablation tests. None means use the lookup.
    vertical_solar_factor: float | None = None
    form_removal_hrs: float = 168.0
    form_color: str = "Red"
    form_type: str = "steel"
    form_orientation: str = "unknown"     # {"south","east","west","north","unknown"}

    side_cure_method: str = "Wet Blanket"
    top_cure_method_1: str = "Wet blanket"
    top_cure_method_2: str = "Clear Plastic"
    top_cure_blanket_time_hrs: float = 2.0
    delay_strip_to_cure_hrs: float = 1.0

    soil_temp_F: float = 80.0
    footing_subbase: str = "Limestone"

    # Barber soil model fields (Sprint 3 / PR 10–11, §2.3.1 coding_passdown_v4.md).
    # PR 17 (R4 disposition): defaults reverted to no-op pair — T_ground(t) = T_amb(t).
    # Empirical justification: §7.5.2 (14-mix library has zero climate variation; all TX
    # Austin; soil params identical across all mixes); Sprint 3 sweep showed damping
    # unidentifiable on MIX-01 (<0.01°F authority). Without identifiability data across
    # climates, defaults that introduce phase offset degrade fit monotonically.
    # Fields retained (not removed) for future cold-climate exports where diurnal
    # ground-air differentials become large enough to be identifiable. See ADR-08.
    soil_lag_hrs: float = 0.0     # PR 17: 5.0 → 0.0 (no-op)
    soil_damping: float = 1.0     # PR 17: 0.7 → 1.0 (no-op); deprecated, retained for cold-climate data


@dataclass
class CWEnvironment:
    """Climate time series and location metadata for the placement.

    Hourly boundary condition arrays (length = 24 × analysis_duration_days,
    starting at placement_date @ placement_hour): T_air_F, RH_pct,
    solar_W_m2, wind_m_s, cloud_cover, pressure_hPa. These are the true
    physics BC values the engine uses for each time step.

    Derived hourly arrays computed at load time from weather columns:
    T_dp_C (dew point, Magnus-Tetens), T_sky_C (sky temperature,
    Berdahl-Martin from T_air + cloud_cover).

    CW UI averages (may differ slightly from raw hourly max due to CW
    internal averaging): cw_ave_max_daily_temp_F, cw_ave_min_daily_temp_F,
    cw_ave_max_solar_W_m2, cw_ave_max_wind_m_s, cw_ave_max_RH_pct,
    cw_ave_min_RH_pct.

    Location metadata: lat_deg, lon_deg, elevation_m, location. Engine v3
    validated envelope is Austin TX summer (location="TX, Austin"). Other
    climates are out-of-envelope; see docs/engine_v3_release_notes.md.
    """
    # Hourly time series, length = 24 * analysis_duration_days
    hours: np.ndarray = field(default_factory=lambda: np.array([]))
    T_air_F: np.ndarray = field(default_factory=lambda: np.array([]))
    RH_pct: np.ndarray = field(default_factory=lambda: np.array([]))
    solar_W_m2: np.ndarray = field(default_factory=lambda: np.array([]))
    wind_m_s: np.ndarray = field(default_factory=lambda: np.array([]))
    cloud_cover: np.ndarray = field(default_factory=lambda: np.array([]))
    pressure_hPa: np.ndarray = field(default_factory=lambda: np.array([]))
    # Computed hourly fields (derived from weather columns at load time)
    T_dp_C: np.ndarray = field(default_factory=lambda: np.array([]))
    T_sky_C: np.ndarray = field(default_factory=lambda: np.array([]))

    # Per-day summary (daily max/min of the dry-bulb air temp)
    daily_max_F: List[float] = field(default_factory=list)
    daily_min_F: List[float] = field(default_factory=list)

    # CW UI-displayed averages (kept separately — may differ from raw hourly max)
    cw_ave_max_daily_temp_F: float = 92.3
    cw_ave_min_daily_temp_F: float = 74.6
    cw_ave_max_solar_W_m2: float = 848.1
    cw_ave_max_wind_m_s: float = 10.5
    cw_ave_max_RH_pct: float = 87.4
    cw_ave_min_RH_pct: float = 44.6

    lat_deg: float = 30.28
    lon_deg: float = -97.68
    elevation_m: float = 90.0
    location: str = "TX, Austin"


@dataclass
class CWValidationSeries:
    """CW 2D output parsed into canonical series."""
    time_hrs: np.ndarray = field(default_factory=lambda: np.array([]))
    T_max_xs_F: np.ndarray = field(default_factory=lambda: np.array([]))
    T_min_xs_F: np.ndarray = field(default_factory=lambda: np.array([]))
    T_diff_xs_F: np.ndarray = field(default_factory=lambda: np.array([]))
    T_ambient_F: np.ndarray = field(default_factory=lambda: np.array([]))

    # Full 2D field (n_time, nD, nW) in °F and the grid coordinates
    T_field_F: Optional[np.ndarray] = None
    depths_m: Optional[np.ndarray] = None
    widths_m: Optional[np.ndarray] = None

    # Centerline-specific convenience series (width index 0)
    T_core_center_F: Optional[np.ndarray] = None
    T_top_center_F: Optional[np.ndarray] = None
    T_bot_center_F: Optional[np.ndarray] = None


@dataclass
class CWScenario:
    """Complete scenario bundle."""
    mix: CWMixDesign
    geometry: CWGeometry
    construction: CWConstruction
    environment: CWEnvironment
    cw_validation: Optional[CWValidationSeries] = None
    raw_dat_values: Dict = field(default_factory=dict)


# ============================================================
# .dat parser
# ============================================================

def parse_cw_dat(path: str) -> Tuple[CWMixDesign, CWGeometry, CWConstruction, Dict]:
    """Parse a CW v2.1.3 .dat input file.

    Returns (mix, geometry, construction, raw_values_dict).
    """
    with open(path, 'r', encoding='latin-1') as f:
        lines = [ln.rstrip('\r\n').strip() for ln in f.readlines()]

    if not lines or 'Concrete Works' not in lines[0]:
        raise ValueError(f"Not a CW .dat file: {path} (first line: {lines[0] if lines else '<empty>'})")

    def _get(idx: int) -> str:
        if idx >= len(lines):
            return ''
        return lines[idx].strip()

    def _getf(idx: int, default: float = 0.0) -> float:
        try:
            return float(_get(idx))
        except (ValueError, IndexError):
            return default

    # Raw passthrough for debugging / downstream use
    raw = {}
    for key, idx in CW_DAT_INDEX.items():
        raw[key] = _get(idx)

    # Parse placement hour from "5 am" / "10 pm" etc.
    placement_hour_str = _get(CW_DAT_INDEX['placement_time_str'])
    placement_hour = _parse_time_of_day(placement_hour_str)

    # Build MixDesign
    mix = CWMixDesign(
        cement_type_I_II_lb_yd3=_getf(CW_DAT_INDEX['cement_lb_yd3']),
        fly_ash_F_lb_yd3=_getf(CW_DAT_INDEX['class_f_fly_ash_lb_yd3']),
        fly_ash_C_lb_yd3=_getf(CW_DAT_INDEX['class_c_fly_ash_lb_yd3']),
        ggbfs_lb_yd3=_getf(CW_DAT_INDEX['ggbfs_lb_yd3']),
        silica_fume_lb_yd3=_getf(CW_DAT_INDEX['silica_fume_lb_yd3']),
        water_lb_yd3=_getf(CW_DAT_INDEX['water_lb_yd3']),
        coarse_agg_lb_yd3=_getf(CW_DAT_INDEX['coarse_agg_lb_yd3']),
        fine_agg_lb_yd3=_getf(CW_DAT_INDEX['fine_agg_lb_yd3']),
        air_content_pct=_getf(CW_DAT_INDEX['air_content_pct']),
        fly_ash_CaO_pct=_getf(CW_DAT_INDEX['fly_ash_CaO_pct']),
        cement_type=_get(CW_DAT_INDEX['cement_type']),
        coarse_agg_type=_get(CW_DAT_INDEX['coarse_agg_type']),
        fine_agg_type=_get(CW_DAT_INDEX['fine_agg_type']),
        activation_energy_J_mol=_getf(CW_DAT_INDEX['activation_energy_J_mol']),
        tau_hrs=_getf(CW_DAT_INDEX['tau_hrs']),
        beta=_getf(CW_DAT_INDEX['beta']),
        alpha_u=_getf(CW_DAT_INDEX['alpha_u']),
        Hu_J_kg=_getf(CW_DAT_INDEX['Hu_J_kg']),
        thermal_conductivity_BTU_hr_ft_F=_getf(CW_DAT_INDEX['thermal_conductivity_BTU_hr_ft_F']),
        aggregate_Cp_BTU_lb_F=_getf(CW_DAT_INDEX['aggregate_Cp_BTU_lb_F']),
        CTE_microstrain_F=_getf(CW_DAT_INDEX['CTE_microstrain_F']),
    )

    # Build Geometry
    geom = CWGeometry(
        depth_ft=_getf(CW_DAT_INDEX['member_depth_ft']),
        width_ft=_getf(CW_DAT_INDEX['member_width_ft']),
        length_ft=_getf(CW_DAT_INDEX['member_length_ft']),
    )

    # Build Construction
    constr = CWConstruction(
        placement_temp_F=_getf(CW_DAT_INDEX['placement_temp_F']),
        placement_hour=placement_hour,
        placement_date=_get(CW_DAT_INDEX['placement_date']),
        blanket_R_value=_getf(CW_DAT_INDEX['blanket_R_value']),
        form_removal_hrs=_getf(CW_DAT_INDEX['form_removal_hrs']),
        form_color=_get(CW_DAT_INDEX['form_color']),
        form_type=_get(CW_DAT_INDEX['form_type']).strip().lower(),  # trust boundary: normalize to lowercase
        top_cure_blanket_time_hrs=_getf(CW_DAT_INDEX['top_cure_blanket_time_hrs']),
        delay_strip_to_cure_hrs=_getf(CW_DAT_INDEX['delay_strip_to_cure_hrs']),
        soil_temp_F=_getf(CW_DAT_INDEX['soil_temp_F']),
        footing_subbase=_get(CW_DAT_INDEX['footing_subbase']),
    )

    return mix, geom, constr, raw


def _parse_time_of_day(s: str) -> int:
    """Parse '5 am' / '12 pm' / '10 PM' → 24-hour int."""
    s = s.strip().lower()
    try:
        parts = s.split()
        h = int(parts[0])
        suffix = parts[1] if len(parts) > 1 else ''
        if suffix == 'pm' and h != 12:
            h += 12
        elif suffix == 'am' and h == 12:
            h = 0
        return h
    except (ValueError, IndexError):
        return 5  # CW's UI default


# ============================================================
# Weather (TX__Austin.dat) parser
# ============================================================

# Column index map for CW weather file (17 columns).
# Decoded by matching values against UI-displayed averages.
WEATHER_COL = {
    'month':           0,
    'day':             1,
    'hour':            2,
    'solar_global':    3,   # W/m² total horizontal
    'solar_diffuse':   4,   # W/m²
    'T_dry_bulb_C':    5,   # °C
    'T_dew_C':         6,   # °C
    'RH_pct':          7,
    'T_aux_C':         8,
    'wind_dir_deg':    9,   # 999 = missing
    'wind_speed_m_s': 10,
    'wind_speed_2':   11,
    'solar_direct':   12,
    'cloud_cover':    13,   # octas
    'sky_temp_C':     14,
    'pressure_hPa':   15,
    'wind_gust':      16,
}


def dew_point_C(T_air_C: np.ndarray, RH_pct: np.ndarray) -> np.ndarray:
    """Magnus-Tetens dew-point approximation.

    Parameters
    ----------
    T_air_C : np.ndarray
        Air temperature in °C.
    RH_pct : np.ndarray
        Relative humidity in percent (0-100).

    Returns
    -------
    np.ndarray
        Dew-point temperature in °C.

    Notes
    -----
    Valid for -45°C < T < +60°C, 1% < RH < 100%. Accuracy ±0.4°C.
    Reference: Alduchov & Eskridge (1996), "Improved Magnus form
    approximation of saturation vapor pressure", J. Appl. Meteor. 35(4).

    Examples
    --------
    >>> dew_point_C(np.array([25.0]), np.array([50.0]))
    array([13.86...])
    """
    a = 17.625
    b = 243.04
    RH_frac = np.clip(np.asarray(RH_pct, dtype=float) / 100.0, 0.01, 1.0)
    T = np.asarray(T_air_C, dtype=float)
    gamma = np.log(RH_frac) + (a * T) / (b + T)
    return (b * gamma) / (a - gamma)


def sky_temp_C(T_air_C: np.ndarray, T_dp_C: np.ndarray,
               cloud_cover_octas: np.ndarray) -> np.ndarray:
    """Berdahl-Martin effective sky temperature with cloud correction.

    Parameters
    ----------
    T_air_C : np.ndarray
        Air temperature in °C.
    T_dp_C : np.ndarray
        Dew-point temperature in °C (e.g., from dew_point_C).
    cloud_cover_octas : np.ndarray
        Cloud cover in octas (0 = clear, 8 = overcast).

    Returns
    -------
    np.ndarray
        Effective sky temperature in °C, typically 5-30°C below T_air
        for clear skies, approaching T_air under overcast.

    Notes
    -----
    Source: Berdahl & Martin (1984), "Emissivity of clear skies",
    Solar Energy 32(5). Cloud correction: eps_sky increases quadratically
    with cloud fraction toward unity at overcast.

    Examples
    --------
    >>> sky_temp_C(np.array([25.0]), np.array([15.0]), np.array([0.0]))
    array([...])  # several degrees below 25°C
    """
    T_dp = np.asarray(T_dp_C, dtype=float)
    T_air = np.asarray(T_air_C, dtype=float)
    N = np.clip(np.asarray(cloud_cover_octas, dtype=float) / 8.0, 0.0, 1.0)
    eps_clear = 0.711 + 0.56 * (T_dp / 100.0) + 0.73 * (T_dp / 100.0) ** 2
    eps_sky = eps_clear + (1.0 - eps_clear) * N ** 2
    T_air_K = T_air + 273.15
    T_sky_K = T_air_K * eps_sky ** 0.25
    return T_sky_K - 273.15


def parse_cw_weather(
    path: str,
    placement_date: str,
    placement_hour: int,
    duration_days: int = 7,
) -> CWEnvironment:
    """Parse a CW weather .dat file and extract hourly arrays for the analysis period.

    Args
    ----
    path : CW weather file path
    placement_date : 'YYYY/M/D' or 'YYYY/MM/DD'
    placement_hour : 0–23
    duration_days : analysis length

    Returns a CWEnvironment with hours[], T_air_F[], RH_pct[], solar_W_m2[], etc.
    Hour 0 in the output = placement moment. Arrays have length 24*duration_days.
    """
    with open(path, 'r', encoding='latin-1') as f:
        lines = [ln.strip() for ln in f.readlines()]

    if not lines or 'weather' not in lines[0].lower():
        raise ValueError(f"Not a CW weather file: {path}")

    location = lines[1].strip()
    lat, lon, elev = [float(x) for x in lines[2].split()]

    data_rows = []
    for ln in lines[3:]:
        parts = ln.split()
        if len(parts) == 17:
            data_rows.append([float(x) for x in parts])
    data = np.array(data_rows)

    if data.size == 0:
        raise ValueError(f"No data rows found in {path}")

    # Parse placement date
    try:
        dt = datetime.strptime(placement_date, '%Y/%m/%d')
    except ValueError:
        # Try '%Y/%#m/%#d' style for single-digit months/days
        dt = datetime.strptime(placement_date.replace('/', '-'), '%Y-%m-%d')

    target_month = dt.month
    target_day = dt.day

    # Find the starting row (matching month, day, hour)
    # Weather file hour is 1–24; placement_hour is 0–23 → map 0→1, 23→24
    target_hour_1based = placement_hour + 1
    mask = ((data[:, WEATHER_COL['month']] == target_month)
            & (data[:, WEATHER_COL['day']] == target_day)
            & (data[:, WEATHER_COL['hour']] == target_hour_1based))
    idx_start = np.where(mask)[0]
    if len(idx_start) == 0:
        raise ValueError(
            f"Could not find placement moment (month={target_month}, day={target_day}, "
            f"hour={target_hour_1based}) in weather file."
        )
    i0 = int(idx_start[0])
    n_hours = 24 * duration_days
    i1 = min(i0 + n_hours, len(data))
    slice_data = data[i0:i1]
    if len(slice_data) < n_hours:
        # Wrap to beginning of year (weather file is annual periodic)
        extra_needed = n_hours - len(slice_data)
        slice_data = np.vstack([slice_data, data[:extra_needed]])

    # Build hourly arrays
    hours = np.arange(n_hours, dtype=float)  # 0, 1, 2, ..., hours since placement
    T_C = slice_data[:, WEATHER_COL['T_dry_bulb_C']]
    T_air_F = T_C * 9.0 / 5.0 + 32.0
    RH = slice_data[:, WEATHER_COL['RH_pct']]
    solar = slice_data[:, WEATHER_COL['solar_global']]
    wind = slice_data[:, WEATHER_COL['wind_speed_m_s']]
    cloud = slice_data[:, WEATHER_COL['cloud_cover']]
    pressure = slice_data[:, WEATHER_COL['pressure_hPa']]

    # Computed hourly fields
    T_dp = dew_point_C(T_C, RH)
    T_sky = sky_temp_C(T_C, T_dp, cloud)

    # Daily max/min of dry-bulb
    n_full_days = n_hours // 24
    T_air_F_days = T_air_F[:n_full_days * 24].reshape(n_full_days, 24)
    daily_max_F = list(T_air_F_days.max(axis=1))
    daily_min_F = list(T_air_F_days.min(axis=1))

    env = CWEnvironment(
        hours=hours,
        T_air_F=T_air_F,
        RH_pct=RH,
        solar_W_m2=solar,
        wind_m_s=wind,
        cloud_cover=cloud,
        pressure_hPa=pressure,
        T_dp_C=T_dp,
        T_sky_C=T_sky,
        daily_max_F=daily_max_F,
        daily_min_F=daily_min_F,
        lat_deg=lat,
        lon_deg=lon,
        elevation_m=elev,
        location=location,
    )
    return env


# ============================================================
# CW temp.txt parser (2D output)
# ============================================================

def parse_cw_temp_output(path: str) -> CWValidationSeries:
    """Parse the `temp.txt` CW export.

    File format:
        Line 0: description ("Distance is from the upper corner...")
        Line 1: header ("Time  w1/d1  w2/d1  ...  Gradient  Ambient")
        Lines 2+: one row per timestep, N spatial cols + Gradient + Ambient

    Returns a CWValidationSeries with per-timestep Max/Min/Gradient/Ambient plus
    the full 2D field.
    """
    with open(path, 'r', encoding='latin-1') as f:
        lines = [ln.strip() for ln in f.readlines() if ln.strip()]

    header = lines[1].split()
    # coords are header[1:-2] (skip 'Time' and trailing 'Gradient', 'Ambient')
    coords = header[1:-2]
    coord_pairs = [tuple(map(float, c.split('/'))) for c in coords]
    widths_m_unique = sorted({c[0] for c in coord_pairs}, reverse=True)  # 6.1 → 0
    depths_m_unique = sorted({c[1] for c in coord_pairs})                # 0 → 2.44
    nW = len(widths_m_unique)
    nD = len(depths_m_unique)
    assert nW * nD == len(coords), f"Grid mismatch: {nW}×{nD} ≠ {len(coords)}"

    # Parse data
    times = []
    fields = []
    gradients = []
    ambients = []
    for ln in lines[2:]:
        toks = ln.split()
        if len(toks) < 3:
            continue
        t = float(toks[0])
        n_spatial = len(toks) - 3
        vals = np.array([float(x) for x in toks[1:1 + n_spatial]])
        if vals.size != nW * nD:
            continue  # malformed line
        # Reshape: file stores all widths at depth 0, then all widths at depth 1, etc.
        field_2d = vals.reshape(nD, nW)
        times.append(t)
        fields.append(field_2d)
        gradients.append(float(toks[-2]))
        ambients.append(float(toks[-1]))

    times = np.array(times)
    field_C = np.array(fields)         # (n_time, nD, nW) in °C
    field_F = field_C * 9.0 / 5.0 + 32.0

    def C_to_F(c):
        return c * 9.0 / 5.0 + 32.0

    T_amb_F = C_to_F(np.array(ambients))
    grad_F = np.array(gradients) * 9.0 / 5.0

    # Per-timestep full-section max/min
    flat = field_F.reshape(field_F.shape[0], -1)
    T_max_xs_F = flat.max(axis=1)
    T_min_xs_F = flat.min(axis=1)
    T_diff_xs_F = T_max_xs_F - T_min_xs_F

    # Centerline series.
    # Width indexing: file stores widths DECREASING (6.1 → 0).
    # "Distance from the upper CORNER in width direction" means w=0 is AT the corner,
    # w=6.1 is furthest from corner = centerline (for 40 ft = 12.2 m half-symmetric mat).
    # Therefore centerline is width index 0 (the FIRST column).
    T_core_center_F = field_F[:, nD // 2, 0]
    T_top_center_F = field_F[:, 0, 0]
    T_bot_center_F = field_F[:, -1, 0]

    return CWValidationSeries(
        time_hrs=times,
        T_max_xs_F=T_max_xs_F,
        T_min_xs_F=T_min_xs_F,
        T_diff_xs_F=T_diff_xs_F,
        T_ambient_F=T_amb_F,
        T_field_F=field_F,
        depths_m=np.array(depths_m_unique),
        widths_m=np.array(widths_m_unique),
        T_core_center_F=T_core_center_F,
        T_top_center_F=T_top_center_F,
        T_bot_center_F=T_bot_center_F,
    )


# ============================================================
# Top-level loader
# ============================================================

def load_cw_scenario(
    input_dat: str,
    weather_dat: Optional[str] = None,
    cw_output_txt: Optional[str] = None,
    cw_ui_overrides: Optional[Dict[str, float]] = None,
) -> CWScenario:
    """Load a full CW scenario from its three native files.

    Args
    ----
    input_dat : path to CW input file (required)
    weather_dat : path to CW weather .dat (optional but recommended)
    cw_output_txt : path to CW temp.txt export (optional, for validation)
    cw_ui_overrides : dict with UI-displayed environment averages that cannot
        be inferred from the weather file alone (e.g. {'ave_max_wind_m_s': 10.5})

    Returns
    -------
    CWScenario with mix, geometry, construction, environment [, cw_validation]
    """
    mix, geom, constr, raw = parse_cw_dat(input_dat)

    if weather_dat is not None:
        env = parse_cw_weather(
            weather_dat,
            placement_date=constr.placement_date,
            placement_hour=constr.placement_hour,
            duration_days=7,
        )
    else:
        env = CWEnvironment()

    # Apply UI overrides (for values CW displays but can't be derived)
    if cw_ui_overrides:
        for k, v in cw_ui_overrides.items():
            if hasattr(env, k):
                setattr(env, k, v)

    validation = None
    if cw_output_txt is not None and os.path.exists(cw_output_txt):
        validation = parse_cw_temp_output(cw_output_txt)

    return CWScenario(
        mix=mix,
        geometry=geom,
        construction=constr,
        environment=env,
        cw_validation=validation,
        raw_dat_values=raw,
    )


# ============================================================
# CLI / Diagnostic
# ============================================================

def summarize_scenario(scn: CWScenario) -> str:
    """Human-readable summary of a loaded scenario."""
    m, g, c, e = scn.mix, scn.geometry, scn.construction, scn.environment
    lines = [
        "=" * 72,
        " CW SCENARIO SUMMARY",
        "=" * 72,
        "",
        "┌─ Project",
        f"│  Location            = {e.location}",
        f"│  Placement           = {c.placement_date} @ {c.placement_hour:02d}:00",
        f"│  Lat/Lon/Elev        = {e.lat_deg:.3f}° / {e.lon_deg:.3f}° / {e.elevation_m:.0f} m",
        "",
        "┌─ Geometry",
        f"│  W × D × L           = {g.width_ft} × {g.depth_ft} × {g.length_ft} ft",
        f"│  Type                = {g.analysis_type}",
        "",
        "┌─ Mix design (lb/yd³)",
        f"│  Cement ({m.cement_type}) = {m.cement_type_I_II_lb_yd3:.0f}",
        f"│  Class F fly ash     = {m.fly_ash_F_lb_yd3:.0f}  ({m.fly_ash_CaO_pct:.0f}% CaO)",
        f"│  GGBFS               = {m.ggbfs_lb_yd3:.0f}",
        f"│  Silica fume         = {m.silica_fume_lb_yd3:.0f}",
        f"│  Water               = {m.water_lb_yd3:.0f}  (w/cm = {m.wcm:.3f})",
        f"│  Coarse agg ({m.coarse_agg_type})  = {m.coarse_agg_lb_yd3:.0f}",
        f"│  Fine agg ({m.fine_agg_type})  = {m.fine_agg_lb_yd3:.0f}",
        f"│  Air content         = {m.air_content_pct:.1f}%",
        f"│  Unit weight (calc)  = {m.concrete_density_lb_ft3:.1f} lb/ft³",
        "",
        "┌─ Hydration (CW Schindler regression)",
        f"│  Ea                  = {m.activation_energy_J_mol:.1f} J/mol",
        f"│  τ                   = {m.tau_hrs:.3f} hr",
        f"│  β                   = {m.beta:.4f}",
        f"│  α_u                 = {m.alpha_u:.4f}",
        f"│  Hu                  = {m.Hu_J_kg:.0f} J/kg",
        "",
        "┌─ Thermal properties",
        f"│  k                   = {m.thermal_conductivity_BTU_hr_ft_F:.3f} BTU/hr·ft·°F"
        f"  = {m.thermal_conductivity_BTU_hr_ft_F * 1.7307:.3f} W/m·K",
        f"│  Cp (agg)            = {m.aggregate_Cp_BTU_lb_F:.4f} BTU/lb·°F"
        f"  = {m.aggregate_Cp_BTU_lb_F * 4186.8:.0f} J/kg·K",
        f"│  CTE                 = {m.CTE_microstrain_F:.2f} µε/°F",
        "",
        "┌─ Construction",
        f"│  Placement temp      = {c.placement_temp_F}°F",
        f"│  Blanket R-value     = {c.blanket_R_value}",
        f"│  Form type / color   = {c.form_type} / {c.form_color}",
        f"│  Form removal        = {c.form_removal_hrs} hrs",
        f"│  Top cure at         = {c.top_cure_blanket_time_hrs} hrs",
        f"│  Side cure method    = {c.side_cure_method}",
        f"│  Subbase             = {c.footing_subbase}",
        f"│  Soil temp           = {c.soil_temp_F}°F",
        "",
        "┌─ Environment",
    ]
    if e.T_air_F.size > 0:
        lines += [
            f"│  Hourly arrays       = {len(e.T_air_F)} hrs loaded",
            f"│  T_air range (raw)   = {e.T_air_F.min():.1f} → {e.T_air_F.max():.1f} °F",
            f"│  Daily max (mean)    = {np.mean(e.daily_max_F):.1f} °F"
            f"  (CW UI: {e.cw_ave_max_daily_temp_F})",
            f"│  Daily min (mean)    = {np.mean(e.daily_min_F):.1f} °F"
            f"  (CW UI: {e.cw_ave_min_daily_temp_F})",
            f"│  Solar max           = {e.solar_W_m2.max():.1f} W/m²"
            f"  (CW UI: {e.cw_ave_max_solar_W_m2})",
            f"│  Wind max (raw)      = {e.wind_m_s.max():.1f} m/s"
            f"  (CW UI: {e.cw_ave_max_wind_m_s})",
            f"│  RH range            = {e.RH_pct.min():.1f} → {e.RH_pct.max():.1f}%"
            f"  (CW UI min/max: {e.cw_ave_min_RH_pct}/{e.cw_ave_max_RH_pct})",
        ]
    else:
        lines.append("│  (no weather file loaded)")

    if scn.cw_validation is not None:
        v = scn.cw_validation
        pk_idx = int(np.argmax(v.T_max_xs_F))
        dT_idx = int(np.argmax(v.T_diff_xs_F))
        lines += [
            "",
            "┌─ CW validation data",
            f"│  Timesteps loaded    = {len(v.time_hrs)}",
            f"│  Peak Max T in xs    = {v.T_max_xs_F[pk_idx]:.1f}°F @ {v.time_hrs[pk_idx]:.1f} hr",
            f"│  Peak ΔT (gradient)  = {v.T_diff_xs_F[dT_idx]:.1f}°F @ {v.time_hrs[dT_idx]:.1f} hr",
            f"│  Peak centerline core= {v.T_core_center_F.max():.1f}°F",
        ]

    lines.append("=" * 72)
    return "\n".join(lines)


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print("Usage: python cw_scenario_loader.py <TEST.dat> [weather.dat] [temp.txt]")
        sys.exit(1)

    dat = sys.argv[1]
    weather = sys.argv[2] if len(sys.argv) > 2 else None
    temp = sys.argv[3] if len(sys.argv) > 3 else None

    # Example UI overrides (what CW shows but weather file can't give us exactly)
    ui_overrides = {
        'cw_ave_max_wind_m_s': 10.5,    # CW UI value differs from raw weather
    }

    scn = load_cw_scenario(dat, weather, temp, cw_ui_overrides=ui_overrides)
    print(summarize_scenario(scn))
