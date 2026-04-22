"""
CalcShore Thermal Engine v2.0
================================
1D Finite-Difference Transient Heat Transfer with Cement Hydration

Improvements over v1:
  Step 1: NumPy vectorization (<1s runtime)
  Step 2: Variable thermal properties (CW Eq 23-25)
  Step 3: Soil nodes below concrete (CW Table 6, Eq 44)
  Step 4: Blanket thermal mass nodes (CW Eq 43)
  Step 5: Menzel evaporation model (CW Eq 38-39)

Physics:
  ρ·Cp·(∂T/∂t) = k·(∂²T/∂x²) + Q(t)
  Hydration: α(te) = αu·exp(-(τ/te)^β)  (Schindler-Folliard)
  Equivalent age: Arrhenius temperature dependence

Author: CalcShore / Qinang Li
"""

import numpy as np
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from typing import Dict, Any, Optional


# ============================================================
# Constants
# ============================================================
R_GAS = 8.314       # J/(mol·K)
T_REF_K = 296.15    # 23°C — Arrhenius reference temperature


# ============================================================
# Data Classes
# ============================================================

@dataclass
class MixDesign:
    """Concrete mix design — all quantities in imperial (lb/yd³)."""
    cement_type_I_II_lb_yd3: float = 350.0
    fly_ash_F_lb_yd3: float = 125.0
    ggbfs_lb_yd3: float = 100.0
    silica_fume_lb_yd3: float = 0.0
    water_lb_yd3: float = 253.0
    coarse_agg_lb_yd3: float = 1800.0
    fine_agg_lb_yd3: float = 1100.0
    air_content_pct: float = 5.0

    # Hydration parameters
    activation_energy_J_mol: float = 26457.9
    tau_hrs: float = 29.401
    beta: float = 0.895
    alpha_u: float = 0.75852
    Hu_J_kg: float = 424143.0

    # Thermal properties
    thermal_conductivity_BTU_hr_ft_F: float = 1.56
    aggregate_Cp_BTU_lb_F: float = 0.20
    concrete_density_lb_ft3: float = 150.0

    @property
    def total_cementitious_lb_yd3(self):
        return (self.cement_type_I_II_lb_yd3 + self.fly_ash_F_lb_yd3 +
                self.ggbfs_lb_yd3 + self.silica_fume_lb_yd3)

    @property
    def wcm(self):
        return self.water_lb_yd3 / self.total_cementitious_lb_yd3


@dataclass
class Geometry:
    """Element geometry in feet."""
    depth_ft: float = 8.0
    width_ft: float = 40.0
    length_ft: float = 60.0


@dataclass
class Construction:
    """Construction and curing parameters."""
    placement_temp_F: float = 60.0
    soil_temp_F: float = 80.0
    blanket_R_value: float = 5.67    # Physical blanket R (hr·ft²·°F/BTU) — mapped to effective R by engine
    form_removal_hrs: float = 168.0
    avg_annual_temp_F: float = 68.9  # Austin TX ≈ 20.5°C


@dataclass
class Environment:
    """Environmental conditions."""
    avg_max_temp_F: float = 92.3
    avg_min_temp_F: float = 74.6
    avg_max_solar_W_m2: float = 848.1
    avg_max_wind_m_s: float = 10.5
    placement_hour: int = 5  # 5 AM
    relative_humidity: float = 0.50  # 50% RH for Austin TX July
    daily_max_F: list = field(default_factory=lambda: [
        91.2, 92.3, 92.7, 92.5, 91.8, 92.5, 92.5, 93.2
    ])
    daily_min_F: list = field(default_factory=lambda: [
        74.7, 74.3, 74.1, 74.8, 75.0, 74.7, 74.1, 74.8
    ])


# ============================================================
# Calibration Parameters
# ============================================================
CALIBRATION = {
    'Hu_factor': 1.06,       # Compensates for 1D vs 2D + variable property effects
    # Calibrated across 15 CW 2D mixes: 14/15 pass ±5°F
    # Peak: mean=+0.1°F, std=1.6°F, max|err|=3.4°F
    # ΔT:   mean=+2.6°F, std=1.5°F, max|err|=4.4°F (excluding MIX-15 cold placement outlier)
}

# v2 domain defaults
V2_DEFAULTS = {
    'n_blanket': 1,              # Single blanket node (thermal mass)
    'blanket_R_effective': 0.32,  # Effective blanket R-value (hr·ft²·°F/BTU)
                                  # Calibrated against CW 2D at R=5.67 input
    'soil_type': 'sand',          # Best match to CW's subbase model
    'soil_depth_m': 3.0,          # 3m soil domain below concrete
    'n_soil': 20,                 # Soil discretization
}

# v1 calibration for backward compatibility / step-by-step validation
CALIBRATION_V1 = {
    'Hu_factor': 0.90,
    'R_top_effective': 0.15,
    'h_bot': 5.0,
    'q_evap': 18.0,
}


# ============================================================
# Soil Properties (CW Manual Table 6)
# ============================================================
SOIL_PROPERTIES = {
    "clay":      {"k": 1.3,  "rho": 1460, "Cp": 880},
    "granite":   {"k": 2.79, "rho": 2630, "Cp": 775},
    "limestone": {"k": 2.15, "rho": 2320, "Cp": 810},
    "marble":    {"k": 2.80, "rho": 2680, "Cp": 830},
    "quartzite": {"k": 5.38, "rho": 2640, "Cp": 1105},
    "sandstone": {"k": 2.90, "rho": 2150, "Cp": 745},
    "sand":      {"k": 0.27, "rho": 1515, "Cp": 800},
    "top_soil":  {"k": 0.52, "rho": 2050, "Cp": 1840},
}


# ============================================================
# Blanket Properties
# ============================================================
BLANKET_DEFAULTS = {
    "thickness_m": 0.02,
    "density": 320,      # kg/m³
    "Cp": 2000,          # J/(kg·K)
    # k derived from R-value: k = thickness / R_SI
}


# ============================================================
# Blanket R-value Mapping
# ============================================================
def blanket_r_effective(r_user):
    """Map user's physical blanket R-value to engine effective R.
    
    Calibrated against ConcreteWorks 2D rectangular footing.
    Log fit from CW R-sweep at R=0.18, 0.28, 0.50, 1.00, 5.67,
    anchored to 15-mix validated R_eff=0.32 at R_user=5.67.
    
    Args:
        r_user: Physical blanket R-value (hr·ft²·°F/BTU), range 0.18–10.0
    Returns:
        Effective R-value for 1D engine (hr·ft²·°F/BTU)
    """
    import numpy as np
    r_eff = 0.0491 * np.log(max(r_user, 0.18)) + 0.2391
    return max(r_eff, 0.15)  # floor to prevent CFL instability


# ============================================================
# Unit Conversions
# ============================================================
def F_to_K(T_F): return (T_F - 32) * 5/9 + 273.15
def K_to_F(T_K): return (T_K - 273.15) * 9/5 + 32
def F_to_C(T_F): return (T_F - 32) * 5/9
def C_to_F(T_C): return T_C * 9/5 + 32


# ============================================================
# Ambient Temperature
# ============================================================
def ambient_temp_F(t_hrs: float, env: Environment) -> float:
    """Sinusoidal diurnal model. Peak at 3 PM, min at ~5 AM."""
    h = (env.placement_hour + t_hrs) % 24.0
    d = min(int((env.placement_hour + t_hrs) / 24.0), len(env.daily_max_F) - 1)
    avg = (env.daily_max_F[d] + env.daily_min_F[d]) / 2
    amp = (env.daily_max_F[d] - env.daily_min_F[d]) / 2
    return avg - amp * np.cos(2 * np.pi * (h - 15.0) / 24.0)


# ============================================================
# Hydration Model (vectorized)
# ============================================================
def arrhenius_vec(T_K, Ea):
    """Arrhenius acceleration factor (vectorized)."""
    return np.exp(-Ea / R_GAS * (1.0/T_K - 1.0/T_REF_K))

def hydration_rate_vec(te, tau, beta, au):
    """dα/d(te) — vectorized hydration rate."""
    te_safe = np.maximum(te, 0.01)
    r = tau / te_safe
    return au * beta * r**beta / te_safe * np.exp(-r**beta)

def hydration_alpha_vec(te, tau, beta, au):
    """α(te) — vectorized degree of hydration."""
    te_safe = np.maximum(te, 0.01)
    return au * np.exp(-(tau / te_safe)**beta)

# Scalar versions for backward compat
def arrhenius(T_K, Ea):
    return np.exp(-Ea / R_GAS * (1.0/T_K - 1.0/T_REF_K))

def hydration_rate(te, tau, beta, au):
    if te < 0.01: return 0.0
    r = tau / te
    return au * beta * r**beta / te * np.exp(-r**beta)

def hydration_alpha(te, tau, beta, au):
    if te < 0.01: return 0.0
    return au * np.exp(-(tau / te)**beta)


# ============================================================
# Menzel Evaporation Model (CW Eq 38-39)
# ============================================================
def saturated_vapor_pressure_mmHg(T_C):
    """ASHRAE saturation vapor pressure, returned in mmHg."""
    T_K = T_C + 273.15
    if T_C >= 0:
        C8, C9, C10, C11, C12, C13 = -5.8002206e3, -5.516256, -4.8640239e-2, 4.1764768e-5, -1.4452093e-8, 6.5459673
        ln_Pws = C8/T_K + C9 + C10*T_K + C11*T_K**2 + C12*T_K**3 + C13*np.log(T_K)
    else:
        C1, C2, C3, C4, C5, C6, C7 = -5.6745359e3, -5.1523058e-1, -9.677843e-3, 6.2215701e-7, 2.0747825e-9, -9.484024e-13, 4.1635019
        ln_Pws = C1/T_K + C2 + C3*T_K + C4*T_K**2 + C5*T_K**3 + C6*T_K**4 + C7*np.log(T_K)
    Pws_kPa = np.exp(ln_Pws)
    return Pws_kPa * 7.50062  # kPa to mmHg

def menzel_evaporation(t_hrs, T_surface_C, T_air_C, RH, wind_m_s):
    """
    Menzel evaporative cooling rate (W/m²).
    CW Eq 38-39.
    """
    e0 = saturated_vapor_pressure_mmHg(T_surface_C)
    ea = saturated_vapor_pressure_mmHg(T_air_C)
    
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


# ============================================================
# Variable Thermal Properties (CW Eq 23-25)
# ============================================================
def thermal_conductivity_variable(k_uc, alpha_node):
    """CW Eq 23: k varies with hydration degree."""
    return k_uc * (1.33 - 0.33 * alpha_node)

def specific_heat_variable(Wc, Wa, Ww, Ca, alpha_node, T_C_node, rho):
    """
    CW Eq 24-25: Cp varies with hydration, mixture, temperature.
    Van Breugel model.
    
    Wc = cement weight kg/m³, Wa = aggregate weight kg/m³, Ww = water weight kg/m³
    Ca = aggregate specific heat J/(kg·K)
    T_C_node = concrete temperature °C per node
    """
    Cw = 4186.0  # water specific heat J/(kg·K)
    c_cef = 8.4 * T_C_node + 339.0  # fictitious specific heat of hydrating cement
    
    Cp_node = (1.0 / rho) * (Wc * alpha_node * c_cef + Wc * (1.0 - alpha_node) * Ca + Wa * Ca + Ww * Cw)
    return Cp_node


# ============================================================
# ACI 207.2R Heat of Hydration
# ============================================================
HEAT_FACTORS = {
    'cement_I_II': 13.0, 'cement_IL': 10.5, 'cement_V': 11.0,
    'fly_ash_F': 4.5, 'fly_ash_C': 9.0, 'ggbfs': 11.0, 'silica_fume': 10.0,
}

DOT_LIMITS = {
    'fly_ash_F_max_pct': 35, 'fly_ash_C_max_pct': 25,
    'ggbfs_max_pct': 70, 'total_SCM_max_pct': 70,
    'wcm_max': 0.45, 'peak_temp_max_F': 160, 'delta_T_max_F': 35,
}

def section_2_2(mix: MixDesign) -> Dict[str, Any]:
    """ACI 207.2R adiabatic temperature rise and SCM compliance."""
    components = {
        'cement_I_II': mix.cement_type_I_II_lb_yd3,
        'fly_ash_F': mix.fly_ash_F_lb_yd3,
        'ggbfs': mix.ggbfs_lb_yd3,
        'silica_fume': mix.silica_fume_lb_yd3,
    }
    contribs = {n: (a/100)*HEAT_FACTORS[n]
                for n, a in components.items() if n in HEAT_FACTORS}
    total = mix.total_cementitious_lb_yd3
    rise = sum(contribs.values())

    scm = {}
    if mix.fly_ash_F_lb_yd3 > 0:
        p = mix.fly_ash_F_lb_yd3 / total * 100
        scm['fly_ash_F_pct'] = (round(p, 1), p <= DOT_LIMITS['fly_ash_F_max_pct'])
    if mix.ggbfs_lb_yd3 > 0:
        p = mix.ggbfs_lb_yd3 / total * 100
        scm['ggbfs_pct'] = (round(p, 1), p <= DOT_LIMITS['ggbfs_max_pct'])
    ts = (mix.fly_ash_F_lb_yd3 + mix.ggbfs_lb_yd3 + mix.silica_fume_lb_yd3) / total * 100
    scm['total_SCM_pct'] = (round(ts, 1), ts <= DOT_LIMITS['total_SCM_max_pct'])
    scm['wcm'] = (round(mix.wcm, 3), mix.wcm <= DOT_LIMITS['wcm_max'])

    return {
        'contributions': contribs,
        'adiabatic_rise_F': round(rise, 1),
        'adiabatic_rise_C': round(rise / 1.8, 1),
        'total_cementitious_lb_yd3': total,
        'effective_rate_F_per_100lb': round(rise / (total / 100), 1),
        'scm_checks': scm,
    }


# ============================================================
# 1D Finite Difference Thermal Solver — v1 VECTORIZED (Step 1)
# ============================================================
def solve_thermal_1d_v1_vectorized(
    mix: MixDesign,
    geom: Geometry,
    constr: Construction,
    env: Environment,
    dt_hrs: float = 0.02,
    duration_hrs: float = 168.0,
    n_nodes: int = 61,
    calibration: Optional[Dict] = None,
) -> Dict[str, Any]:
    """
    v1 physics with NumPy vectorization. Produces identical results to v1.
    Used to validate that vectorization introduces no accuracy change.
    """
    cal = calibration or CALIBRATION_V1

    # --- Material properties (SI) ---
    depth_m = geom.depth_ft * 0.3048
    dx = depth_m / (n_nodes - 1)

    rho = mix.concrete_density_lb_ft3 * 16.0185
    Cp = mix.aggregate_Cp_BTU_lb_F * 4186.8
    k = mix.thermal_conductivity_BTU_hr_ft_F * 1.7307
    alpha_d = k / (rho * Cp)

    Cc = mix.total_cementitious_lb_yd3 * 0.593276
    Ea = mix.activation_energy_J_mol
    tau = mix.tau_hrs
    beta_h = mix.beta
    au = mix.alpha_u
    Hu = mix.Hu_J_kg * cal['Hu_factor']

    # --- Time step (CFL stability) ---
    dt_s = dt_hrs * 3600.0
    dt_max = 0.4 * dx**2 / alpha_d
    if dt_s > dt_max:
        dt_hrs = dt_max / 3600.0
        dt_s = dt_max
    n_steps = int(duration_hrs / dt_hrs) + 1

    # --- Boundary coefficients ---
    R_top_SI = cal['R_top_effective'] * 0.1761
    wind = env.avg_max_wind_m_s * 0.4
    h_outer = 5.6 + 3.5 * wind + 5.5
    h_top = 1.0 / (R_top_SI + 1.0 / h_outer)

    h_bot = cal['h_bot']
    q_evap = cal['q_evap']

    # --- Initialize ---
    T = np.full(n_nodes, F_to_K(constr.placement_temp_F))
    te = np.full(n_nodes, 0.01)
    T_soil_K = F_to_K(constr.soil_temp_F)
    mid = n_nodes // 2

    out_int = max(1, int(0.5 / dt_hrs))
    times, cores, tops, bots, ambs, dts = [], [], [], [], [], []

    # --- Time-stepping loop (VECTORIZED) ---
    for step in range(n_steps):
        t = step * dt_hrs
        T_amb_K = F_to_K(ambient_temp_F(t, env))

        # Hydration — vectorized
        af = arrhenius_vec(T, Ea)
        te += dt_hrs * af
        da_dt = hydration_rate_vec(te, tau, beta_h, au) * af
        Q = Hu * Cc * da_dt / 3600.0

        # FD stencil — vectorized interior
        Tn = T.copy()
        d2T = (T[2:] - 2*T[1:-1] + T[:-2]) / dx**2
        Tn[1:-1] = T[1:-1] + (k * d2T + Q[1:-1]) / (rho * Cp) * dt_s

        # Top BC
        flux_top = (k * (T[1] - T[0]) / dx
                    + h_top * (T_amb_K - T[0])
                    - q_evap
                    + Q[0] * dx / 2)
        Tn[0] = T[0] + flux_top / (rho * Cp * dx / 2) * dt_s

        # Bottom BC
        flux_bot = (k * (T[-2] - T[-1]) / dx
                    + h_bot * (T_soil_K - T[-1])
                    + Q[-1] * dx / 2)
        Tn[-1] = T[-1] + flux_bot / (rho * Cp * dx / 2) * dt_s

        T = Tn

        if step % out_int == 0:
            tc = K_to_F(T[mid])
            tt = K_to_F(T[0])
            tb = K_to_F(T[-1])
            ta = K_to_F(T_amb_K)
            dT_val = max(tc - min(tt, tb), 0)
            times.append(t)
            cores.append(tc)
            tops.append(tt)
            bots.append(tb)
            ambs.append(ta)
            dts.append(dT_val)

    times = np.array(times)
    cores = np.array(cores)
    tops = np.array(tops)
    bots = np.array(bots)
    ambs = np.array(ambs)
    dts = np.array(dts)

    pi = np.argmax(cores)
    di = np.argmax(dts)

    peak_ok = cores[pi] <= DOT_LIMITS['peak_temp_max_F']
    dt_ok = dts[di] <= DOT_LIMITS['delta_T_max_F']

    return {
        'time_hrs': times, 'T_core_F': cores, 'T_top_F': tops,
        'T_bot_F': bots, 'T_ambient_F': ambs, 'delta_T_F': dts,
        'peak_core_F': round(float(cores[pi]), 1),
        'peak_time_hrs': round(float(times[pi]), 1),
        'max_delta_T_F': round(float(dts[di]), 1),
        'max_delta_T_time_hrs': round(float(times[di]), 1),
        'peak_limit_F': DOT_LIMITS['peak_temp_max_F'],
        'delta_T_limit_F': DOT_LIMITS['delta_T_max_F'],
        'peak_pass': bool(peak_ok),
        'delta_T_pass': bool(dt_ok),
        'overall_pass': bool(peak_ok and dt_ok),
    }


# ============================================================
# 1D Finite Difference Thermal Solver — v2 FULL (Steps 1-5)
# ============================================================
def solve_thermal_1d(
    mix: MixDesign,
    geom: Geometry,
    constr: Construction,
    env: Environment,
    dt_hrs: float = 0.02,
    duration_hrs: float = 168.0,
    n_concrete: int = 61,
    n_soil: int = 20,
    n_blanket: int = 1,
    soil_depth_m: float = 3.0,
    soil_type: str = "sand",
    calibration: Optional[Dict] = None,
    use_variable_properties: bool = True,
    use_menzel_evap: bool = False,
    blanket_evap_W_m2: float = 0.0,
) -> Dict[str, Any]:
    """
    v2 solver with blanket nodes + concrete nodes + soil nodes.
    
    Domain: [blanket(top)] → [concrete] → [soil(bottom)]
    
    Steps implemented:
      1. Vectorized NumPy operations
      2. Variable k and Cp (CW Eq 23-25)
      3. Soil nodes with physical properties
      4. Blanket thermal mass nodes
      5. Menzel evaporation model
    """
    cal = calibration or CALIBRATION

    # ─── Material properties (SI) ───
    depth_m = geom.depth_ft * 0.3048
    
    rho_c = mix.concrete_density_lb_ft3 * 16.0185      # kg/m³
    Cp_c_base = mix.aggregate_Cp_BTU_lb_F * 4186.8     # J/(kg·K) base value
    k_uc = mix.thermal_conductivity_BTU_hr_ft_F * 1.7307  # W/(m·K) ultimate hardened
    
    # Weight fractions for variable Cp (Step 2)
    Cc = mix.total_cementitious_lb_yd3 * 0.593276  # kg/m³
    Wc = Cc  # cementitious weight per m³
    Ww = mix.water_lb_yd3 * 0.593276  # water weight per m³
    # Total aggregate = fine + coarse, convert lb/yd³ to kg/m³
    Wa = (mix.coarse_agg_lb_yd3 + mix.fine_agg_lb_yd3) * 0.593276
    Ca = mix.aggregate_Cp_BTU_lb_F * 4186.8  # aggregate Cp in J/(kg·K)
    
    # Hydration parameters
    Ea = mix.activation_energy_J_mol
    tau = mix.tau_hrs
    beta_h = mix.beta
    au = mix.alpha_u
    Hu = mix.Hu_J_kg * cal['Hu_factor']

    # ─── Soil properties ───
    soil = SOIL_PROPERTIES[soil_type]
    k_soil = soil["k"]
    rho_soil = soil["rho"]
    Cp_soil = soil["Cp"]
    
    # Deep ground temperature (CW Eq 44)
    T_aat_C = F_to_C(constr.avg_annual_temp_F)
    T_gw_C = 0.83 * T_aat_C + 3.7
    T_gw_K = T_gw_C + 273.15

    # ─── Blanket properties ───
    blanket_thickness = BLANKET_DEFAULTS["thickness_m"]
    rho_blanket = BLANKET_DEFAULTS["density"]
    Cp_blanket = BLANKET_DEFAULTS["Cp"]
    # Map user's physical R-value to effective R for 1D engine
    R_eff = blanket_r_effective(constr.blanket_R_value)
    R_blanket_SI = R_eff * 0.1761  # m²·K/W
    k_blanket = blanket_thickness / R_blanket_SI

    # ─── Grid ───
    dx_c = depth_m / (n_concrete - 1)
    dx_s = soil_depth_m / n_soil
    dx_b = blanket_thickness / n_blanket

    # Total nodes
    n_total = n_blanket + n_concrete + n_soil
    
    # Index boundaries
    ib_start = 0                           # blanket start
    ib_end = n_blanket                     # blanket end (exclusive) = concrete start
    ic_start = n_blanket
    ic_end = n_blanket + n_concrete
    is_start = n_blanket + n_concrete
    is_end = n_total

    # ─── CFL stability ───
    alpha_d_c = k_uc * 1.33 / (rho_c * Cp_c_base)  # worst case (early age, high k)
    alpha_d_s = k_soil / (rho_soil * Cp_soil)
    alpha_d_b = k_blanket / (rho_blanket * Cp_blanket)
    
    dt_s = dt_hrs * 3600.0
    dt_max_c = 0.4 * dx_c**2 / alpha_d_c
    dt_max_s = 0.4 * dx_s**2 / alpha_d_s
    dt_max_b = 0.4 * dx_b**2 / alpha_d_b
    dt_max = min(dt_max_c, dt_max_s, dt_max_b)
    
    if dt_s > dt_max:
        dt_hrs = dt_max / 3600.0
        dt_s = dt_max
    n_steps = int(duration_hrs / dt_hrs) + 1

    # ─── Outer convection coefficient ───
    wind = env.avg_max_wind_m_s * 0.4
    h_outer = 5.6 + 3.5 * wind + 5.5  # W/(m²·K)

    # ─── Initialize temperature arrays ───
    T_place_K = F_to_K(constr.placement_temp_F)
    T_soil_surface_K = F_to_K(constr.soil_temp_F)
    
    T = np.zeros(n_total)
    # Blanket: starts at placement temp (just placed on fresh concrete)
    T[ib_start:ib_end] = T_place_K
    # Concrete: uniform at placement temp
    T[ic_start:ic_end] = T_place_K
    # Soil: linear interpolation from soil_temp at surface to T_gw at depth
    for i in range(n_soil):
        frac = i / max(n_soil - 1, 1)
        T[is_start + i] = T_soil_surface_K + frac * (T_gw_K - T_soil_surface_K)

    # Equivalent age (only concrete nodes)
    te = np.full(n_concrete, 0.01)

    # Index for concrete core and surfaces
    mid_concrete = ic_start + n_concrete // 2

    # ─── Output sampling ───
    out_int = max(1, int(0.5 / dt_hrs))
    times, cores, tops, bots, ambs, dts = [], [], [], [], [], []

    # ─── Time-stepping loop ───
    for step in range(n_steps):
        t = step * dt_hrs
        T_amb_K = F_to_K(ambient_temp_F(t, env))

        # ── Hydration heat (concrete nodes only) ──
        T_concrete = T[ic_start:ic_end]
        af = arrhenius_vec(T_concrete, Ea)
        te += dt_hrs * af
        da_dt = hydration_rate_vec(te, tau, beta_h, au) * af
        Q_concrete = Hu * Cc * da_dt / 3600.0  # W/m³

        # Degree of hydration for variable properties (Step 2)
        if use_variable_properties:
            alpha_node = hydration_alpha_vec(te, tau, beta_h, au)
            T_C_concrete = T_concrete - 273.15
            k_c_nodes = thermal_conductivity_variable(k_uc, alpha_node)
            Cp_c_nodes = specific_heat_variable(Wc, Wa, Ww, Ca, alpha_node, T_C_concrete, rho_c)
        else:
            k_c_nodes = np.full(n_concrete, k_uc)
            Cp_c_nodes = np.full(n_concrete, Cp_c_base)

        # ── Build new temperature array ──
        Tn = T.copy()

        # ── BLANKET interior nodes ──
        if n_blanket > 2:
            d2T_b = (T[ib_start+2:ib_end] - 2*T[ib_start+1:ib_end-1] + T[ib_start:ib_end-2]) / dx_b**2
            Tn[ib_start+1:ib_end-1] = T[ib_start+1:ib_end-1] + (k_blanket * d2T_b) / (rho_blanket * Cp_blanket) * dt_s

        # ── CONCRETE interior nodes (Steps 1+2: vectorized with variable properties) ──
        if n_concrete > 2:
            d2T_c = (T[ic_start+2:ic_end] - 2*T[ic_start+1:ic_end-1] + T[ic_start:ic_end-2]) / dx_c**2
            k_interior = k_c_nodes[1:-1]
            Cp_interior = Cp_c_nodes[1:-1]
            Tn[ic_start+1:ic_end-1] = (T[ic_start+1:ic_end-1] 
                + (k_interior * d2T_c + Q_concrete[1:-1]) / (rho_c * Cp_interior) * dt_s)

        # ── SOIL interior nodes ──
        if n_soil > 2:
            d2T_s = (T[is_start+2:is_end] - 2*T[is_start+1:is_end-1] + T[is_start:is_end-2]) / dx_s**2
            Tn[is_start+1:is_end-1] = T[is_start+1:is_end-1] + (k_soil * d2T_s) / (rho_soil * Cp_soil) * dt_s

        # ── INTERFACE: Blanket top (node 0) — convection + evaporation to ambient ──
        # Two evaporation components:
        # 1. Menzel bleed water evaporation (decays in first ~6 hours)
        # 2. Constant wet-blanket evaporation (blanket stays moist, evaporates steadily)
        q_evap_total = blanket_evap_W_m2  # constant wet-blanket component
        if use_menzel_evap:
            T_surface_C = T[ib_start] - 273.15
            T_air_C = F_to_C(ambient_temp_F(t, env))
            q_evap_total += menzel_evaporation(t, T_surface_C, T_air_C, 
                                                env.relative_humidity, env.avg_max_wind_m_s * 0.4)
        
        # Blanket top node: conduction from below + convection to ambient - evaporation
        flux_b_top = (k_blanket * (T[ib_start+1] - T[ib_start]) / dx_b
                      + h_outer * (T_amb_K - T[ib_start])
                      - q_evap_total)
        Tn[ib_start] = T[ib_start] + flux_b_top / (rho_blanket * Cp_blanket * dx_b / 2) * dt_s

        # ── INTERFACE: Blanket-Concrete junction ──
        # Blanket bottom node (ib_end - 1) and Concrete top node (ic_start)
        # Use harmonic mean of k at interface
        k_bc = 2 * k_blanket * k_c_nodes[0] / (k_blanket + k_c_nodes[0])
        
        # Blanket bottom: conduction from blanket above + conduction to concrete below
        flux_b_bot = (k_blanket * (T[ib_end-2] - T[ib_end-1]) / dx_b
                      + k_bc * (T[ic_start] - T[ib_end-1]) / ((dx_b + dx_c) / 2))
        Tn[ib_end-1] = T[ib_end-1] + flux_b_bot / (rho_blanket * Cp_blanket * dx_b / 2) * dt_s
        
        # Concrete top: conduction from blanket above + conduction from concrete below + hydration
        flux_c_top = (k_bc * (T[ib_end-1] - T[ic_start]) / ((dx_b + dx_c) / 2)
                      + k_c_nodes[0] * (T[ic_start+1] - T[ic_start]) / dx_c
                      + Q_concrete[0] * dx_c / 2)
        Tn[ic_start] = T[ic_start] + flux_c_top / (rho_c * Cp_c_nodes[0] * dx_c / 2) * dt_s

        # ── INTERFACE: Concrete-Soil junction ──
        k_cs = 2 * k_c_nodes[-1] * k_soil / (k_c_nodes[-1] + k_soil)
        
        # Concrete bottom: conduction from concrete above + conduction to soil below + hydration
        flux_c_bot = (k_c_nodes[-1] * (T[ic_end-2] - T[ic_end-1]) / dx_c
                      + k_cs * (T[is_start] - T[ic_end-1]) / ((dx_c + dx_s) / 2)
                      + Q_concrete[-1] * dx_c / 2)
        Tn[ic_end-1] = T[ic_end-1] + flux_c_bot / (rho_c * Cp_c_nodes[-1] * dx_c / 2) * dt_s
        
        # Soil top: conduction from concrete above + conduction to soil below (no heat gen)
        flux_s_top = (k_cs * (T[ic_end-1] - T[is_start]) / ((dx_c + dx_s) / 2)
                      + k_soil * (T[is_start+1] - T[is_start]) / dx_s)
        Tn[is_start] = T[is_start] + flux_s_top / (rho_soil * Cp_soil * dx_s / 2) * dt_s

        # ── Bottom BC: Fixed deep ground temperature ──
        Tn[is_end - 1] = T_gw_K

        T = Tn

        # Sample output
        if step % out_int == 0:
            tc = K_to_F(T[mid_concrete])
            tt = K_to_F(T[ic_start])       # concrete top surface
            tb = K_to_F(T[ic_end - 1])     # concrete bottom surface
            ta = K_to_F(T_amb_K)
            dT_val = max(tc - min(tt, tb), 0)
            times.append(t)
            cores.append(tc)
            tops.append(tt)
            bots.append(tb)
            ambs.append(ta)
            dts.append(dT_val)

    times = np.array(times)
    cores = np.array(cores)
    tops = np.array(tops)
    bots = np.array(bots)
    ambs = np.array(ambs)
    dts = np.array(dts)

    pi = np.argmax(cores)
    di = np.argmax(dts)

    peak_ok = cores[pi] <= DOT_LIMITS['peak_temp_max_F']
    dt_ok = dts[di] <= DOT_LIMITS['delta_T_max_F']

    return {
        'time_hrs': times, 'T_core_F': cores, 'T_top_F': tops,
        'T_bot_F': bots, 'T_ambient_F': ambs, 'delta_T_F': dts,
        'peak_core_F': round(float(cores[pi]), 1),
        'peak_time_hrs': round(float(times[pi]), 1),
        'max_delta_T_F': round(float(dts[di]), 1),
        'max_delta_T_time_hrs': round(float(times[di]), 1),
        'peak_limit_F': DOT_LIMITS['peak_temp_max_F'],
        'delta_T_limit_F': DOT_LIMITS['delta_T_max_F'],
        'peak_pass': bool(peak_ok),
        'delta_T_pass': bool(dt_ok),
        'overall_pass': bool(peak_ok and dt_ok),
    }


# ============================================================
# TCP Section 3.2 Table Generator
# ============================================================
def generate_section_3_2(results: Dict, mix: MixDesign, constr: Construction) -> str:
    """Generate TCP Section 3.2 predicted results table."""
    s22 = section_2_2(mix)
    r = results
    lines = [
        "SECTION 3.2 — PREDICTED THERMAL ANALYSIS RESULTS",
        "=" * 55,
        "",
        f"Analysis Software:     CalcShore Thermal Engine v2.0",
        f"Analysis Method:       1D Finite Difference (explicit)",
        f"Placement Temperature: {constr.placement_temp_F}°F",
        f"Soil Temperature:      {constr.soil_temp_F}°F",
        f"Analysis Duration:     {r['time_hrs'][-1]:.0f} hours",
        "",
        "Mix Design Summary:",
        f"  Total Cementitious:  {s22['total_cementitious_lb_yd3']} lb/yd³",
        f"  w/cm Ratio:          {mix.wcm:.2f}",
        f"  ACI 207.2R Rise:     {s22['adiabatic_rise_F']}°F",
        "",
        "Predicted Results:",
        f"  Peak Core Temp:      {r['peak_core_F']}°F at {r['peak_time_hrs']} hours",
        f"  Max Differential:    {r['max_delta_T_F']}°F at {r['max_delta_T_time_hrs']} hours",
        "",
        "DOT Compliance:",
        f"  Peak Temperature:    {r['peak_core_F']}°F ≤ {r['peak_limit_F']}°F  "
        f"→ {'PASS ✅' if r['peak_pass'] else 'FAIL ❌'}",
        f"  Temperature Diff:    {r['max_delta_T_F']}°F ≤ {r['delta_T_limit_F']}°F  "
        f"→ {'PASS ✅' if r['delta_T_pass'] else 'FAIL ❌'}",
        "",
        f"  OVERALL:             {'PASS ✅' if r['overall_pass'] else 'FAIL ❌'}",
        "",
        "Temperature at Key Times:",
        f"  {'Time':>6}  {'Core':>7}  {'Top':>7}  {'Bottom':>7}  {'ΔT':>6}",
        f"  {'─'*40}",
    ]
    for t_target in [12, 24, 48, 72, 96, 120, 168]:
        idx = np.argmin(np.abs(r['time_hrs'] - t_target))
        lines.append(
            f"  {t_target:4d} hr  {r['T_core_F'][idx]:6.1f}°F"
            f"  {r['T_top_F'][idx]:6.1f}°F"
            f"  {r['T_bot_F'][idx]:6.1f}°F"
            f"  {r['delta_T_F'][idx]:5.1f}°F"
        )
    return "\n".join(lines)


# ============================================================
# Plotting
# ============================================================
def plot_results(results: Dict, title: str = "", save_path: Optional[str] = None):
    """Generate ConcreteWorks-style Max-Min temperature plot."""
    fig, ax1 = plt.subplots(figsize=(12, 7))
    t = results['time_hrs']
    ax1.plot(t, results['T_core_F'], 'b-', lw=2.5, label='Core Temperature')
    ax1.plot(t, results['T_top_F'], 'g-', lw=1.5, label='Top Surface')
    ax1.plot(t, results['T_bot_F'], 'm--', lw=1.2, label='Bottom Surface')
    ax1.plot(t, results['T_ambient_F'], color='#DAA520', lw=1, alpha=0.7, label='Ambient')
    ax1.axhline(160, color='r', ls=':', lw=1, alpha=0.4, label='Peak Limit (160°F)')
    ax1.set_xlabel('Time Since Placement (hours)', fontsize=12)
    ax1.set_ylabel('Temperature (°F)', fontsize=12, color='#333')
    ax1.set_ylim(40, 170)
    ax1.set_xlim(0, max(t))
    ax1.grid(True, alpha=0.2)

    ax2 = ax1.twinx()
    ax2.fill_between(t, results['delta_T_F'], alpha=0.08, color='red')
    ax2.plot(t, results['delta_T_F'], 'r-', lw=2, label='ΔT (Core − Surface)')
    ax2.axhline(35, color='r', ls='--', lw=1.5, alpha=0.6, label='ΔT Limit (35°F)')
    ax2.set_ylabel('Temperature Differential ΔT (°F)', fontsize=12, color='r')
    ax2.set_ylim(0, 60)

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc='upper left', fontsize=9, framealpha=0.9, edgecolor='#ccc')

    r = results
    sp = 'PASS' if r['peak_pass'] else 'FAIL'
    sd = 'PASS' if r['delta_T_pass'] else 'FAIL'
    txt = (f"Peak Core: {r['peak_core_F']}°F @ {r['peak_time_hrs']} hrs\n"
           f"Max dT:    {r['max_delta_T_F']}°F @ {r['max_delta_T_time_hrs']} hrs\n"
           f"Peak: {sp}  (<=160°F)\ndT:   {sd}  (<=35°F)")
    ax1.text(0.98, 0.97, txt, transform=ax1.transAxes, fontsize=10,
             va='top', ha='right', family='monospace',
             bbox=dict(boxstyle='round,pad=0.5', fc='#FFFBEF', ec='#DAA520', alpha=0.9))

    plt.title(title or 'Mass Concrete Thermal Analysis', fontsize=14, fontweight='bold', pad=15)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    return fig
