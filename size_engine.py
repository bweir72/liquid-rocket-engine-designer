"""
Engine sizing tool: inverse design from thrust target to geometry.

Given a thrust target, design altitude, propellant, chamber pressure, and O/F,
returns the throat area, expansion ratio, mass flow, tank pressures, and Isp
for a fully-closed pressure-fed engine design.

This is the inverse of nozzle.solve(): instead of "given geometry, what's the
performance?", it answers "given the mission, what geometry?"
"""

import math
from dataclasses import dataclass
from scipy.optimize import brentq

from propellants import Propellant
from nozzle import IsentropicNozzle
from atmosphere import pressure_at_altitude


G0 = 9.80665   # standard gravity, m/s^2


@dataclass
class EngineDesign:
    """Complete pressure-fed engine design at one operating point."""
    thrust_target_n: float
    altitude_m: float
    pc_pa: float
    of_ratio: float

    area_throat_m2: float
    diameter_throat_m: float
    area_exit_m2: float
    diameter_exit_m: float
    expansion_ratio: float

    thrust_actual_n: float
    isp_s: float
    mdot_kgs: float
    mdot_ox_kgs: float
    mdot_fuel_kgs: float

    tank_pressure_pa: float

    tc_k: float          # chamber temperature, K
    cstar_ms: float      # characteristic velocity, m/s
    gamma: float         # ratio of specific heats

    mu_kgms: float
    cp_j_kg_k: float
    prandtl: float


def size_engine(
    thrust_target_n: float,
    altitude_m: float,
    pc_pa: float,
    of_ratio: float,
    ox: str = "LOX",
    fuel: str = "Ethanol",
    injector_dp_frac: float = 0.25,
    line_losses_frac: float = 0.05,
    ullage_margin_pa: float = 2e5,
) -> EngineDesign:
    """
    Size a pressure-fed engine to hit the thrust target at the design altitude.
    """

    # ---------- setup ----------
    # Get ambient pressure at design altitude.
    # Build a Propellant, get the chamber CombustionState.
    # Build an IsentropicNozzle from that chamber state.

    ambient_pa = pressure_at_altitude(altitude_m)
    prop = Propellant(ox, fuel)
    chamber = prop.get_combustion(pc_pa = pc_pa , of = of_ratio , eps = 8.0) # eps here only affects Isp field; not used downstream
    nozzle = IsentropicNozzle( gamma = chamber.gamma, mw = chamber.mw, tc_k = chamber.tc_k, pc_pa = pc_pa)


    # ---------- STEP 1: matched-exit expansion ratio ----------
    # Ask the nozzle for the ε that matches exit pressure to ambient.

    epsilon = nozzle.expansion_ratio_for_matched_exit(ambient_pa)


    # ---------- STEP 2: find throat area for target thrust ----------
    # Define a residual function that takes a candidate throat area 
    # that uses brentq to find the root where ther residual = 0
    
    def thrust_residual(at_m2):         
        nozzle_state = nozzle.solve(area_throat_m2 = at_m2, expansion_ratio = epsilon, ambient_pa = ambient_pa, cstar_ms = chamber.cstar_ms)
        return nozzle_state.thrust_n - thrust_target_n
        
    area_lower_bound = 1e-6   # 
    area_upper_bound = 1e-2   # 
    area_throat_m2 = brentq(thrust_residual, area_lower_bound, area_upper_bound)

    # ---------- STEP 3: re-run solve at the final geometry ----------

    final_state = nozzle.solve(area_throat_m2 = area_throat_m2, expansion_ratio = epsilon, ambient_pa = ambient_pa, cstar_ms = chamber.cstar_ms)


    # ---------- STEP 4: derived quantities ----------
    # From area_throat_m2, epsilon, final_state, chamber, and the inputs,
    # compute:
    #   diameter_throat_m  — throat area is a circle, convert A to D
    #   area_exit_m2       — ε times throat area
    #   diameter_exit_m    — same circle formula, for exit area
    #   isp_s              — thrust ÷ (mdot · g₀)
    #   mdot_ox_kgs        — total mdot split by O/F ratio (oxidizer share)
    #   mdot_fuel_kgs      — fuel share (remainder)
    #   tank_pressure_pa   — Pc scaled up for injector ΔP + line losses,
    #                        plus the ullage margin

    diameter_throat_m = 2 * math.sqrt(area_throat_m2 / math.pi)
    area_exit_m2 = epsilon * area_throat_m2
    diameter_exit_m = 2 * math.sqrt(area_exit_m2 / math.pi)
    isp_s = final_state.thrust_n / (final_state.mdot_kgs * G0)
    mdot_ox_kgs = final_state.mdot_kgs * of_ratio / (1 + of_ratio)
    mdot_fuel_kgs = final_state.mdot_kgs / (1 + of_ratio)
    tank_pressure_pa = pc_pa * (1 + injector_dp_frac + line_losses_frac) + ullage_margin_pa


    # ---------- assemble and return ----------

    return EngineDesign(
        thrust_target_n = thrust_target_n,
        altitude_m = altitude_m,
        pc_pa = pc_pa,
        of_ratio = of_ratio,
        area_throat_m2 = area_throat_m2,
        diameter_throat_m = diameter_throat_m,
        area_exit_m2 = area_exit_m2,
        diameter_exit_m = diameter_exit_m,
        expansion_ratio = epsilon,
        thrust_actual_n = final_state.thrust_n,
        isp_s = isp_s,
        mdot_kgs = final_state.mdot_kgs,
        mdot_ox_kgs = mdot_ox_kgs,
        mdot_fuel_kgs = mdot_fuel_kgs,
        tank_pressure_pa = tank_pressure_pa,
        tc_k = chamber.tc_k,
        cstar_ms = chamber.cstar_ms,
        gamma = chamber.gamma,
        mu_kgms=chamber.mu_kgms,
        cp_j_kg_k=chamber.cp_j_kg_k,
        prandtl=chamber.prandtl,
    )


# ---------- smoke test ----------

if __name__ == "__main__":
    design = size_engine(
        thrust_target_n=1000,
        altitude_m=0,
        pc_pa=20e5,
        of_ratio=1.5,
    )

    print(f"Engine sizing — {design.thrust_target_n:.0f} N target at {design.altitude_m/1000:.0f} km altitude:")
    print(f"  ε               = {design.expansion_ratio:.2f}")
    print(f"  A_t             = {design.area_throat_m2*1e6:.2f} mm²")
    print(f"  D_t             = {design.diameter_throat_m*1000:.2f} mm")
    print(f"  A_e             = {design.area_exit_m2*1e6:.2f} mm²")
    print(f"  D_e             = {design.diameter_exit_m*1000:.2f} mm")
    print(f"  Thrust (actual) = {design.thrust_actual_n:.1f} N")
    print(f"  Isp             = {design.isp_s:.1f} s")
    print(f"  mdot            = {design.mdot_kgs*1000:.1f} g/s")
    print(f"  ox              = {design.mdot_ox_kgs*1000:.1f} g/s")
    print(f"  fuel            = {design.mdot_fuel_kgs*1000:.1f} g/s")
    print(f"  Tank pressure   = {design.tank_pressure_pa/1e5:.1f} bar")