"""
Altitude sweep: engine performance vs altitude.

Fixes the engine geometry and chamber conditions, sweeps altitude from
sea level to 30 km, and plots how thrust and Isp change with ambient
pressure.

This is the first plot for the Week 2 deliverable and the first script
where propellants.py, nozzle.py, and atmosphere.py talk to each other.
"""

import numpy as np
import matplotlib.pyplot as plt

from propellants import Propellant
from nozzle import IsentropicNozzle
from atmosphere import pressure_at_altitude

# ---------- design point (fixed for this sweep) ----------

PC_PA = 20e5              # 20 bar chamber pressure
OF_RATIO = 1.5            # LOX/Ethanol mixture ratio
EXPANSION_RATIO = 8.0     # geometric area ratio A_e/A_t
THROAT_AREA_M2 = 1e-4     # 100 mm^2 throat, placeholder

# Altitude sweep parameters
H_MIN_M = 0
H_MAX_M = 30000
H_STEP_M = 500

G0 = 9.80665              # standard gravity, m/s^2 (for Isp)

# ---------- setup ----------
# Build a Propellant, get the chamber CombustionState.
# Build an IsentropicNozzle from that chamber state.
prop = Propellant(ox = "LOX", fuel = "Ethanol")
chamber = prop.get_combustion(pc_pa = PC_PA, of = OF_RATIO, eps = EXPANSION_RATIO)
nozzle = IsentropicNozzle(gamma = chamber.gamma, mw = chamber.mw, tc_k = chamber.tc_k, pc_pa = PC_PA)

# Altitude array: sea level to 30 km in 500 m steps.
altitude_m = np.arange(H_MIN_M, H_MAX_M + H_STEP_M, H_STEP_M)

# Collect results across the sweep. Kept as plain lists during the loop
thrusts_n = []
isps_s = []
ambient_pressure_pa = []
exit_pressure_pa = []

# Main sweep
for h in altitude_m:
    ambient_pa = pressure_at_altitude(h)
    state = nozzle.solve( area_throat_m2 = THROAT_AREA_M2, expansion_ratio = EXPANSION_RATIO,
              ambient_pa = ambient_pa , cstar_ms = chamber.cstar_ms)
    isp = state.thrust_n / (state.mdot_kgs * G0)

    thrusts_n.append(state.thrust_n)
    isps_s.append(isp)
    ambient_pressure_pa.append(ambient_pa)
    exit_pressure_pa.append(state.pe_pa)

# Convert to km for x-axis readability.
altitudes_km = np.array(altitude_m) / 1000.0

# Graphing altitude vs thrust
fig, ax = plt.subplots()
ax.plot(altitudes_km, thrusts_n)
ax.set_xlabel("Altitude (km)")
ax.set_ylabel("Thrust (N)")
ax.set_title(f"Thrust vs Altitude — LOX/Ethanol, Pc={PC_PA/1e5:.0f} bar, ε={EXPANSION_RATIO:.0f}")
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig("plots/thrust_vs_altitude.png", dpi=150)

#Graphing isp vs altitude
fig, ax = plt.subplots()
ax.plot(altitudes_km, isps_s)
ax.set_xlabel("Altitude (km)")
ax.set_ylabel("Specific Impulse (s)")
ax.set_title(f"Isp vs Altitude — LOX/Ethanol, Pc={PC_PA/1e5:.0f} bar, ε={EXPANSION_RATIO:.0f}")
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig("plots/isp_vs_altitude.png", dpi=150)

#Graphing ambient/exit pressure vs altitude
fig, ax = plt.subplots()
ax.plot(altitudes_km, ambient_pressure_pa, label = "Ambient")
ax.plot(altitudes_km, exit_pressure_pa, label = "Exit")
ax.legend()
ax.set_xlabel("Altitude (km)")
ax.set_ylabel("Pressure (Pa)")
ax.set_title("Exit and Ambient Pressure vs Altitude")
ax.set_yscale("log")
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig("plots/pressure_vs_altitude.png", dpi=150)

plt.show()
