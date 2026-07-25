"""
US Standard Atmosphere 1976 model.

Provides ambient pressure, temperature, and density as a function of geometric
altitude. This is what your nozzle model calls to know what pressure the exhaust
is fighting against at any given altitude.

Model scope:
    - Valid from sea level (0 m) up to 86 km
    - We only implement the first two layers (troposphere + stratosphere lower)
      because a 1 kN LOX/Ethanol pressure-fed engine will never fly above ~30 km.
    - Layer 1 (troposphere):  0 to 11,000 m,  linear temperature lapse
    - Layer 2 (stratosphere): 11,000 to 20,000 m,  isothermal
    - Above 20 km: extrapolate the isothermal layer (adequate accuracy up to 25-30 km)

All quantities in SI:
    - Altitude: m (meters above sea level)
    - Pressure: Pa
    - Temperature: K
    - Density: kg/m^3

Reference: US Standard Atmosphere 1976, NASA/NOAA/USAF.
Sutton Appendix 2 tabulates the same values.
"""
import math
from dataclasses import dataclass


# ---------- constants ----------

# Sea-level reference values (US Std Atm 1976)
P_SEA_LEVEL_PA = 101325.0
T_SEA_LEVEL_K = 288.15
RHO_SEA_LEVEL = 1.225

# Gas properties of air
R_AIR = 287.053         # specific gas constant for dry air, J/(kg*K)
G0 = 9.80665            # standard gravity, m/s^2

# Layer 1 (troposphere)
LAPSE_RATE_TROPO = -0.0065   # K per meter (temperature decreases with altitude)
H_TROPOPAUSE = 11000.0       # altitude where troposphere ends, m

# Layer 2 (lower stratosphere) — isothermal
T_STRATOSPHERE = 216.65      # constant temperature 11-20 km, K
P_TROPOPAUSE_PA = 22632.06   # pressure at 11 km, Pa (from Layer 1 formula)


# ---------- data container ----------

@dataclass
class AtmosphereState:
    """Snapshot of ambient conditions at one altitude."""
    altitude_m: float
    pressure_pa: float
    temperature_k: float
    density_kgm3: float


# ---------- core functions ----------

def temperature_at_altitude(altitude_m: float) -> float:
    """
    Temperature (K) as a function of geometric altitude (m).

    Layer 1 (0 to 11 km):   T = T0 + lapse_rate * h
    Layer 2 (11+ km):        T = T_stratosphere (constant)
    """
    if altitude_m <= H_TROPOPAUSE:
        return T_SEA_LEVEL_K + LAPSE_RATE_TROPO * altitude_m
    else:
        return T_STRATOSPHERE


def pressure_at_altitude(altitude_m: float) -> float:
    """
    Ambient pressure (Pa) as a function of geometric altitude (m).

    Layer 1 (troposphere, linear lapse):
        P = P0 * (T / T0) ^ (-g / (R * lapse_rate))

    Layer 2 (stratosphere, isothermal):
        P = P_tropopause * exp( -g * (h - h_tropopause) / (R * T_stratosphere) )
    """
    temperature = temperature_at_altitude(altitude_m)

    if altitude_m <= H_TROPOPAUSE:
        pressure = P_SEA_LEVEL_PA * ( (temperature / T_SEA_LEVEL_K) ** (-G0 / (R_AIR * LAPSE_RATE_TROPO)))
        return pressure
    else:
        pressure = P_TROPOPAUSE_PA * math.exp( -G0 * (altitude_m - H_TROPOPAUSE) / (R_AIR * T_STRATOSPHERE ))
        return pressure

def density_at_altitude(altitude_m: float, pressure_pa: float, temperature_k: float) -> float:
    """
    Density (kg/m^3) from ideal gas law: rho = P / (R * T).
    """
    return pressure_pa / (R_AIR * temperature_k)


def get_state(altitude_m: float):
    """Convenience: return all three at once."""
    t = temperature_at_altitude(altitude_m)
    p = pressure_at_altitude(altitude_m)
    rho = density_at_altitude(altitude_m, p, t)
    return AtmosphereState(altitude_m = altitude_m, pressure_pa = p, temperature_k = t, density_kgm3 = rho)


# ---------- Test ----------

if __name__ == "__main__":
    print(f"{'Altitude':>10} {'T (K)':>10} {'P (kPa)':>12} {'rho (kg/m^3)':>14}")
    for h in [0, 1000, 5000, 11000, 15000, 20000, 25000, 30000]:
        state = get_state(h)
        print(f"{state.altitude_m:>10.0f} "
              f"{state.temperature_k:>10.2f} "
              f"{state.pressure_pa/1000:>12.3f} "
              f"{state.density_kgm3:>14.5f}")