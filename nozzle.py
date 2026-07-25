"""
Isentropic nozzle flow calculator.

Given a combustion chamber state (Tc, gamma, MW, c*) and nozzle geometry
(throat area, expansion ratio), compute exit conditions and performance:
    - Exit Mach number
    - Exit pressure
    - Exit velocity
    - Thrust coefficient (Cf)
    - Mass flow rate
    - Thrust

All quantities in SI:
    - Pressure: Pa
    - Temperature: K
    - Area: m^2
    - Velocity: m/s
    - Mass flow: kg/s
    - Thrust: N

Assumes:
    - Ideal, quasi-1D, steady, isentropic flow
    - Choked at the throat (Pc / P_ambient large enough for supersonic exit)
    - Constant gamma through the nozzle (frozen composition — a simplification)
    - Perfect gas
"""

from dataclasses import dataclass
from scipy.optimize import brentq
import math

# Universal gas constant, J / (kmol * K)
R_UNIVERSAL = 8314.462

# ---------- data container ----------
@dataclass
class NozzleState:
    """Snapshot of nozzle exit conditions and performance at one operating point."""
    mach_e: float        # exit Mach number
    pe_pa: float         # exit pressure, Pa
    te_k: float          # exit temperature, K
    ve_ms: float         # exit velocity, m/s
    cf: float            # thrust coefficient (dimensionless)
    mdot_kgs: float      # mass flow rate, kg/s
    thrust_n: float      # thrust, N
    
class IsentropicNozzle:
    """
    Quasi-1D isentropic nozzle model.

    Parameters
    ----------
    gamma : float
        Ratio of specific heats of combustion products (dimensionless).
    mw : float
        Molecular weight of combustion products, kg/kmol.
    tc_k : float
        Chamber (stagnation) temperature, K.
    pc_pa : float
        Chamber (stagnation) pressure, Pa.
    """

    def __init__(self, gamma: float, mw: float, tc_k: float, pc_pa: float):
        self.gamma = gamma
        self.mw = mw
        self.tc = tc_k
        self.pc = pc_pa
        self.R = R_UNIVERSAL / mw   # specific gas constant, J/(kg*K)

    # ---------- core isentropic relations ----------

    def area_ratio_from_mach(self, mach: float) -> float:
        """
        Area ratio A/A* as a function of local Mach number.
        Sutton Eq. 3-14:
            A/A* = (1/M) * [ (2/(k+1)) * (1 + (k-1)/2 * M^2) ] ^ ((k+1)/(2(k-1)))
        """
        k = self.gamma
        inner = 1.0 + 0.5 * (k - 1.0) * mach**2            # (1 + (k-1)/2 * M^2)
        bracket = (2.0 / (k + 1.0)) * inner                # [ 2/(k+1) * inner ]
        exponent = (k + 1.0) / (2.0 * (k - 1.0))           # (k+1) / (2(k-1))
        ratio = (1.0 / mach) * bracket**exponent           # A/A*
        return ratio

    def mach_from_area_ratio(self, area_ratio: float) -> float:
        """
        Solve the area-Mach relation numerically for the *supersonic* branch.
        """
        if area_ratio < 1.0:
            raise ValueError(f"Area ratio must be >= 1, got {area_ratio}")

        # brentq finds the M where (area_ratio_from_mach(M) - area_ratio) = 0
        def residual(m):
            return self.area_ratio_from_mach(m) - area_ratio

        mach_lower_bound = 1.0001   # just above sonic (supersonic branch)
        mach_upper_bound = 50.0     # comfortably above any realistic exit Mach
        mach = brentq(residual, mach_lower_bound, mach_upper_bound)
        return mach

    def pressure_ratio(self, mach: float) -> float:
        """
        Static-to-stagnation pressure ratio P/Pc.
        Sutton Eq. 3-13:
            P/Pc = (1 + (k-1)/2 * M^2) ^ (-k/(k-1))
        """
        k = self.gamma
        inner = 1.0 + 0.5 * (k - 1.0) * mach**2   # (1 + (k-1)/2 * M^2)
        exponent = -k / (k - 1.0)                 # -k/(k-1)
        ratio = inner ** exponent                 # P/Pc
        return ratio

    def temperature_ratio(self, mach: float) -> float:
        """
        Static-to-stagnation temperature ratio T/Tc.
        Sutton Eq. 3-12:
            T/Tc = 1 / (1 + (k-1)/2 * M^2)
        """
        k = self.gamma
        inner = 1.0 + 0.5 * (k - 1.0) * mach**2   # (1 + (k-1)/2 * M^2)
        ratio = 1.0 / inner                       # T/Tc
        return ratio
    
    # ---------- solver: geometry + chamber state → performance ----------

    def solve(self, area_throat_m2: float, expansion_ratio: float,
              ambient_pa: float, cstar_ms: float) -> NozzleState:
        """
        Given throat area, expansion ratio, ambient pressure, and c*,
        compute all exit conditions and thrust.

        Steps:
          1. Solve area-Mach relation for exit Mach (M_e)
          2. Compute exit pressure P_e and exit temperature T_e
          3. Compute exit velocity V_e
          4. Compute mass flow rate mdot
          5. Compute thrust coefficient Cf and thrust F
        """
        mach_e = self.mach_from_area_ratio(expansion_ratio)
        pe_pa = self.pc * self.pressure_ratio(mach_e)
        te_k = self.tc * self.temperature_ratio(mach_e)
        ve_ms = mach_e * math.sqrt(self.gamma * self.R * te_k)
        mdot_kgs = (self.pc * area_throat_m2) / cstar_ms
        area_exit_m2 = area_throat_m2 * expansion_ratio
        thrust_n = mdot_kgs * ve_ms + (pe_pa - ambient_pa) * area_exit_m2
        cf = thrust_n / (self.pc * area_throat_m2)  # dimensionless

        return NozzleState(
            mach_e = mach_e,
            pe_pa = pe_pa,
            te_k = te_k,
            ve_ms = ve_ms,
            cf = cf,
            mdot_kgs = mdot_kgs,
            thrust_n =  thrust_n,
        )
    

if __name__ == "__main__":
    # Smoke test — compute performance for a plausible 1 kN LOX/Ethanol design point.
    nozzle = IsentropicNozzle(
        gamma=1.135,
        mw=22.17,
        tc_k=3189.1,
        pc_pa=20e5,   # 20 bar chamber pressure
    )

    state = nozzle.solve(
        area_throat_m2=1e-4,   # 100 mm^2 throat (~11.3 mm diameter)
        expansion_ratio=8.0,
        ambient_pa=101325,     # sea level
        cstar_ms=1720.3,
    )

    print("Nozzle exit conditions (sea level, eps=8):")
    print(f"  M_e    = {state.mach_e:.3f}")
    print(f"  P_e    = {state.pe_pa/1000:.2f} kPa   (ambient = {101325/1000:.2f} kPa)")
    print(f"  T_e    = {state.te_k:.1f} K")
    print(f"  V_e    = {state.ve_ms:.1f} m/s")
    print(f"  mdot   = {state.mdot_kgs*1000:.2f} g/s")
    print(f"  Cf     = {state.cf:.3f}")
    print(f"  Thrust = {state.thrust_n:.1f} N")
