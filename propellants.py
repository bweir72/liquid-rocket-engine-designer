"""
Propellant chemistry wrapper around NASA CEA (via rocketcea).

Everything in this module works in SI units:
  - Pressure: Pa
  - Temperature: K
  - Velocity: m/s

CEA internally uses imperial (psia, Rankine, ft/s). We convert at the boundary
so the rest of the codebase never has to think about it.
"""

from dataclasses import dataclass
from rocketcea.cea_obj import CEA_Obj

# Unit conversions
PA_PER_PSI = 6894.757
FTS_PER_MS = 1.0 / 0.3048  # divide m/s by this to get ft/s; multiply ft/s by this... wait
                            # actually: 1 m = 3.28084 ft, so 1 m/s = 3.28084 ft/s
                            # so m/s = ft/s / 3.28084
FT_PER_M = 3.28084


@dataclass
class CombustionState:
    """Snapshot of chamber combustion products at a given operating point."""
    tc_k: float          # chamber temperature, Kelvin
    cstar_ms: float      # characteristic velocity, m/s
    gamma: float         # ratio of specific heats (dimensionless)
    mw: float            # molecular weight of combustion gases, kg/kmol
    isp_ideal_s: float   # ideal vacuum-ish Isp from CEA, seconds


class Propellant:
    """LOX/Ethanol combustion model wrapping rocketcea."""

    def __init__(self, ox: str = "LOX", fuel: str = "Ethanol"):
        self.ox = ox
        self.fuel = fuel
        self._cea = CEA_Obj(oxName=ox, fuelName=fuel)

    def get_combustion(self, pc_pa: float, of: float, eps: float = 8.0) -> CombustionState:
        """
        Returns combustion chamber state at the given operating point.

        pc_pa : chamber pressure in Pa
        of    : oxidizer-to-fuel mass ratio (dimensionless)
        eps   : nozzle area expansion ratio (dimensionless). Only affects Isp,
                not the chamber-side properties.
        """
        pc_psia = pc_pa / PA_PER_PSI

        tc_rankine = self._cea.get_Tcomb(Pc=pc_psia, MR=of)
        cstar_fts = self._cea.get_Cstar(Pc=pc_psia, MR=of)
        gamma = self._cea.get_Chamber_MolWt_gamma(Pc=pc_psia, MR=of, eps=eps)[1]
        mw = self._cea.get_Chamber_MolWt_gamma(Pc=pc_psia, MR=of, eps=eps)[0]
        isp = self._cea.get_Isp(Pc=pc_psia, MR=of, eps=eps)

        return CombustionState(
            tc_k=tc_rankine * 5.0 / 9.0,
            cstar_ms=cstar_fts / FT_PER_M,
            gamma=gamma,
            mw=mw,
            isp_ideal_s=isp,
        )


if __name__ == "__main__":
    # Quick smoke test — run this file directly to sanity-check CEA
    prop = Propellant()
    state = prop.get_combustion(pc_pa=20e5, of=1.5, eps=8.0)
    print(f"LOX/Ethanol at Pc=20 bar, O/F=1.5, eps=8:")
    print(f"  Tc     = {state.tc_k:.1f} K")
    print(f"  c*     = {state.cstar_ms:.1f} m/s")
    print(f"  gamma  = {state.gamma:.3f}")
    print(f"  MW     = {state.mw:.2f} kg/kmol")
    print(f"  Isp    = {state.isp_ideal_s:.1f} s (ideal)")