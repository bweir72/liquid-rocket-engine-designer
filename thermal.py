"""
Thermal analysis: peak convective heat flux at the throat.

Computes the gas-side heat transfer coefficient via the Bartz correlation
and the resulting heat flux from combustion gas to the chamber wall. This
is the peak-flux location; if a chamber survives here, it survives everywhere.

Assumes:
    - Ablative or heat-sink chamber wall (uncooled)
    - Wall temperature is a design input (fixed operating condition)
    - Throat radius of curvature approximated as r_c ≈ D_t / 2
    - Local station factor (A_t/A)^0.9 = 1 (throat-only analysis)
"""

from dataclasses import dataclass
from size_engine import EngineDesign

TWALL_K = 500.0

@dataclass
class ThermalState:
    """Snapshot of thermal conditions at the throat."""
    hg_wm2k: float          # gas-side heat transfer coefficient, W/(m²·K)
    taw_k: float            # adiabatic wall temperature, K
    twall_k: float          # assumed wall temperature (input), K
    q_wm2: float            # heat flux, W/m²
    sigma_correction: float # Bartz wall-temp correction, dimensionless

def compute_throat_heat_flux(engine: EngineDesign, twall_k = TWALL_K) -> ThermalState:
    """
    Compute peak convective heat flux at the throat via Bartz.
    
    twall_k: assumed wall temperature. 500 K is representative of an
             ablative or heat-sink chamber at burn onset.
    """
    pc = engine.pc_pa
    tc = engine.tc_k
    gamma = engine.gamma
    cstar = engine.cstar_ms
    dt = engine.diameter_throat_m
    mu = engine.mu_kgms
    cp = engine.cp_j_kg_k
    pr = engine.prandtl
    
    r = pr ** (1.0 / 3.0)                                           # Recovery factor and adiabatic wall temperature at throat (M=1)
    taw = tc * (1 + r * (gamma - 1) / 2) / (1 + (gamma - 1) / 2)    # Wall-temperature correction sigma

    #inner1 = 0.5 * (twall_k / tc) * (1 + (gamma - 1) / 2) + 0.5
    #inner2 = 1 + (gamma - 1) / 2
    sigma = 1.0 / (0.5 * (twall_k / tc + 1)) ** 0.68                  # Sigma (Wall temperature correction)
    
    A_factor = 0.026 / (dt ** 0.2)
    B_factor = cp * (mu ** 0.2) / (pr ** 0.6)
    C_factor = (pc / cstar) ** 0.8
    hg = A_factor * B_factor * C_factor * sigma                     # Bartz h_g (SI form, curvature term dropped, σ handles wall correction)

    
    q = hg * (taw - twall_k)                                        # Heat flux
    
    return ThermalState(
        hg_wm2k = hg,
        taw_k = taw,
        twall_k = twall_k,
        q_wm2 = q,
        sigma_correction = sigma
    )


if __name__ == "__main__":
    from size_engine import size_engine
    
    engine = size_engine(
        thrust_target_n=1000,
        altitude_m=0,
        pc_pa=20e5,
        of_ratio=1.5,
    )
    
    thermal = compute_throat_heat_flux(engine, twall_k=500.0)
    
    print(f"Thermal analysis at throat (T_wall = {thermal.twall_k:.0f} K):")
    print(f"  σ (Bartz corr.) = {thermal.sigma_correction:.3f}")
    print(f"  T_aw            = {thermal.taw_k:.1f} K")
    print(f"  h_g             = {thermal.hg_wm2k/1000:.2f} kW/(m²·K)")
    print(f"  q               = {thermal.q_wm2/1e6:.2f} MW/m²")