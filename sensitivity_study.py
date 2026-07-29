"""
Sensitivity study: how the engine design responds to each design knob.

For a 1 kN pressure-fed LOX/Ethanol engine sized to sea level, this script
sweeps each of the main design parameters (chamber pressure, O/F ratio) one
at a time while holding the others fixed at the reference design point.
Each sweep produces three plots showing the tradeoffs that justify the
final design choices in the report.
"""


#------------REFACTOR THIS LATER IF THERE'S TIME-------------------


import numpy as np
import matplotlib.pyplot as plt
from size_engine import size_engine

def sweep_chamber_pressure():
    #---------- Fixed Parameters ---------- 
    thrust_target_n = 1000
    altitude_m = 0
    of_ratio = 1.5

    #---------- Varied Parameter ---------- 
    pc_range_pa = np.linspace(10e5, 40e5, 31)

    # ---------- Setup ----------
    diameter_throat_m = []
    isps_s = []
    tank_pressure_pa = []

    # Main sweep
    for pc in pc_range_pa:
        engine = size_engine(thrust_target_n = thrust_target_n,
                            altitude_m = altitude_m,
                            pc_pa = pc,
                            of_ratio = of_ratio
        )
        diameter_throat_m.append(engine.diameter_throat_m)
        isps_s.append(engine.isp_s)
        tank_pressure_pa.append(engine.tank_pressure_pa)

    # Converting pressure (pa) to bar
    pc_range_bar = pc_range_pa / 1e5

    # Converting to diamter to mm 
    diameter_throat_mm = np.array(diameter_throat_m) * 1000

    # Converting tank pressure to bar
    tank_pressure_bar = np.array(tank_pressure_pa) / 1e5

    # Plotting throat pressure vs throat diameter
    fig, ax = plt.subplots()
    ax.plot(pc_range_bar, diameter_throat_mm)
    ax.set_xlabel("Chamber Pressure (bar)")
    ax.set_ylabel("Throat Diameter (mm)")
    ax.set_title("Chamber Pressure vs Throat Diameter")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig("plots/sensitivity_pc/throat_dia_vs_pc.png", dpi=150)

    # Plotting throat pressure vs Specific Impulse
    fig, ax = plt.subplots()
    ax.plot(pc_range_bar, isps_s)
    ax.set_xlabel("Chamber Pressure (bar)")
    ax.set_ylabel("Specific Impulse (s)")
    ax.set_title("Chamber Pressure vs Specific Impulse")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig("plots/sensitivity_pc/isp_vs_pc.png", dpi=150)

    # Plotting throat pressure vs Tank Pressure
    fig, ax = plt.subplots()
    ax.plot(pc_range_bar, tank_pressure_bar)
    ax.set_xlabel("Chamber Pressure (bar)")
    ax.set_ylabel("Tank Pressure (bar)")
    ax.set_title("Chamber Pressure vs Tank Pressure")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig("plots/sensitivity_pc/tank_pressure_vs_pc.png", dpi=150)



def sweep_of_ratio():
#---------- Fixed Parameters ---------- 
    thrust_target_n = 1000
    altitude_m = 0
    pc_pa = 20e5        # Fixed at 20 bar because it was the sweet spot in chamber pressure sweep

    #---------- Varied Parameter ---------- 
    of_range = np.linspace(1.0, 2.0, 31)

    # ---------- Setup ----------
    diameter_throat_m = []
    isps_s = []
    tc_k = []

    # Main sweep
    for of in of_range:
        engine = size_engine(thrust_target_n = thrust_target_n,
                            altitude_m = altitude_m,
                            pc_pa = pc_pa,
                            of_ratio = of
        )
        diameter_throat_m.append(engine.diameter_throat_m)
        isps_s.append(engine.isp_s)
        tc_k.append(engine.tc_k)

    # Converting to diamter to mm 
    diameter_throat_mm = np.array(diameter_throat_m) * 1000
    
    # Plotting o/f ratio vs throat diameter
    fig, ax = plt.subplots()
    ax.plot(of_range, diameter_throat_mm)
    ax.set_xlabel("O/F Ratio ")
    ax.set_ylabel("Throat Diameter (mm)")
    ax.set_title("O/F Ratio vs Throat Diameter")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig("plots/sensitivity_of/throat_dia_vs_of.png", dpi=150)

    # Plotting o/f ratio vs Specific Impulse
    fig, ax = plt.subplots()
    ax.plot(of_range, isps_s)
    ax.set_xlabel("O/F Ratio")
    ax.set_ylabel("Specific Impulse (s)")
    ax.set_title("O/F Ratio vs Specific Impulse")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig("plots/sensitivity_of/isp_vs_of.png", dpi=150)

    # Plotting o/f ratio vs chamber temp
    fig, ax = plt.subplots()
    ax.plot(of_range, tc_k)
    ax.set_xlabel("O/F Ratio")
    ax.set_ylabel("Chamber Temp (K)")
    ax.set_title("O/F Ratio vs Chamber Temp")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig("plots/sensitivity_of/tc_vs_of.png", dpi=150)

def sweep_altitude():
#---------- Fixed Parameters ---------- 
    thrust_target_n = 1000
    of_ratio = 1.5
    pc_pa = 20e5        # Fixed at 20 bar because it was the sweet spot in chamber pressure sweep

    #---------- Varied Parameter ---------- 
    altitude_range_m = np.linspace(0.0, 3.0e4, 31)

    # ---------- Setup ----------
    diameter_exit_m = []
    isps_s = []
    expansion_ratio = []

    # Main sweep
    for h in altitude_range_m:
        engine = size_engine(thrust_target_n = thrust_target_n,
                            altitude_m = h,
                            pc_pa = pc_pa,
                            of_ratio = of_ratio
        )
        diameter_exit_m.append(engine.diameter_exit_m)
        isps_s.append(engine.isp_s)
        expansion_ratio.append(engine.expansion_ratio)

    # Converting diamter to mm 
    diameter_exit_mm = np.array(diameter_exit_m) * 1000

    # Converting altitude to km
    altitude_range_km = np.array(altitude_range_m) / 1000
    
    # Plotting altitude vs exit diameter
    fig, ax = plt.subplots()
    ax.plot(altitude_range_km, diameter_exit_mm)
    ax.set_xlabel("Altitude (km)")
    ax.set_ylabel("Exit Diameter (mm)")
    ax.set_title("Altitude vs Exit Diameter")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig("plots/sensitivity_altitude/exit_dia_vs_altitude.png", dpi=150)

    # Plotting altitude vs Specific Impulse
    fig, ax = plt.subplots()
    ax.plot(altitude_range_km, isps_s)
    ax.set_xlabel("Altitude (km)")
    ax.set_ylabel("Specific Impulse (s)")
    ax.set_title("Altitude vs Specific Impulse")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig("plots/sensitivity_altitude/isp_vs_altitude.png", dpi=150)

    # Plotting altitude vs expansion ratio
    fig, ax = plt.subplots()
    ax.plot(altitude_range_km, expansion_ratio)
    ax.set_xlabel("Altitude (km)")
    ax.set_ylabel("Expansion Ratio")
    ax.set_title("Altitude vs Expansion Ratio")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig("plots/sensitivity_altitude/expansion_ratio_vs_altitude.png", dpi=150)


if __name__ == "__main__":
    sweep_chamber_pressure()
    sweep_of_ratio()
    sweep_altitude()
    plt.show()