import numpy as np
import matplotlib.pyplot as plt

from size_engine import size_engine

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

# Plotting pressure vs throat diameter
fig, ax = plt.subplots()
ax.plot(pc_range_bar, diameter_throat_mm)
ax.set_xlabel("Chamber Pressure (bar)")
ax.set_ylabel("Throat Diameter (mm)")
ax.set_title("Chamber Pressure vs Throat Diameter")
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig("plots/sensitivity_throat_vs_pc.png", dpi=150)

# Plotting pressure vs Specific Impulse
fig, ax = plt.subplots()
ax.plot(pc_range_bar, isps_s)
ax.set_xlabel("Chamber Pressure (bar)")
ax.set_ylabel("Specific Impulse (s)")
ax.set_title("Chamber Pressure vs Specific Impulse")
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig("plots/sensitivity_isp_vs_pc.png", dpi=150)

# Plotting pressure vs Tank Pressure
fig, ax = plt.subplots()
ax.plot(pc_range_bar, tank_pressure_bar)
ax.set_xlabel("Chamber Pressure (bar)")
ax.set_ylabel("Tank Pressure (bar)")
ax.set_title("Chamber Pressure vs Tank Pressure")
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig("plots/sensitivity_tank_vs_pc.png", dpi=150)

plt.show()