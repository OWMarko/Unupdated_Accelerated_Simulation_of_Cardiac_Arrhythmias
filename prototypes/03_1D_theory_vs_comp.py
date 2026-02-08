import numpy as np
import matplotlib.pyplot as plt

# The goal here is to move from 0D (local cell kinetics) to 1D (wave propagation)
# We want to simulate a cable of cardiac tissue and we want to prove that our numerical solver is correct
# We'ill compare the simulated wave speed with the theoretical formula derived in our research trail

# 1) Parameters and discretisation

# This is the coupling strength between cells. If D is hig, the wave travels fast
D = 1.0       # cm^2/ms

# As seen in 0D this controls how fast u goes 0 -> 1
k = 8.0       # 1/ms

# Length of the tissue 
L = 10.0      # cm

# Time traval
T_max = 40.0  # ms

# Discretization
dx = 0.05             # Space step st
dt = 0.0001           # Time step ms
Nx = int(L / dx) + 1  # Number of spatial points
Nt = int(T_max / dt)  # Number of time steps

x = np.linspace(0, L, Nx)

# CFL
# Since we use an Explicit Euler scheme information cannot travel faster than one grid cell per time step
# If sigma > 0.5 the simulation will explode causing numerical instability.

sigma = D * dt / (dx**2)
if sigma > 0.5:
    raise ValueError("No CLF.")


# 2) The Finite Difference
# We ignore the variable v temporarily to strictly check the theoretical speed formula
        # Reminder : 
        # u[2:]    is u_{i+1}
        # u[1:-1]  is u_{i}
        # u[:-2]   is u_{i-1}

def solve_1d(alpha_val):
    u = np.zeros(Nx)
    u_new = np.zeros(Nx)

    # Stimulation by applying a shock on the left side
    u[0:10] = 1.0 #I_stim

    # We will track the position of the wave front to calculate velocity later
    wave_pos = []
    wave_times = []

    # Time Loop (n)
    # We use our scheme detailled in our research trail
    for n in range(Nt):
        
        laplacian = (u[2:] - 2*u[1:-1] + u[:-2]) / (dx**2)
        
        u_inner = u[1:-1]
        reaction = k * u_inner * (1.0 - u_inner) * (u_inner - alpha_val)
        u_new[1:-1] = u_inner + dt * (D * laplacian + reaction)

        # Left Boundary (i=0)
        lap_left = (2*u[1] - 2*u[0]) / (dx**2)
        reac_left = k * u[0] * (1 - u[0]) * (u[0] - alpha_val)
        u_new[0] = u[0] + dt * (D * lap_left + reac_left)

        # Right Boundary (i=Nx-1)
        lap_right = (2*u[-2] - 2*u[-1]) / (dx**2)
        reac_right = k * u[-1] * (1 - u[-1]) * (u[-1] - alpha_val)
        u_new[-1] = u[-1] + dt * (D * lap_right + reac_right)
        
        u[:] = u_new[:]

        # We don't measure every step because too slow.
        if n % 500 == 0:
            # We define the wavefront location as the point where u crosses 0.5
            # We search from right to left to find the leading edge
            
            idx_front = np.where(u > 0.5)[0]
            if len(idx_front) > 0:
                
                # We take the furthest point to the right
                front_x = x[idx_front[-1]]
                
                # We only record if the wave is traveling
                if 1.0 < front_x < L - 1.0:
                    wave_pos.append(front_x)
                    wave_times.append(n * dt)

    
    # Velocity Calculation via Linear Regression
    # Position = Velocity * Time + Constant
    if len(wave_pos) > 5:
        velocity = np.polyfit(wave_times, wave_pos, 1)[0]
    else:
        velocity = 0.0 # Wave died or didn't propagate
    return velocity

# 3) Sensitivity Analysis
# We loop over different values of alpha (Threshold).

alpha_list = np.linspace(0.0, 0.5, 15)
velocities_sim = []
velocities_theory = []

for val in alpha_list:
    v_sim = solve_1d(val)
    velocities_sim.append(v_sim)

    # The theoretical prediction
    v_th = np.sqrt(2 * D * k) * (0.5 - val)
    velocities_theory.append(v_th)
    

# 4) Plot
plt.figure(figsize=(10, 6))
plt.plot(alpha_list, velocities_theory, 'r-', linewidth=2, label='Analytical Theory')
plt.plot(alpha_list, velocities_sim, 'bo--', markersize=6, label='Finite Difference Simulation')

plt.title(r"Validation : Wave Velocity vs. Excitation Threshold $\alpha$", fontsize=14)
plt.xlabel(r"Threshold $\alpha$", fontsize=12)
plt.ylabel("Conduction Velocity $c$ (cm/ms)", fontsize=12)
plt.legend()
plt.grid(True, alpha=0.3)
plt.axvline(0.5, color='gray', linestyle=':')
plt.text(0.4, 1, "Propagation Failure Region", color='red')

plt.tight_layout()
plt.savefig("chap1_CFL_1D.png", dpi=300)
plt.show()
