import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Now we want to show that if we change some parameters and we change de physiology of our equation we encounter tresholds
# As we are in prototypes stage, we don't do complexe algos, imports or what else, we do basic coding to test our understanding and intuitions

def aliev_panfilov(y, t, k, alpha, eps0, mu1, mu2):
    u, v = y

    epsilon = eps0 + (mu1 * v) / (u + mu2)
    
    dudt = k * u * (1 - u) * (u - alpha) - u * v

    dvdt = epsilon * (-v - k * u * (u - alpha - 1))

    return [dudt, dvdt]


t = np.linspace(0, 600, 2000)
y0_stim = [0.16, 0.0]   

params_dog = {
    'k': 8.0, 
    'alpha': 0.15, 
    'eps0': 0.002, 
    'mu1': 0.2, 
    'mu2': 0.3
}

# Our first intuition is to see thge sensitivity of alpha
# Hypothesis : As alpha increases the system becomes harder to excite
# If alpha exceeds our stimulus (0.16) the wave should die

a_values = np.linspace(0.05, 0.25, 50)
peak_voltages = []

for val in a_values:
    sol = odeint(aliev_panfilov_dynamic, y0_stim, t, 
                 args=(params_dog['k'], val, params_dog['eps0'], params_dog['mu1'], params_dog['mu2']))
    
    peak_voltages.append(np.max(sol[:, 0]))

plt.figure(figsize=(10, 6))
plt.plot(a_values, peak_voltages, 'o-', color='darkred', linewidth=2, markersize=4)
plt.axvline(x=0.15, color='black', linestyle='--', label='Selected $a=0.15$')
plt.axhline(y=0.8, color='grey', linestyle=':', label='Excitation Threshold')
plt.title(r"Sensitivity of Excitation to Parameter $alpha$", fontsize=14)
plt.xlabel(r"Parameter $alpha$ (Threshold)", fontsize=12)
plt.ylabel("Peak Voltage ($u_{max}$)", fontsize=12)
plt.text(0.06, 0.9, "Region : Spontaneous/Easy Fire", color='green', fontsize=10)
plt.text(0.17, 0.2, "Region : Wave Death", color='red', fontsize=10)
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()


# The second intuition is, what happens when k is high ? In Mathematical Physiology II: Systems Physiology
# we can read at page 537 that Na+ Ion concentrations in most cardiac cells is at : 60 mV near potential, 15 nM intracellular, 145 nM extracellular
# Hypothesis : High k => Fast Na+ channel opening i.e Vertical upstroke

# High K = 8 (Healthy dog intracellular concentration)
sol_k_high = odeint(aliev_panfilov_dynamic, y0_stim, t, 
                    args=(8.0, 0.15, 0.002, 0.2, 0.3))

# Low K = 2 (Pathological dog intracellular concentration)
sol_k_low = odeint(aliev_panfilov_dynamic, y0_stim, t, 
                   args=(2.0, 0.15, 0.002, 0.2, 0.3))

plt.figure(figsize=(10, 6))
plt.plot(t, sol_k_high[:, 0], 'b-', linewidth=2.5, label=r'$k=8.0$ (Healthy Sodium Current)')
plt.plot(t, sol_k_low[:, 0], 'r--', linewidth=2.5, label=r'$k=2.0$ (Slow/Pathological)')
plt.xlim(0, 50) 
plt.title(r"Effect of $k$ on Depolarization Velocity (Upstroke)", fontsize=14)
plt.xlabel("Time (dimensionless)", fontsize=12)
plt.ylabel("Voltage $u$", fontsize=12)
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()


# Since the begining we focused the dog case but what about Human ? Is it the same ?
# Hypothesis: Lowering epsilon stretches the plateau creating a "stiff" system
# According the biologist and meds, we lower eps, more details in report (spoiler : slower recovery)

# Dog 
sol_dog = sol_k_high 

# Human
sol_human = odeint(aliev_panfilov_dynamic, y0_stim, t, 
                   args=(8.0, 0.15, 0.0005, 0.2, 0.3))

plt.figure(figsize=(10, 6))
plt.plot(t, sol_dog[:, 0], 'k--', linewidth=2, label='Dog Model (Baseline)')
plt.plot(t, sol_human[:, 0], 'b-', linewidth=2.5, label='Human Model (Stiff)')

plt.arrow(100, 0.5, 200, 0, head_width=0.05, head_length=10, fc='k', ec='k')
plt.text(120, 0.55, "Extended Plateau = High Stiffness", fontsize=11, fontweight='bold')

plt.title(r"The Stiffness Problem", fontsize=14)
plt.xlabel("Time (dimensionless)", fontsize=12)
plt.ylabel("Voltage $u$", fontsize=12)
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
