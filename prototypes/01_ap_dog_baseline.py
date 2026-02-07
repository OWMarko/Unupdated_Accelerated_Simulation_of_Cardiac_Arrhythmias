import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# The goal is to understand the physical phenomen and biological behind these equations. First we look at nullclines
# which are fixed points in a way. We must define parameters, our epsilon, Aliev Panfilov, our roots and all kind of maths 
# explained in our recherch trail

# 1)
# These parameters are phenomenological, not constants. These parameters are chosen by Aliev and Panfilov to force the mathematic equation
# to match experimental data recorded from dog ventricles. Each parameter has his physical meaning.


# Maximum Conductance of the fast sodium channels, these value controls how aggressive the depolarization (transition from  negative value to a positive)
# Mathematicly it scales the reaction term relative to the diffusion term. Thanks to k = 8.0 the shape of the wave looks like a real heartbeat
k = 8.0

# Excitation treshold, our value means "push the system 15% of the way to the max voltage to trigger a self-sustaining wave"
# We must make sure that alpha isn't too high because if alpha > 1/2 the tissue becomes hard to excite and if alpha is near 0
# the heart becomes unstable and fires randomly
alpha = 0.15 

# The Stiffness/Restitution Parameters, it's the tuned part of the model done by Aliev and Panfilov
# This sets the maximum length of the action potential, a smaller parameter means a longer plateau at beging
epsilon_0 = 0.002

# These parameters were fitted to make the Action Potential Duration Restituation Curve match with the experimental data of a dog heart
mu_1 = 0.2
mu_2 = 0.3


# It makes the value of epsilon change depending on the current voltage u and the memory of the previous beat v.
# This ensures that if you stimulate the cell before it has fully recovered the next beat will be shorter.
def epsilon(u, v):
    return epsilon_0 + (mu_1 * v) / (u + mu_2)


# 2)
# Now we can code our 0D equation. It's simple in python, there isn't any difficulty. Remark : we wrote the equation without the - 

def aliev_panfilov(y, t, k, alpha, eps0, mu1, mu2):
    u, v = y

    epsilon = eps0 + (mu1 * v) / (u + mu2)
    
    dudt = k * u * (1 - u) * (u - alpha) - u * v

    dvdt = epsilon * (-v - k * u * (u - alpha - 1))

    return [dudt, dvdt]

#3)
# Now we can do our analyse (phase plan), we create a grid of points wich is our playground

u_range = np.linspace(-0.2, 1.2, 25)
v_range = np.linspace(-0.1, 0.5, 25)
U, V = np.meshgrid(u_range, v_range)

# And we ensure that we calculate epsilon for every point on the grid
EPS_Grid = epsilon_0 + (mu_1 * V) / (U + mu_2)

# Compute derivatives at each point for the vector field 
dU = k * U * (1 - U) * (U - alpha) - U * V
dV = EPS_Grid * (-V - k * U * (U - alpha - 1))

# And we need to normalize arrows to see the physics and the direction (and we avoid division by 0)
M = np.hypot(dU, dV)
M[M == 0] = 1.
    
dU_norm = dU / M
dV_norm = dV / M

# 4) 
# Now we can see create our nullclines
u_vals = np.linspace(-0.1, 1.1, 400)

# u-nullcline 1
v_null_u = k * (1 - u_vals) * (u_vals - alpha)

# u-nullcline 2 since it's juste the u = 0 line we'ill plot it later

# v-nullcline 1
v_null_v = k * u_vals * (1 + alpha - u_vals)

# 5)
# Don't forget to simulate our equation and we start alpha above 0.15 to capture the excitation and to trigger it
t = np.linspace(0, 100, 2000) # Increased time to see full cycle
y0 = [0.16, 0.0] 
sol = odeint(aliev_panfilov, y0, t, args=(k, alpha, epsilon_0, mu_1, mu_2))

# 6)
# We plot
plt.figure(figsize=(12, 10))
plt.title("Phase Plane : Aliev-Panfilov", fontsize=16)


plt.quiver(U, V, dU_norm, dV_norm, color='gray', alpha=0.4, pivot='mid', scale=25)


plt.plot(u_vals, v_null_u, 'r-', linewidth=2, label=r'$u$-nullcline (Parabola)')

plt.axvline(x=0, color='r', linewidth=2, linestyle='-', alpha=0.6, label=r'$u$-nullcline ($u=0$)')


plt.plot(u_vals, v_null_v, 'b--', linewidth=2, label=r'$v$-nullcline ($\dot{v}=0$)')
plt.plot(0, 0, 'ko', markersize=10, markerfacecolor='black', label='Stable Fixed Point (Rest)')
plt.plot(sol[:, 0], sol[:, 1], 'k-', linewidth=3, label='Limit Cycle')

plt.text(0.5, 0.4, 'Vector Field normalized', color='gray', fontsize=10, style='italic')
plt.xlabel('Voltage $u$ (dimensionless)', fontsize=14)
plt.ylabel('Recovery $v$ (dimensionless)', fontsize=14)
plt.ylim(-0.15,2.7) 
plt.xlim(-0.15, 1.1)
plt.legend(fontsize=12, loc='upper right')
plt.grid(True, linestyle=':', alpha=0.6)

plt.show()
