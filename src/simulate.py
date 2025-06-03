import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# -----------------------------
# Define the ship dynamics
# -----------------------------

def ship_dynamics(t, y):
    # State variables
    u, r, psi, x, y_pos, delta = y

    # Ship parameters
    m = 500.0       # mass [kg]
    Iz = 1000.0     # moment of inertia [kg m^2]
    X_u = -50.0     # surge damping [Ns/m]
    N_r = -30.0     # yaw damping [Nms/rad]
    L = 5.0         # rudder moment arm [m]
    delta_max = np.deg2rad(30)  # rudder angle limit
    K_rudder = 10.0             # rudder gain

    # Desired rudder input: step input of 10 degrees
    delta_input = np.deg2rad(10)

    # Rudder angle dynamics (1st-order lag)
    delta_dot = (delta_input - delta) * 1.0  # time constant = 1 s

    # Equations of motion
    du = (X_u * u) / m
    dr = (N_r * r + L * delta) / Iz
    dpsi = r
    dx = u * np.cos(psi)
    dy = u * np.sin(psi)

    return [du, dr, dpsi, dx, dy, delta_dot]

# -----------------------------
# Initial conditions
# -----------------------------

u0 = 5.0       # Initial forward velocity [m/s]
r0 = 0.0       # Initial yaw rate [rad/s]
psi0 = 0.0     # Initial heading angle [rad]
x0 = 0.0       # Initial x position [m]
y0 = 0.0       # Initial y position [m]
delta0 = 0.0   # Initial rudder angle [rad]

y_init = [u0, r0, psi0, x0, y0, delta0]

# Time span for simulation
t_span = (0, 100)  # simulate from t=0 to t=100s
t_eval = np.linspace(t_span[0], t_span[1], 1000)

# -----------------------------
# Solve the system
# -----------------------------

solution = solve_ivp(ship_dynamics, t_span, y_init, t_eval=t_eval, rtol=1e-8, atol=1e-8)

# Extract solution variables
u = solution.y[0]
r = solution.y[1]
psi = solution.y[2]
x = solution.y[3]
y_pos = solution.y[4]
delta = solution.y[5]
time = solution.t

# -----------------------------
# Plotting the results
# -----------------------------

plt.figure(figsize=(12, 10))

# 1. Ship trajectory in x-y plane
plt.subplot(3, 2, 1)
plt.plot(x, y_pos, label='Trajectory')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title('Ship Trajectory (x-y)')
plt.axis('equal')
plt.grid(True)

# 2. Velocity vs time
plt.subplot(3, 2, 2)
plt.plot(time, u, label='Velocity u')
plt.xlabel('Time [s]')
plt.ylabel('Velocity [m/s]')
plt.title('Forward Velocity vs Time')
plt.grid(True)

# 3. Turn rate vs time
plt.subplot(3, 2, 3)
plt.plot(time, r, label='Turn rate r')
plt.xlabel('Time [s]')
plt.ylabel('Yaw Rate [rad/s]')
plt.title('Yaw Rate vs Time')
plt.grid(True)

# 4. Rudder angle vs time
plt.subplot(3, 2, 4)
plt.plot(time, np.rad2deg(delta), label='Rudder angle')
plt.xlabel('Time [s]')
plt.ylabel('Rudder Angle [deg]')
plt.title('Rudder Angle vs Time')
plt.grid(True)

# 5. Yaw angle vs time
plt.subplot(3, 2, 5)
plt.plot(time, np.rad2deg(psi), label='Yaw angle')
plt.xlabel('Time [s]')
plt.ylabel('Yaw Angle [deg]')
plt.title('Yaw Angle vs Time')
plt.grid(True)

plt.tight_layout()
plt.savefig("task9_ship_simulation.png")  # Saves the figure
# plt.show()  # Uncomment if you're using a graphical environment

print("✅ Task 9 simulation completed. Output saved as 'task9_ship_simulation.png'.")

