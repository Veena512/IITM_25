import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

class Vessel:
    def __init__(self):
        # Any vessel-specific parameters can be initialized here
        pass

    def vessel_ode(self, t, state):
        """
        Returns the time derivative of the state vector.
        State is a 13x1 vector:
        [u, v, w, p, q, r, x, y, z, phi, theta, psi, delta]
        """
        # Unpack state variables
        u, v, w = state[0], state[1], state[2]       # Linear velocities
        p, q, r = state[3], state[4], state[5]       # Angular velocities
        x, y, z = state[6], state[7], state[8]       # Position
        phi, theta, psi = state[9], state[10], state[11]  # Orientation
        delta = state[12]                            # Rudder angle

        # For now, use a simple constant model (placeholder dynamics)
        # You should replace this with your vessel's dynamics
        du = 0
        dv = 0
        dw = 0
        dp = 0
        dq = 0
        dr = 0

        dx = u * np.cos(psi) - v * np.sin(psi)
        dy = u * np.sin(psi) + v * np.cos(psi)
        dz = 0  # Assume no vertical motion

        dphi = p
        dtheta = q
        dpsi = r

        ddelta = 0  # Assuming constant rudder for now

        # Assemble derivative vector
        dstate = np.array([
            du, dv, dw, dp, dq, dr,
            dx, dy, dz, dphi, dtheta, dpsi,
            ddelta
        ]).reshape((13, 1))

        return dstate

    def simulate(self):
        # Initial state vector [u, v, w, p, q, r, x, y, z, phi, theta, psi, delta]
        state0 = np.zeros((13, 1))
        state0[0] = 2.0  # Initial forward speed (u)
        state0[12] = 0.1  # Rudder angle (delta)

        state0 = state0.flatten()

        t0 = 0.0
        tf = 100.0
        t_eval = np.linspace(t0, tf, 1000)

        # Define a wrapper for solve_ivp
        def ode_wrapper(t, state):
            return self.vessel_ode(t, state).flatten()

        # Run the simulation
        sol = solve_ivp(ode_wrapper, (t0, tf), state0, t_eval=t_eval, method='RK45')

        # Extract position for plotting
        x = sol.y[6]
        y = sol.y[7]

        # Plot trajectory
        plt.figure(figsize=(8, 5))
        plt.plot(x, y, label='Vessel Trajectory')
        plt.xlabel("X Position (m)")
        plt.ylabel("Y Position (m)")
        plt.title("Vessel Trajectory Simulation")
        plt.axis("equal")
        plt.grid(True)
        plt.legend()
        plt.show()


# --- Run simulation ---
if __name__ == "__main__":
    vessel = Vessel()
    vessel.simulate()

