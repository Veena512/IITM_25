import numpy as np

class Vessel:
    def __init__(self):
        # Non-dimensional hydrodynamic coefficients
        self.C_Xu = -0.1
        self.C_Yv = -0.2
        self.C_Nr = -0.3
        self.C_Xu_dot = -0.05
        self.C_Yv_dot = -0.06
        self.C_Nr_dot = -0.07
        
        # Placeholders for dimensional forms
        self.X_u = 0
        self.Y_v = 0
        self.N_r = 0
        self.X_u_dot = 0
        self.Y_v_dot = 0
        self.N_r_dot = 0
        
        # Mass matrix placeholder
        self.M = np.eye(6)

    # ----------------------
    # Task 6: Generate Mass Matrix
    # ----------------------
    def _generate_mass_matrix(self):
        m = 500.0       # Mass in kg
        Izz = 1000.0    # Yaw inertia
        
        # Simplified 3DOF mass matrix (surge, sway, yaw)
        self.M = np.array([
            [m - self.X_u_dot, 0, 0],
            [0, m - self.Y_v_dot, 0],
            [0, 0, Izz - self.N_r_dot]
        ])

    # ----------------------
    # Task 7: Dimensionalize Coefficients
    # ----------------------
    def _dimensionalize_coefficients(self, rho, L, U):
        q = 0.5 * rho * U**2
        self.X_u = q * L * self.C_Xu
        self.Y_v = q * L * self.C_Yv
        self.N_r = q * L**2 * self.C_Nr
        self.X_u_dot = rho * L**2 * self.C_Xu_dot
        self.Y_v_dot = rho * L**2 * self.C_Yv_dot
        self.N_r_dot = rho * L**4 * self.C_Nr_dot

    # ----------------------
    # Task 8: Vessel ODE
    # ----------------------
    def vessel_ode(self, t, state):
        u, v, w, p, q, r, x, y, z, phi, theta, psi, delta = state.flatten()
        
        # Forces (placeholder)
        X = self.X_u * u
        Y = self.Y_v * v
        N = self.N_r * r
        
        tau = np.array([X, Y, N])
        
        # Accelerations
        nu_dot = np.linalg.inv(self.M) @ tau
        
        # Velocities (placeholder for Euler angle rates)
        eul = np.array([[phi], [theta], [psi]])
        J = eul_rate_matrix(eul)
        eul_dot = J @ np.array([[p], [q], [r]])
        
        # Construct derivative of state vector
        state_dot = np.zeros((13, 1))
        state_dot[0:3] = nu_dot.reshape((3, 1))        # Linear acceleration
        state_dot[3:6] = np.zeros((3, 1))              # Angular acceleration (zero for now)
        state_dot[6:9] = np.array([[u], [v], [w]])     # Position change
        state_dot[9:12] = eul_dot                      # Euler angle change
        state_dot[12] = 0                              # Rudder angle derivative
        
        return state_dot

