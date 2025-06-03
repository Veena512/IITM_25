import numpy as np

class Vessel:
    def __init__(self, m, Iz, X_u_dot, Y_v_dot, N_r_dot):
        self.m = m
        self.Iz = Iz
        self.X_u_dot = X_u_dot
        self.Y_v_dot = Y_v_dot
        self.N_r_dot = N_r_dot
        self.M = None  # to be set in _generate_mass_matrix

    def _generate_mass_matrix(self):
        # Rigid-body mass matrix (3 DOF: surge, sway, yaw)
        M_RB = np.array([
            [self.m,     0,      0],
            [0,     self.m,      0],
            [0,         0,   self.Iz]
        ])

        # Added mass matrix (symmetric, negative due to added inertia)
        M_A = -np.array([
            [self.X_u_dot,     0,          0],
            [0,           self.Y_v_dot,    0],
            [0,               0,      self.N_r_dot]
        ])

        # Total mass matrix
        self.M = M_RB + M_A

    def __init__(self):
        # Example non-dimensional coefficients
        self.C_Xu = -0.05
        self.C_Yv = -0.1
        self.C_Nr = -0.2
        self.C_Xu_dot = -0.01
        self.C_Yv_dot = -0.02
        self.C_Nr_dot = -0.03

        # Dimensional coefficients (to be computed)
        self.X_u = 0
        self.Y_v = 0
        self.N_r = 0
        self.X_u_dot = 0
        self.Y_v_dot = 0
        self.N_r_dot = 0

    def _dimensionalize_coefficients(self, rho, L, U):
        # Precompute constants
        q = 0.5 * rho * U**2
        qL = q * L
        qL2 = q * L**2

        # Linear damping derivatives
        self.X_u = qL * self.C_Xu
        self.Y_v = qL * self.C_Yv
        self.N_r = qL2 * self.C_Nr

        # Added mass derivatives
        self.X_u_dot = rho * L**2 * self.C_Xu_dot
        self.Y_v_dot = rho * L**2 * self.C_Yv_dot
        self.N_r_dot = rho * L**4 * self.C_Nr_dot
# Example usage:
if __name__ == "__main__":
    vessel = Vessel(
        m=500.0,
        Iz=2000.0,
        X_u_dot=-50.0,
        Y_v_dot=-200.0,
        N_r_dot=-100.0
    )
    vessel._generate_mass_matrix()
    print("Mass matrix M:\n", vessel.M)

