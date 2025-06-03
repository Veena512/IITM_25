import numpy as np

def Smat(vec):
    x, y, z = vec.flatten()
    return np.array([
        [0,   -z,  y],
        [z,    0, -x],
        [-y,   x,  0]
    ])

def eul_to_rotm(eul):
    roll, pitch, yaw = eul.flatten()

    # Rotation about X-axis (roll)
    Rx = np.array([
        [1, 0, 0],
        [0, np.cos(roll), -np.sin(roll)],
        [0, np.sin(roll),  np.cos(roll)]
    ])

    # Rotation about Y-axis (pitch)
    Ry = np.array([
        [ np.cos(pitch), 0, np.sin(pitch)],
        [0, 1, 0],
        [-np.sin(pitch), 0, np.cos(pitch)]
    ])

    # Rotation about Z-axis (yaw)
    Rz = np.array([
        [np.cos(yaw), -np.sin(yaw), 0],
        [np.sin(yaw),  np.cos(yaw), 0],
        [0, 0, 1]
    ])

    # Final rotation matrix: R = Rz * Ry * Rx
    R = Rz @ Ry @ Rx
    return R
    
def eul_rate_matrix(eul):
    roll, pitch, yaw = eul.flatten()

    c_phi = np.cos(roll)
    s_phi = np.sin(roll)
    c_theta = np.cos(pitch)
    s_theta = np.sin(pitch)

    J2 = np.array([
        [1, 0, -s_theta],
        [0, c_phi, s_phi * c_theta],
        [0, -s_phi, c_phi * c_theta]
    ])

    return J2
# Example usage
if __name__ == "__main__":
    # Input vector (3x1)
    v = np.array([[1], [2], [3]])
    
    # Create skew-symmetric matrix
    S = Smat(v)
    print("Skew-symmetric matrix (S):")
    print(S)
    
    # Multiply with another vector
    u = np.array([[4], [5], [6]])
    cross_product = S @ u
    print("\nCross product result (v × u):")
    print(cross_product)
    
    eul = np.array([[np.pi/6], [np.pi/4], [np.pi/3]])  # 30°, 45°, 60°

    # Get rotation matrix
    R = eul_to_rotm(eul)
    print("Rotation matrix:")
    print(R)

    # Example vector in body coordinates
    v_body = np.array([[1], [0], [0]])
    v_global = R @ v_body
    print("\nTransformed vector in global coordinates:")
    print(v_global)
    
    eul = np.array([[np.pi/6], [np.pi/4], [np.pi/3]])  # roll, pitch, yaw in radians
    J2 = eul_rate_matrix(eul)

    print("Euler Rate Matrix (J2):")
    print(J2)

    # Example angular velocity in body frame
    omega = np.array([[0.5], [0.2], [0.1]])  # rad/s
    eul_dot = np.linalg.inv(J2) @ omega     # Euler angle rates
    print("\nEuler angle rates (rad/s):")
    print(eul_dot)

