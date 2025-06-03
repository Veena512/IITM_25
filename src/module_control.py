import numpy as np

def fixed_rudder(t, state):
    delta_c = np.deg2rad(10)  # Fixed rudder angle of 10 degrees
    return delta_c

def switching_rudder(t, state):
    delta_max = np.deg2rad(10)  # maximum rudder angle (10 degrees)
    omega = 0.2  # switching frequency in rad/s

    # Switch direction every half-period based on sine wave
    delta_c = delta_max if np.sin(omega * t) >= 0 else -delta_max

    return delta_c
# Example usage
if __name__ == "__main__":
    t = 5.0  # arbitrary time
    state = np.zeros((13, 1))  # dummy vessel state
    rudder_cmd = fixed_rudder(t, state)
    print(f"Commanded rudder angle δ_c: {rudder_cmd:.4f} rad")
    
    for t in np.linspace(0, 40, 9):
        state = np.zeros((13, 1))
        rudder_cmd = switching_rudder(t, state)
        print(f"t={t:4.1f} s -> δ_c = {rudder_cmd:.4f} rad")

