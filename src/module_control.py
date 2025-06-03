import numpy as np

def fixed_rudder(t, state):
    return 0.0

def switching_rudder(t, state):
    return 0.1 if t < 50 else -0.1
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

