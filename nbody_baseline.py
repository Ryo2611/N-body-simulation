import pandas as pd
import matplotlib.pyplot as plt

# Load the C++ output
df = pd.read_csv('nbody_output.csv')

plt.figure(figsize=(8, 8))
plt.plot(df['b1_x'], df['b1_y'], label='Body 1')
plt.plot(df['b2_x'], df['b2_y'], label='Body 2')
plt.plot(df['b3_x'], df['b3_y'], label='Body 3')
plt.legend()
plt.title("N-Body Simulation Trajectories")
plt.xlabel("X")
plt.ylabel("Y")
plt.grid(True)
plt.show()