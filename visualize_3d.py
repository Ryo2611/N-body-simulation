import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

def main():
    print("Loading simulation output (output.csv)...")
    try:
        df = pd.read_csv('output.csv')
    except Exception as e:
        print(f"Error reading output.csv. Please run './nbody_sim --visualize' first. {e}")
        return

    steps = sorted(df['step'].unique())
    print(f"Loaded {len(steps)} steps.")

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Add a cool dark background theme for the portfolio
    plt.style.use('dark_background')
    fig.patch.set_facecolor('black')
    ax.set_facecolor('black')

    # Remove axes to make it look like space
    ax.axis('off')
    ax.grid(False)

    # Calculate global max/min to fix camera perspective
    max_val = max(df['x'].max(), df['y'].max(), df['z'].max())
    min_val = min(df['x'].min(), df['y'].min(), df['z'].min())
    
    # Add a bit of padding
    padding = (max_val - min_val) * 0.1
    ax.set_xlim(min_val - padding, max_val + padding)
    ax.set_ylim(min_val - padding, max_val + padding)
    ax.set_zlim(min_val - padding, max_val + padding)

    # Use a vibrant color for particles
    scatter = ax.scatter([], [], [], s=2, c='cyan', alpha=0.8, edgecolors='none')

    def update(frame):
        current_df = df[df['step'] == frame]
        scatter._offsets3d = (current_df['x'], current_df['y'], current_df['z'])
        ax.set_title(f"N-Body 3D Octree Simulation\nStep: {frame} / {steps[-1]}", color='white', pad=20)
        return scatter,

    print("Generating animation...")
    ani = animation.FuncAnimation(fig, update, frames=steps, blit=False, interval=30)

    out_file = "simulation_3d.gif"
    print(f"Saving animation to {out_file}. This might take a minute...")
    try:
        ani.save(out_file, writer='pillow', fps=30)
        print(f"Successfully saved {out_file}!")
    except Exception as e:
        print(f"Error saving animation. Saving as mp4 instead... {e}")
        ani.save("simulation_3d.mp4", fps=30)
        print("Successfully saved simulation_3d.mp4!")

if __name__ == '__main__':
    main()
