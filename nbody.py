import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_trajectory(filename='output.csv'):
    # 1. Load the data
    # Pandas is the standard tool for time-series data in Finance (OHLC data, etc.)
    try:
        df = pd.read_csv(filename)
    except FileNotFoundError:
        print(f"Error: Could not find {filename}. Run the C++ code first!")
        return

    # 2. Setup the 3D Plot
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # 3. Plot the trajectory
    # We use the 't' column for color to show time evolution (Start=Blue, End=Red)
    scatter = ax.scatter(df['x'], df['y'], df['z'], 
                         c=df['t'], cmap='coolwarm', s=10, alpha=0.6)

    # Add a line to connect the points for clarity
    ax.plot(df['x'], df['y'], df['z'], color='gray', alpha=0.3, linewidth=1)

    # 4. Styling (Crucial for presentation)
    ax.set_title(f'Particle Trajectory (ID: {df["id"].iloc[0]})', fontsize=14)
    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    ax.set_zlabel('Z Position')
    
    # Add a color bar to indicate time
    cbar = plt.colorbar(scatter, ax=ax, pad=0.1)
    cbar.set_label('Time')

    # 5. Show and Save
    plt.tight_layout()
    plt.savefig('trajectory_plot.png', dpi=300) # Save high-quality image
    print("Plot saved to 'trajectory_plot.png'")
    plt.show()

if __name__ == "__main__":
    plot_trajectory()