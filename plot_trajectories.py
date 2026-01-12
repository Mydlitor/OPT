import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import os

# Function definition for plotting contours
def ff4T(x1, x2):
    """Test function: f(x1,x2) = (1/6)x1^6 - 1.05*x1^4 + 2*x1^2 + x2^2 + x1*x2"""
    return (1.0/6.0) * x1**6 - 1.05 * x1**4 + 2.0 * x1**2 + x2**2 + x1 * x2

# Create meshgrid for contour plots
x1_range = np.linspace(-3, 3, 300)
x2_range = np.linspace(-3, 3, 300)
X1, X2 = np.meshgrid(x1_range, x2_range)
Z = ff4T(X1, X2)

# Function to plot trajectory on contour
def plot_trajectory_on_contour(ax, x1_range, x2_range, Z, trajectory_files, labels, title, colors=None):
    """Plot contour with multiple trajectories"""
    # Plot contour
    contour = ax.contour(x1_range, x2_range, Z, levels=30, cmap='viridis', alpha=0.6)
    ax.clabel(contour, inline=True, fontsize=8)
    
    # Default colors if not provided
    if colors is None:
        colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown']
    
    # Plot each trajectory
    for i, (file, label) in enumerate(zip(trajectory_files, labels)):
        if os.path.exists(file):
            try:
                data = pd.read_csv(file)
                if len(data) > 0:
                    ax.plot(data['x1'], data['x2'], 'o-', 
                           color=colors[i % len(colors)], 
                           label=label, 
                           linewidth=2, 
                           markersize=4,
                           alpha=0.8)
                    # Mark start and end points
                    ax.plot(data['x1'].iloc[0], data['x2'].iloc[0], 'o', 
                           color=colors[i % len(colors)], 
                           markersize=10, 
                           markeredgecolor='black', 
                           markeredgewidth=1.5)
                    ax.plot(data['x1'].iloc[-1], data['x2'].iloc[-1], '*', 
                           color=colors[i % len(colors)], 
                           markersize=15, 
                           markeredgecolor='black', 
                           markeredgewidth=1.5)
            except Exception as e:
                print(f"Error loading {file}: {e}")
        else:
            print(f"File not found: {file}")
    
    ax.set_xlabel('x1', fontsize=12)
    ax.set_ylabel('x2', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(loc='best', fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim([-3, 3])
    ax.set_ylim([-3, 3])

# Create figure with 6 subplots
fig = plt.figure(figsize=(18, 12))

# ==============================================================
# PLOT 1: All methods with h=0.05
# ==============================================================
ax1 = plt.subplot(2, 3, 1)
plot_trajectory_on_contour(
    ax1, x1_range, x2_range, Z,
    trajectory_files=[
        'trajectory_SD_h005.csv',
        'trajectory_CG_h005.csv',
        'trajectory_Newton_h005.csv'
    ],
    labels=['SD (h=0.05)', 'CG (h=0.05)', 'Newton (h=0.05)'],
    title='Wykres 1: Wszystkie metody dla h=0.05',
    colors=['red', 'blue', 'green']
)

# ==============================================================
# PLOT 2: All methods with h=0.25
# ==============================================================
ax2 = plt.subplot(2, 3, 2)
plot_trajectory_on_contour(
    ax2, x1_range, x2_range, Z,
    trajectory_files=[
        'trajectory_SD_h025.csv',
        'trajectory_CG_h025.csv',
        'trajectory_Newton_h025.csv'
    ],
    labels=['SD (h=0.25)', 'CG (h=0.25)', 'Newton (h=0.25)'],
    title='Wykres 2: Wszystkie metody dla h=0.25',
    colors=['red', 'blue', 'green']
)

# ==============================================================
# PLOT 3: All methods with line search (variable step)
# ==============================================================
ax3 = plt.subplot(2, 3, 3)
plot_trajectory_on_contour(
    ax3, x1_range, x2_range, Z,
    trajectory_files=[
        'trajectory_SD_linesearch.csv',
        'trajectory_CG_linesearch.csv',
        'trajectory_Newton_linesearch.csv'
    ],
    labels=['SD (line search)', 'CG (line search)', 'Newton (line search)'],
    title='Wykres 3: Wszystkie metody - wersja zmiennokrokowa',
    colors=['red', 'blue', 'green']
)

# ==============================================================
# PLOT 4: Steepest Descent with all step sizes
# ==============================================================
ax4 = plt.subplot(2, 3, 4)
plot_trajectory_on_contour(
    ax4, x1_range, x2_range, Z,
    trajectory_files=[
        'trajectory_SD_h005.csv',
        'trajectory_SD_h025.csv',
        'trajectory_SD_linesearch.csv'
    ],
    labels=['SD (h=0.05)', 'SD (h=0.25)', 'SD (line search)'],
    title='Wykres 4: Metoda najszybszego spadku\ndla różnych kroków',
    colors=['red', 'orange', 'darkred']
)

# ==============================================================
# PLOT 5: Conjugate Gradient with all step sizes
# ==============================================================
ax5 = plt.subplot(2, 3, 5)
plot_trajectory_on_contour(
    ax5, x1_range, x2_range, Z,
    trajectory_files=[
        'trajectory_CG_h005.csv',
        'trajectory_CG_h025.csv',
        'trajectory_CG_linesearch.csv'
    ],
    labels=['CG (h=0.05)', 'CG (h=0.25)', 'CG (line search)'],
    title='Wykres 5: Metoda gradientów sprzężonych\ndla różnych kroków',
    colors=['blue', 'cyan', 'darkblue']
)

# ==============================================================
# PLOT 6: Newton method with all step sizes
# ==============================================================
ax6 = plt.subplot(2, 3, 6)
plot_trajectory_on_contour(
    ax6, x1_range, x2_range, Z,
    trajectory_files=[
        'trajectory_Newton_h005.csv',
        'trajectory_Newton_h025.csv',
        'trajectory_Newton_linesearch.csv'
    ],
    labels=['Newton (h=0.05)', 'Newton (h=0.25)', 'Newton (line search)'],
    title='Wykres 6: Metoda Newtona\ndla różnych kroków',
    colors=['green', 'lime', 'darkgreen']
)

plt.tight_layout()
plt.savefig('trajectories_all_plots.png', dpi=300, bbox_inches='tight')
print("Saved: trajectories_all_plots.png")

# Also save individual plots for better visibility
for i, (ax_num, ax_title) in enumerate([
    (1, 'plot1_methods_h005'),
    (2, 'plot2_methods_h025'),
    (3, 'plot3_methods_linesearch'),
    (4, 'plot4_SD_all_steps'),
    (5, 'plot5_CG_all_steps'),
    (6, 'plot6_Newton_all_steps')
]):
    fig_single = plt.figure(figsize=(10, 8))
    
    # Recreate the specific plot
    if i == 0:  # Plot 1
        plot_trajectory_on_contour(
            plt.gca(), x1_range, x2_range, Z,
            trajectory_files=['trajectory_SD_h005.csv', 'trajectory_CG_h005.csv', 'trajectory_Newton_h005.csv'],
            labels=['SD (h=0.05)', 'CG (h=0.05)', 'Newton (h=0.05)'],
            title='Wykres 1: Wszystkie metody dla h=0.05',
            colors=['red', 'blue', 'green']
        )
    elif i == 1:  # Plot 2
        plot_trajectory_on_contour(
            plt.gca(), x1_range, x2_range, Z,
            trajectory_files=['trajectory_SD_h025.csv', 'trajectory_CG_h025.csv', 'trajectory_Newton_h025.csv'],
            labels=['SD (h=0.25)', 'CG (h=0.25)', 'Newton (h=0.25)'],
            title='Wykres 2: Wszystkie metody dla h=0.25',
            colors=['red', 'blue', 'green']
        )
    elif i == 2:  # Plot 3
        plot_trajectory_on_contour(
            plt.gca(), x1_range, x2_range, Z,
            trajectory_files=['trajectory_SD_linesearch.csv', 'trajectory_CG_linesearch.csv', 'trajectory_Newton_linesearch.csv'],
            labels=['SD (line search)', 'CG (line search)', 'Newton (line search)'],
            title='Wykres 3: Wszystkie metody - wersja zmiennokrokowa',
            colors=['red', 'blue', 'green']
        )
    elif i == 3:  # Plot 4
        plot_trajectory_on_contour(
            plt.gca(), x1_range, x2_range, Z,
            trajectory_files=['trajectory_SD_h005.csv', 'trajectory_SD_h025.csv', 'trajectory_SD_linesearch.csv'],
            labels=['SD (h=0.05)', 'SD (h=0.25)', 'SD (line search)'],
            title='Wykres 4: Metoda najszybszego spadku dla różnych kroków',
            colors=['red', 'orange', 'darkred']
        )
    elif i == 4:  # Plot 5
        plot_trajectory_on_contour(
            plt.gca(), x1_range, x2_range, Z,
            trajectory_files=['trajectory_CG_h005.csv', 'trajectory_CG_h025.csv', 'trajectory_CG_linesearch.csv'],
            labels=['CG (h=0.05)', 'CG (h=0.25)', 'CG (line search)'],
            title='Wykres 5: Metoda gradientów sprzężonych dla różnych kroków',
            colors=['blue', 'cyan', 'darkblue']
        )
    elif i == 5:  # Plot 6
        plot_trajectory_on_contour(
            plt.gca(), x1_range, x2_range, Z,
            trajectory_files=['trajectory_Newton_h005.csv', 'trajectory_Newton_h025.csv', 'trajectory_Newton_linesearch.csv'],
            labels=['Newton (h=0.05)', 'Newton (h=0.25)', 'Newton (line search)'],
            title='Wykres 6: Metoda Newtona dla różnych kroków',
            colors=['green', 'lime', 'darkgreen']
        )
    
    plt.savefig(f'{ax_title}.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {ax_title}.png")
    plt.close(fig_single)

plt.show()

print("\nAll plots generated successfully!")
print("Generated files:")
print("  - trajectories_all_plots.png (all 6 plots combined)")
print("  - plot1_methods_h005.png")
print("  - plot2_methods_h025.png")
print("  - plot3_methods_linesearch.png")
print("  - plot4_SD_all_steps.png")
print("  - plot5_CG_all_steps.png")
print("  - plot6_Newton_all_steps.png")
