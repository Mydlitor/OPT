#!/usr/bin/env python3
"""
Script to generate 6 contour plots showing optimization paths for different methods and step sizes.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import cm

# Test function ff4T
def ff4T(x1, x2):
    """Test function: (1/6)*x1^6 - 1.05*x1^4 + 2*x1^2 + x2^2 + x1*x2"""
    return (1.0/6.0) * x1**6 - 1.05 * x1**4 + 2.0 * x1**2 + x2**2 + x1 * x2

# Create grid for contour plots
x1 = np.linspace(-2, 2, 300)
x2 = np.linspace(-2, 2, 300)
X1, X2 = np.meshgrid(x1, x2)
Z = ff4T(X1, X2)

# Color scheme for methods
colors = {'SD': 'red', 'CG': 'blue', 'Newton': 'green'}
labels = {
    'SD': 'Metoda najszybszego spadku', 
    'CG': 'Metoda gradientów sprzężonych', 
    'Newton': 'Metoda Newtona'
}

def load_history(filename):
    """Load iteration history from CSV file, filtering out NaN and inf values"""
    try:
        df = pd.read_csv(filename)
        hist = df[['x1', 'x2']].values
        
        # Filter out rows with NaN or inf values
        valid_mask = np.all(np.isfinite(hist), axis=1)
        hist_filtered = hist[valid_mask]
        
        # If we filtered out some points, report it
        if len(hist_filtered) < len(hist):
            print(f"  Note: {filename} - filtered {len(hist) - len(hist_filtered)} invalid points, kept {len(hist_filtered)} valid points")
        
        # Return None if no valid points remain
        if len(hist_filtered) == 0:
            print(f"Warning: {filename} contains no valid points")
            return None
            
        return hist_filtered
    except FileNotFoundError:
        print(f"Warning: File {filename} not found")
        return None
    except Exception as e:
        print(f"Warning: Error loading {filename}: {e}")
        return None

# WYKRES 1: h=0.05, wszystkie metody
print("Generating Plot 1: h=0.05, all methods...")
fig, ax = plt.subplots(figsize=(10, 8))
contour = ax.contour(X1, X2, Z, levels=30, cmap='viridis', alpha=0.6)
ax.clabel(contour, inline=True, fontsize=8, fmt='%.1f')

for method in ['SD', 'CG', 'Newton']:
    hist = load_history(f'history_{method}_0.05.csv')
    if hist is not None:
        ax.plot(hist[:, 0], hist[:, 1], 'o-', color=colors[method], 
                label=labels[method], markersize=4, linewidth=1.5, alpha=0.8)
        # Mark start and end points
        ax.plot(hist[0, 0], hist[0, 1], 'o', color=colors[method], 
                markersize=8, markeredgewidth=2, markeredgecolor='black')
        ax.plot(hist[-1, 0], hist[-1, 1], '*', color=colors[method], 
                markersize=15, markeredgewidth=1, markeredgecolor='black')

ax.set_xlabel('x1', fontsize=12)
ax.set_ylabel('x2', fontsize=12)
ax.set_title('Wykres 1: Optymalizacja z krokiem h=0.05', fontsize=14, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
plt.savefig('wykres1_h005.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: wykres1_h005.png")

# WYKRES 2: h=0.25, wszystkie metody
print("Generating Plot 2: h=0.25, all methods...")
fig, ax = plt.subplots(figsize=(10, 8))
contour = ax.contour(X1, X2, Z, levels=30, cmap='viridis', alpha=0.6)
ax.clabel(contour, inline=True, fontsize=8, fmt='%.1f')

for method in ['SD', 'CG', 'Newton']:
    hist = load_history(f'history_{method}_0.25.csv')
    if hist is not None:
        ax.plot(hist[:, 0], hist[:, 1], 'o-', color=colors[method], 
                label=labels[method], markersize=4, linewidth=1.5, alpha=0.8)
        # Mark start and end points
        ax.plot(hist[0, 0], hist[0, 1], 'o', color=colors[method], 
                markersize=8, markeredgewidth=2, markeredgecolor='black')
        ax.plot(hist[-1, 0], hist[-1, 1], '*', color=colors[method], 
                markersize=15, markeredgewidth=1, markeredgecolor='black')

ax.set_xlabel('x1', fontsize=12)
ax.set_ylabel('x2', fontsize=12)
ax.set_title('Wykres 2: Optymalizacja z krokiem h=0.25', fontsize=14, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
plt.savefig('wykres2_h025.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: wykres2_h025.png")

# WYKRES 3: wersja zmiennokrokowa (h=0), wszystkie metody
print("Generating Plot 3: variable step, all methods...")
fig, ax = plt.subplots(figsize=(10, 8))
contour = ax.contour(X1, X2, Z, levels=30, cmap='viridis', alpha=0.6)
ax.clabel(contour, inline=True, fontsize=8, fmt='%.1f')

for method in ['SD', 'CG', 'Newton']:
    hist = load_history(f'history_{method}_variable.csv')
    if hist is not None:
        ax.plot(hist[:, 0], hist[:, 1], 'o-', color=colors[method], 
                label=labels[method], markersize=4, linewidth=1.5, alpha=0.8)
        # Mark start and end points
        ax.plot(hist[0, 0], hist[0, 1], 'o', color=colors[method], 
                markersize=8, markeredgewidth=2, markeredgecolor='black')
        ax.plot(hist[-1, 0], hist[-1, 1], '*', color=colors[method], 
                markersize=15, markeredgewidth=1, markeredgecolor='black')

ax.set_xlabel('x1', fontsize=12)
ax.set_ylabel('x2', fontsize=12)
ax.set_title('Wykres 3: Optymalizacja ze zmiennym krokiem (line search)', fontsize=14, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
plt.savefig('wykres3_variable.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: wykres3_variable.png")

# WYKRES 4: Metoda SD, wszystkie kroki
print("Generating Plot 4: SD method, all step sizes...")
fig, ax = plt.subplots(figsize=(10, 8))
contour = ax.contour(X1, X2, Z, levels=30, cmap='viridis', alpha=0.6)
ax.clabel(contour, inline=True, fontsize=8, fmt='%.1f')

step_configs = [('0.05', 'h=0.05'), ('0.25', 'h=0.25'), ('variable', 'zmienny krok')]
line_styles = ['-', '--', '-.']
for (step, label_step), style in zip(step_configs, line_styles):
    hist = load_history(f'history_SD_{step}.csv')
    if hist is not None:
        ax.plot(hist[:, 0], hist[:, 1], 'o' + style, label=f'SD {label_step}', 
                markersize=4, linewidth=1.5, alpha=0.8)
        # Mark start and end points
        ax.plot(hist[0, 0], hist[0, 1], 'o', markersize=8, 
                markeredgewidth=2, markeredgecolor='black')
        ax.plot(hist[-1, 0], hist[-1, 1], '*', markersize=15, 
                markeredgewidth=1, markeredgecolor='black')

ax.set_xlabel('x1', fontsize=12)
ax.set_ylabel('x2', fontsize=12)
ax.set_title('Wykres 4: Metoda najszybszego spadku - różne długości kroku', fontsize=14, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
plt.savefig('wykres4_SD_all.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: wykres4_SD_all.png")

# WYKRES 5: Metoda CG, wszystkie kroki
print("Generating Plot 5: CG method, all step sizes...")
fig, ax = plt.subplots(figsize=(10, 8))
contour = ax.contour(X1, X2, Z, levels=30, cmap='viridis', alpha=0.6)
ax.clabel(contour, inline=True, fontsize=8, fmt='%.1f')

step_configs = [('0.05', 'h=0.05'), ('0.25', 'h=0.25'), ('variable', 'zmienny krok')]
line_styles = ['-', '--', '-.']
for (step, label_step), style in zip(step_configs, line_styles):
    hist = load_history(f'history_CG_{step}.csv')
    if hist is not None:
        ax.plot(hist[:, 0], hist[:, 1], 'o' + style, label=f'CG {label_step}', 
                markersize=4, linewidth=1.5, alpha=0.8)
        # Mark start and end points
        ax.plot(hist[0, 0], hist[0, 1], 'o', markersize=8, 
                markeredgewidth=2, markeredgecolor='black')
        ax.plot(hist[-1, 0], hist[-1, 1], '*', markersize=15, 
                markeredgewidth=1, markeredgecolor='black')

ax.set_xlabel('x1', fontsize=12)
ax.set_ylabel('x2', fontsize=12)
ax.set_title('Wykres 5: Metoda gradientów sprzężonych - różne długości kroku', fontsize=14, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
plt.savefig('wykres5_CG_all.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: wykres5_CG_all.png")

# WYKRES 6: Metoda Newtona, wszystkie kroki
print("Generating Plot 6: Newton method, all step sizes...")
fig, ax = plt.subplots(figsize=(10, 8))
contour = ax.contour(X1, X2, Z, levels=30, cmap='viridis', alpha=0.6)
ax.clabel(contour, inline=True, fontsize=8, fmt='%.1f')

step_configs = [('0.01', 'h=0.01'), ('0.0001', 'h=0.0001'), ('variable', 'zmienny krok')]
line_styles = ['-', '--', '-.']
for (step, label_step), style in zip(step_configs, line_styles):
    hist = load_history(f'history_Newton_{step}.csv')
    if hist is not None:
        ax.plot(hist[:, 0], hist[:, 1], 'o' + style, label=f'Newton {label_step}', 
                markersize=4, linewidth=1.5, alpha=0.8)
        # Mark start and end points
        ax.plot(hist[0, 0], hist[0, 1], 'o', markersize=8, 
                markeredgewidth=2, markeredgecolor='black')
        ax.plot(hist[-1, 0], hist[-1, 1], '*', markersize=15, 
                markeredgewidth=1, markeredgecolor='black')

ax.set_xlabel('x1', fontsize=12)
ax.set_ylabel('x2', fontsize=12)
ax.set_title('Wykres 6: Metoda Newtona - różne długości kroku', fontsize=14, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
plt.savefig('wykres6_Newton_all.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: wykres6_Newton_all.png")

print("\n✓ Wszystkie 6 wykresów zostały wygenerowane!")
print("  - wykres1_h005.png")
print("  - wykres2_h025.png")
print("  - wykres3_variable.png")
print("  - wykres4_SD_all.png")
print("  - wykres5_CG_all.png")
print("  - wykres6_Newton_all.png")
