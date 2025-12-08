#!/usr/bin/env python3
"""
Plotting script for Lab 4 results
Generates visualizations for gradient-based optimization experiments
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Set style
plt.style.use('seaborn-v0_8-darkgrid')

def plot_test_function_results():
    """Plot results from test function experiments"""
    df = pd.read_csv('table1_lab4.csv')
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))
    
    # Plot 1: Average function values
    methods = []
    labels = []
    values = []
    colors = []
    
    color_map = {'SD': 'blue', 'CG': 'green', 'Newton': 'red'}
    
    for _, row in df.iterrows():
        method = row['Method']
        step = row['StepSize']
        avg_f = row['AvgFinalF']
        
        if not np.isnan(avg_f):
            methods.append(f"{method}\n(h={step})")
            values.append(avg_f)
            colors.append(color_map.get(method, 'gray'))
    
    ax1.bar(range(len(methods)), values, color=colors, alpha=0.7)
    ax1.set_xticks(range(len(methods)))
    ax1.set_xticklabels(methods, rotation=45, ha='right')
    ax1.set_ylabel('Average Final f(x)')
    ax1.set_title('Test Function: Average Final Values')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Function calls
    methods = []
    fcalls = []
    colors = []
    
    for _, row in df.iterrows():
        method = row['Method']
        step = row['StepSize']
        avg_fc = row['AvgFCalls']
        
        methods.append(f"{method}\n(h={step})")
        fcalls.append(avg_fc)
        colors.append(color_map.get(method, 'gray'))
    
    ax2.bar(range(len(methods)), fcalls, color=colors, alpha=0.7)
    ax2.set_xticks(range(len(methods)))
    ax2.set_xticklabels(methods, rotation=45, ha='right')
    ax2.set_ylabel('Average Function Calls')
    ax2.set_title('Test Function: Computational Cost')
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Convergence rates
    methods = []
    conv_rates = []
    colors = []
    
    for _, row in df.iterrows():
        method = row['Method']
        step = row['StepSize']
        convs = row['Convergences']
        
        methods.append(f"{method}\n(h={step})")
        conv_rates.append(convs)
        colors.append(color_map.get(method, 'gray'))
    
    ax3.bar(range(len(methods)), conv_rates, color=colors, alpha=0.7)
    ax3.set_xticks(range(len(methods)))
    ax3.set_xticklabels(methods, rotation=45, ha='right')
    ax3.set_ylabel('Successful Convergences (out of 100)')
    ax3.set_title('Test Function: Convergence Success Rate')
    ax3.set_ylim([0, 105])
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Summary table
    ax4.axis('off')
    table_data = []
    for _, row in df.iterrows():
        table_data.append([
            row['Method'],
            f"{row['StepSize']:.4f}" if not np.isnan(row['StepSize']) else 'N/A',
            f"{row['AvgFinalF']:.4f}" if not np.isnan(row['AvgFinalF']) else 'NaN',
            f"{row['AvgFCalls']:.1f}",
            f"{row['Convergences']}/100"
        ])
    
    table = ax4.table(cellText=table_data,
                     colLabels=['Method', 'Step', 'Avg f', 'Avg Calls', 'Conv.'],
                     cellLoc='center',
                     loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2)
    ax4.set_title('Test Function: Summary Table', pad=20)
    
    plt.tight_layout()
    plt.savefig('test_function_results.png', dpi=150, bbox_inches='tight')
    print("Saved: test_function_results.png")
    plt.close()

def plot_decision_boundary():
    """Plot decision boundary for logistic regression"""
    df = pd.read_csv('decision_boundary_lab4.csv')
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Separate by actual class
    class0 = df[df['y'] == 0]
    class1 = df[df['y'] == 1]
    
    # Plot actual classes
    ax.scatter(class0['x1'], class0['x2'], c='red', marker='o', 
              s=50, alpha=0.6, label='Class 0 (actual)')
    ax.scatter(class1['x1'], class1['x2'], c='blue', marker='s', 
              s=50, alpha=0.6, label='Class 1 (actual)')
    
    # Mark misclassifications
    misclass = df[df['y'] != df['predicted']]
    ax.scatter(misclass['x1'], misclass['x2'], 
              facecolors='none', edgecolors='black', 
              s=200, linewidths=2, label='Misclassified')
    
    # Try to plot decision boundary
    # Decision boundary: theta0 + theta1*x1 + theta2*x2 = 0
    # From table2_lab4.csv, best model (SD h=0.05): theta = [-36.94, 2.09, -0.56]
    try:
        theta0, theta1, theta2 = -36.94474329, 2.094255867, -0.5558911932
        
        # Create grid
        x1_min, x1_max = df['x1'].min() - 5, df['x1'].max() + 5
        x2_min, x2_max = df['x2'].min() - 5, df['x2'].max() + 5
        
        # Decision boundary: x2 = -(theta0 + theta1*x1) / theta2
        if abs(theta2) > 1e-10:
            x1_line = np.linspace(x1_min, x1_max, 100)
            x2_line = -(theta0 + theta1 * x1_line) / theta2
            ax.plot(x1_line, x2_line, 'g-', linewidth=2, label='Decision Boundary')
            
            # Shade regions
            ax.fill_between(x1_line, x2_line, x2_max, alpha=0.1, color='blue')
            ax.fill_between(x1_line, x2_min, x2_line, alpha=0.1, color='red')
    except Exception as e:
        print(f"Could not plot decision boundary: {e}")
    
    ax.set_xlabel('Feature x1', fontsize=12)
    ax.set_ylabel('Feature x2', fontsize=12)
    ax.set_title('Logistic Regression: Decision Boundary (Best Model: SD h=0.05)', 
                fontsize=14)
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('decision_boundary.png', dpi=150, bbox_inches='tight')
    print("Saved: decision_boundary.png")
    plt.close()

def plot_logistic_regression_results():
    """Plot logistic regression comparison"""
    df = pd.read_csv('table2_lab4.csv')
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Filter out N/A rows
    df_valid = df[df['Theta0'] != 'N/A'].copy()
    
    # Extract method names for x-axis
    methods = []
    for _, row in df_valid.iterrows():
        method = row['Method'].replace('_', ' ')
        step = row['StepSize']
        if step == 0:
            methods.append(f"{method}\n(line search)")
        else:
            methods.append(f"{method}\n(h={step})")
    
    # Plot 1: Accuracy comparison
    accuracies = df_valid['Accuracy'].values
    colors = ['green' if a == max(accuracies) else 'skyblue' for a in accuracies]
    
    ax1.bar(range(len(methods)), accuracies, color=colors, alpha=0.7)
    ax1.set_xticks(range(len(methods)))
    ax1.set_xticklabels(methods, rotation=45, ha='right', fontsize=9)
    ax1.set_ylabel('Classification Accuracy', fontsize=12)
    ax1.set_title('Logistic Regression: Accuracy Comparison', fontsize=14)
    ax1.set_ylim([0, 1])
    ax1.axhline(y=0.5, color='r', linestyle='--', alpha=0.3, label='Random guess')
    ax1.grid(True, alpha=0.3, axis='y')
    ax1.legend()
    
    # Add value labels on bars
    for i, (m, a) in enumerate(zip(methods, accuracies)):
        ax1.text(i, a + 0.02, f'{a:.2f}', ha='center', va='bottom', fontsize=9)
    
    # Plot 2: Function calls comparison
    fcalls = df_valid['FCalls'].values
    
    ax2.bar(range(len(methods)), fcalls, color='coral', alpha=0.7)
    ax2.set_xticks(range(len(methods)))
    ax2.set_xticklabels(methods, rotation=45, ha='right', fontsize=9)
    ax2.set_ylabel('Function Calls', fontsize=12)
    ax2.set_title('Logistic Regression: Computational Cost', fontsize=14)
    ax2.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig('logistic_regression_comparison.png', dpi=150, bbox_inches='tight')
    print("Saved: logistic_regression_comparison.png")
    plt.close()

def plot_test_function_contour():
    """Plot contour of test function"""
    # f(x1,x2) = (1/6)*x1^6 - 1.05*x1^4 + 2*x1^2 + x2^2 + x1*x2
    
    x1 = np.linspace(-2, 2, 200)
    x2 = np.linspace(-2, 2, 200)
    X1, X2 = np.meshgrid(x1, x2)
    
    F = (1/6)*X1**6 - 1.05*X1**4 + 2*X1**2 + X2**2 + X1*X2
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Contour plot
    levels = np.logspace(-1, 2, 20)
    contour = ax.contour(X1, X2, F, levels=levels, cmap='viridis', alpha=0.6)
    ax.clabel(contour, inline=True, fontsize=8)
    
    # Fill contour
    contourf = ax.contourf(X1, X2, F, levels=levels, cmap='viridis', alpha=0.3)
    
    # Colorbar
    cbar = plt.colorbar(contourf, ax=ax)
    cbar.set_label('f(x1, x2)', fontsize=12)
    
    # Mark some special points
    # Origin
    ax.plot(0, 0, 'r*', markersize=15, label='Origin')
    
    ax.set_xlabel('x1', fontsize=12)
    ax.set_ylabel('x2', fontsize=12)
    ax.set_title('Test Function: f(x1,x2) = (1/6)x1⁶ - 1.05x1⁴ + 2x1² + x2² + x1·x2', 
                fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xlim([-2, 2])
    ax.set_ylim([-2, 2])
    
    plt.tight_layout()
    plt.savefig('test_function_contour.png', dpi=150, bbox_inches='tight')
    print("Saved: test_function_contour.png")
    plt.close()

if __name__ == '__main__':
    print("Generating plots for Lab 4...")
    
    try:
        plot_test_function_results()
    except Exception as e:
        print(f"Error plotting test function results: {e}")
    
    try:
        plot_test_function_contour()
    except Exception as e:
        print(f"Error plotting test function contour: {e}")
    
    try:
        plot_decision_boundary()
    except Exception as e:
        print(f"Error plotting decision boundary: {e}")
    
    try:
        plot_logistic_regression_results()
    except Exception as e:
        print(f"Error plotting logistic regression results: {e}")
    
    print("\nAll plots generated successfully!")
