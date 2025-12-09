# Generating Contour Plots for Optimization Methods

This document describes how to generate contour plots showing the optimization paths for different gradient-based methods (Steepest Descent, Conjugate Gradients, and Newton's method) with various step sizes.

## Overview

The implementation generates 6 contour plots showing how different optimization methods navigate the test function `ff4T`:

```
f(x1, x2) = (1/6)*x1^6 - 1.05*x1^4 + 2*x1^2 + x2^2 + x1*x2
```

## Files Generated

### CSV History Files

When you run the C++ program, it generates 11 CSV files containing iteration histories:

- `history_SD_0.05.csv` - Steepest Descent with step size 0.05
- `history_SD_0.25.csv` - Steepest Descent with step size 0.25
- `history_SD_variable.csv` - Steepest Descent with line search (variable step)
- `history_CG_0.05.csv` - Conjugate Gradients with step size 0.05
- `history_CG_0.25.csv` - Conjugate Gradients with step size 0.25
- `history_CG_variable.csv` - Conjugate Gradients with line search (variable step)
- `history_Newton_0.05.csv` - Newton's method with step size 0.05
- `history_Newton_0.25.csv` - Newton's method with step size 0.25
- `history_Newton_0.01.csv` - Newton's method with step size 0.01
- `history_Newton_0.0001.csv` - Newton's method with step size 0.0001
- `history_Newton_variable.csv` - Newton's method with line search (variable step)

Each CSV file has the format:
```
iteration,x1,x2
0,<x1_value>,<x2_value>
1,<x1_value>,<x2_value>
...
```

### Plot Files

After running the Python script, 6 PNG image files are created:

1. **wykres1_h005.png** - All three methods with fixed step h=0.05
2. **wykres2_h025.png** - All three methods with fixed step h=0.25
3. **wykres3_variable.png** - All three methods with variable step (line search)
4. **wykres4_SD_all.png** - Steepest Descent with all step sizes
5. **wykres5_CG_all.png** - Conjugate Gradients with all step sizes
6. **wykres6_Newton_all.png** - Newton's method with all step sizes

## How to Use

### Step 1: Compile the C++ Program

If using g++:
```bash
g++ -std=c++11 -o opt_program main.cpp matrix.cpp solution.cpp opt_alg.cpp ode_solver.cpp user_funs.cpp
```

Or if using Visual Studio, open the solution file `OPT.sln` and build the project.

### Step 2: Run the Program

```bash
./opt_program
```

or on Windows:
```bash
opt_program.exe
```

The program will:
1. Perform 100 optimization experiments for each method/step combination
2. Save results to `table1_lab4.csv` and `table2_lab4.csv`
3. Generate iteration history CSV files for plotting
4. Display progress messages

### Step 3: Generate the Plots

Make sure you have the required Python packages installed:
```bash
pip install numpy matplotlib pandas
```

Then run the plotting script:
```bash
python3 plot_contours.py
```

The script will:
1. Load the iteration history CSV files
2. Generate 6 contour plots with optimization paths
3. Save the plots as PNG files
4. Display confirmation messages

## Understanding the Plots

- **Contour lines** show the values of the objective function
- **Colored paths** show the trajectory of each optimization method
- **Circle (○)** marks the starting point
- **Star (★)** marks the final converged point
- **Line colors:**
  - Red: Steepest Descent (SD)
  - Blue: Conjugate Gradients (CG)
  - Green: Newton's Method

## Implementation Details

### Code Changes Made

1. **opt_alg.h**: Added overloaded function signatures with history tracking parameter
2. **opt_alg.cpp**: Implemented overloaded versions of SD, CG, and Newton that record position at each iteration
3. **main.cpp**: Added section in `lab4()` to generate history data using a single consistent starting point
4. **plot_contours.py**: Created Python script to visualize the optimization paths

### Key Design Decisions

- **Backward Compatibility**: Original function signatures remain unchanged; history tracking is an optional feature via function overloading
- **Consistent Starting Point**: All methods use the same random starting point (the first from the generated list) to enable fair comparison
- **History Recording**: Positions are recorded at the start and after each iteration update

## Troubleshooting

**Problem**: CSV files not generated
- **Solution**: Make sure the C++ program ran successfully and completed without errors

**Problem**: Python script fails with "File not found"
- **Solution**: Ensure you're running the Python script from the same directory as the CSV files

**Problem**: Missing Python packages
- **Solution**: Install required packages: `pip install numpy matplotlib pandas`

**Problem**: Plots look empty or incomplete
- **Solution**: Check that the optimization converged and generated iteration data in the CSV files

## Notes

- The starting point is randomly generated with seed 42 for reproducibility
- All methods use epsilon = 1e-3 as convergence tolerance
- Maximum iterations is set to 10,000 to prevent infinite loops
- Line search uses golden section method for variable step size optimization
