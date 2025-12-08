# Lab 4: Gradient-Based Optimization

This lab implements gradient-based optimization methods and applies them to test functions and logistic regression.

## Implementation

### Core Algorithms

1. **Golden Section Line Search** (`golden()` in `opt_alg.cpp`)
   - Implements golden section search for 1D optimization
   - Used for line search in gradient-based methods

2. **Steepest Descent** (`SD()` in `opt_alg.cpp`)
   - Gradient descent with fixed step size or line search
   - Direction: d = -∇f(x)

3. **Conjugate Gradients** (`CG()` in `opt_alg.cpp`)
   - Fletcher-Reeves formula: β = ||∇f^(i)||² / ||∇f^(i-1)||²
   - Direction update: d = -∇f + β*d_prev

4. **Newton's Method** (`Newton()` in `opt_alg.cpp`)
   - Uses Hessian matrix
   - Direction: d = -H^(-1)∇f(x)

### Test Function

**Function**: f(x1,x2) = (1/6)x1^6 - 1.05*x1^4 + 2*x1^2 + x2^2 + x1*x2

Implemented in:
- `ff4T()`: Objective function
- `gf4T()`: Gradient
- `Hf4T()`: Hessian

### Logistic Regression

**Hypothesis**: h_θ(x) = 1/(1 + e^(-θ^T x))

**Cost Function**: J(θ) = -(1/m) Σ[y ln h + (1-y) ln(1-h)]

**Gradient**: ∂J/∂θ_j = (1/m) Σ (h(x^(i)) - y^(i)) x_j^(i)

Implemented in:
- `ff4R_cost()`: Cost function
- `gf4R_grad()`: Gradient

## Building and Running

### Compilation
```bash
g++ -std=c++17 -o opt_program main.cpp matrix.cpp solution.cpp ode_solver.cpp opt_alg.cpp user_funs.cpp
```

### Execution
```bash
./opt_program
```

### Output Files

1. **table1_lab4.csv**: Test function experiment results
   - Columns: Method, StepSize, AvgFinalF, AvgIterations, AvgFCalls, Convergences
   - 100 random starts in [-2,2]×[-2,2] for each method/step combination

2. **table2_lab4.csv**: Logistic regression results
   - Columns: Method, StepSize, Theta0, Theta1, Theta2, J, Accuracy, FCalls
   - Starting point: θ = [0,0,0]^T

3. **decision_boundary_lab4.csv**: Classification data
   - Columns: x1, x2, y, predicted
   - Used for plotting decision boundary

## Experimental Setup

### Test Function Experiments
- **Number of runs**: 100 per configuration
- **Starting points**: Random in [-2,2]×[-2,2]
- **Convergence criterion**: ||∇f|| < 10^(-3)
- **Max function calls**: 10,000
- **Random seed**: 42 (for reproducibility)

**Step sizes tested**:
- Steepest Descent: 0.05, 0.25
- Conjugate Gradients: 0.05, 0.25
- Newton's Method: 0.01, 0.0001

### Logistic Regression
- **Dataset**: 100 samples, 2 features
- **Data files**: XData.txt, YData.txt
- **Starting point**: θ = [0,0,0]^T
- **Methods**: SD, CG (fixed step and line search variants)

## Results

### Test Function
- **SD (h=0.05)**: Best performer with 100% convergence, avg f=0.131
- **Newton (h=0.01)**: Good convergence rate, avg f=0.554
- **CG**: Numerical instability issues with some starting points

### Logistic Regression
- **Best model**: SD with h=0.05
- **Accuracy**: 66%
- **Parameters**: θ = [-36.94, 2.09, -0.56]

## Known Issues

1. **CG Numerical Stability**: Conjugate gradients can produce NaN values for some starting points due to numerical instability in the Fletcher-Reeves update.

2. **CG Line Search**: CG with line search can be slow or fail to converge on the logistic regression problem.

3. **Newton Step Sizes**: Very small step sizes (0.0001) can prevent convergence within the iteration limit.

## Generating Plots

A Python script (`plot_lab4.py`) can be created to visualize:
- Convergence trajectories for test function
- Decision boundaries for logistic regression
- Classification accuracy comparison

## Future Improvements

1. Implement additional CG variants (Polak-Ribière, Hestenes-Stiefel)
2. Add Hessian computation for logistic regression to enable Newton's method
3. Implement adaptive step size strategies
4. Add restart mechanisms for CG to improve numerical stability
5. Implement BFGS or L-BFGS for large-scale problems
