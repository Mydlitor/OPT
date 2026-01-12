# Lab5 Implementation Summary

## Overview
This implementation adds multi-criteria optimization capabilities using Powell's method to the OPT repository.

## Components Implemented

### 1. Powell's Method (`opt_alg.cpp`)
- **Function**: `solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)`
- **Algorithm**: Follows the pseudocode from lab materials
- **Features**:
  - Iterative optimization along conjugate directions
  - Uses golden section search for line minimization
  - Adaptive direction updates
  - Convergence based on position norm and function call limit
- **Optimizations**: Fixed memory leak, appropriate line search intervals

### 2. Test Functions (`user_funs.cpp`)

#### Test Problem Functions
- **`ff5_f1`**: First criterion `f₁(x) = a((x₁-3)² + (x₂-3)²)`
- **`ff5_f2`**: Second criterion `f₂(x) = (1/a)((x₁+3)² + (x₂+3)²)`  
- **`ff5T`**: Weighted test function `f(x) = w·f₁(x) + (1-w)·f₂(x)`
- **Parameters**: 
  - `a ∈ {1, 10, 100}` - scaling parameter
  - `w ∈ [0, 1]` - weight for multi-criteria optimization

#### Real-World Problem Function
- **`ff5R`**: Cantilever beam optimization
- **Objectives**:
  - Minimize mass: `m = ρπ(d/2)²l`
  - Minimize deflection: `u = (64Pl³)/(3Eπd⁴)`
- **Constraints**:
  - Deflection: `u ≤ 2.5 mm`
  - Stress: `σ = (32Pl)/(πd³) ≤ 300 MPa`
  - Bounds: `d ∈ [0.01, 1000] mm`, `l ∈ [0.2, 1000] mm`
- **Implementation**: Penalty function approach with balanced coefficients

### 3. Main Program (`main.cpp`)

#### `lab5()` Function
- **Test Problem**: 
  - 101 optimizations for each value of `a` (1, 10, 100)
  - Weight `w` varies from 0.00 to 1.00 in steps of 0.01
  - Random starting points in `[-10, 10] × [-10, 10]`
  - Outputs: `lab5_test_a1.csv`, `lab5_test_a10.csv`, `lab5_test_a100.csv`

- **Real Problem**:
  - 101 optimizations for cantilever beam
  - Weight `w` varies from 0.00 to 1.00 in steps of 0.01
  - Random starting points in feasible region
  - Output: `lab5_beam.csv`

## Testing Results

### Validation Tests
1. **Simple quadratic**: Converges to `x* = [2, 3]` with `f* ≈ 0` ✓
2. **Multi-criteria test**: Finds Pareto optimal solutions for all `a` values ✓
3. **Beam problem**: Converges without numerical overflow/NaN ✓

### Code Quality
- **Compilation**: Clean compilation with no errors
- **Code Review**: All issues addressed (memory leaks, numerical stability, variable handling)
- **Security**: No vulnerabilities detected by CodeQL

## Usage

### Running Lab5
```cpp
// In main.cpp, uncomment:
lab5();
```

### Compilation
```bash
g++ -std=c++17 -o opt_program main.cpp opt_alg.cpp user_funs.cpp solution.cpp matrix.cpp ode_solver.cpp -lm
./opt_program
```

### Output Files
- `lab5_test_a1.csv` - Pareto front for a=1
- `lab5_test_a10.csv` - Pareto front for a=10
- `lab5_test_a100.csv` - Pareto front for a=100
- `lab5_beam.csv` - Pareto front for beam problem

Each CSV contains:
- Weight parameter `w`
- Starting point
- Optimal point
- Individual criteria values
- Combined objective value
- Function call count

## Implementation Notes

### Numerical Stability
- Penalty coefficients tuned to avoid overflow (1e5 for bounds, 1e4 for constraints)
- Line search intervals adjusted to problem scale ([-10, 10])
- Proper variable clamping in penalty functions

### Memory Management
- All dynamically allocated arrays properly deallocated
- No memory leaks in Powell method

### Reproducibility
- Random seed (42) used for consistent results
- Can be changed to `time(NULL)` for true randomness

## Known Behavior

### Beam Optimization
The beam problem tends to converge toward minimum length (lower bound) for most weight values. This is expected because:
- Both mass and deflection decrease with decreasing length
- The optimizer finds the minimum feasible length that satisfies constraints
- This represents valid Pareto-optimal solutions

## Future Enhancements

Potential improvements for users:
1. Add gradient information to Powell method for faster convergence
2. Implement constraint handling methods beyond penalties
3. Add visualization tools for Pareto fronts
4. Support for more than 2 objectives
