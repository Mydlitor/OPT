# Lab5 Implementation Summary

## Overview
This implementation adds multi-criteria optimization capabilities using Powell's method to the OPT repository, following specifications from K5.pdf.

## Components Implemented

### 1. Powell's Method (`opt_alg.cpp`)
- **Function**: `solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)`
- **Algorithm**: Follows the pseudocode from K5.pdf lab materials
- **Features**:
  - Iterative optimization along conjugate directions
  - Uses golden section search for line minimization
  - Adaptive direction updates
  - Convergence based on position norm and function call limit
- **Line Search**: Golden section method with expansion for initial interval determination

### 2. Test Functions (`user_funs.cpp`)

#### Test Problem Functions
- **`ff5_f1`**: First criterion `f₁(x) = a((x₁-3)² + (x₂-3)²)`
- **`ff5_f2`**: Second criterion `f₂(x) = (1/a)((x₁+3)² + (x₂+3)²)`  
- **`ff5T`**: Weighted test function `f(x) = w·f₁(x) + (1-w)·f₂(x)`
- **Parameters**: 
  - `a ∈ {1, 10, 100}` - scaling parameter
  - `w ∈ [0, 1]` - weight for multi-criteria optimization

#### Real-World Problem Function (Cantilever Beam)
- **`ff5R`**: Cantilever beam optimization
- **Physical Parameters** (from K5.pdf):
  - Force: `P = 2 kN` (2000 N)
  - Young's modulus: `E = 120 GPa`
  - Density: `ρ = 8920 kg/m³`
- **Objectives**:
  - Minimize mass: `m = ρπ(d/2)²l`
  - Minimize deflection: `u = (64Pl³)/(3Eπd⁴)`
- **Constraints**:
  - Deflection: `u ≤ 2.5 mm`
  - Stress: `σ = (32Pl)/(πd³) ≤ 300 MPa`
  - Bounds: `d ∈ [10, 50] mm`, `l ∈ [200, 1000] mm`
- **Implementation**: External penalty function approach with coefficient 1e6
- **Validation** (l=500mm, d=25mm):
  - Mass: 2.189 kg ≈ 2.19 kg ✓
  - Deflection: 36.217 mm ≈ 36.22 mm ✓
  - Stress: 651.9 MPa ✓

### 3. Main Program (`main.cpp`)

#### `lab5()` Function
- **Test Problem**: 
  - 101 optimizations for each value of `a` (1, 10, 100)
  - Weight `w` varies from 0.00 to 1.00 in steps of 0.01
  - Random starting points in `[-10, 10] × [-10, 10]`
  - Outputs: `lab5_test_a1.csv`, `lab5_test_a10.csv`, `lab5_test_a100.csv`

- **Real Problem** (Cantilever Beam):
  - 101 optimizations for cantilever beam
  - Weight `w` varies from 0.00 to 1.00 in steps of 0.01
  - Random starting points in feasible region `[10, 50] × [200, 1000]`
  - Output: `lab5_beam.csv`

## Algorithm Details

### Weighted Criterion Method
Converts multi-criteria problem to single-criterion:
```
f(x) = w·f₁(x) + (1-w)·f₂(x)
```
where `w ∈ [0, 1]`

### Powell's Method
- Optimizes along n mutually conjugate directions
- Uses golden section search for 1D line minimization
- Updates direction vectors iteratively
- Converges when `||x_{i+1} - x_i|| < ε`

### External Penalty Function
For constraints `g_j(x) ≤ 0`:
```
Φ(x) = f(x) + c·∑max(0, g_j(x))²
```
where `c = 1e6` (penalty coefficient)

## Testing Results

### Validation Tests
1. **Beam validation** (d=25mm, l=500mm):
   - Mass: 2.189 kg (expected 2.19 kg) ✓
   - Deflection: 36.217 mm (expected 36.22 mm) ✓
   - Stress: 651.9 MPa (expected 651.9 MPa) ✓

2. **Multi-criteria test**: Finds Pareto optimal solutions for all `a` values ✓
3. **Convergence**: All optimizations complete within Nmax ✓

### Code Quality
- **Compilation**: Clean compilation with no errors
- **Code Review**: All issues addressed (no duplicate variables, named constants used)
- **Security**: No vulnerabilities detected by CodeQL
- **Memory Management**: No memory leaks

## Usage

### Running Lab5
```cpp
// In main.cpp:
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
- Starting point coordinates
- Optimal point coordinates
- Individual criteria values (f1, f2 or mass, deflection)
- Combined objective value
- Additional metrics (stress for beam problem)
- Function call count

## Implementation Notes

### Numerical Stability
- Penalty coefficient 1e6 provides good constraint enforcement without overflow
- Line search intervals [-10, 10] match problem scale
- Proper variable clamping in penalty functions
- Stress normalized to MPa scale in penalty calculation

### Memory Management
- All dynamically allocated arrays properly deallocated in Powell method
- No memory leaks detected

### Reproducibility
- Random seed (42) used for consistent results across runs
- Can be changed to `time(NULL)` for true randomness

### Constants and Bounds
- All constants defined as named variables (no magic numbers)
- Bounds centralized in ff5R function
- Values match K5.pdf specifications exactly

## Known Behavior

### Beam Optimization
- Most weight values converge toward minimum length (200mm)
- This is expected because both mass and deflection decrease with decreasing length
- The optimizer finds minimum feasible length satisfying all constraints
- Represents valid Pareto-optimal solutions along the constraint boundary

### Pareto Front
- Test problems show smooth transition from f1 minimum to f2 minimum
- Beam problem shows trade-off between mass and deflection
- All 101 weight values produce distinct solutions forming the Pareto front

## Compliance with K5.pdf

✅ **All requirements implemented**:
1. Test function with a = 1, 10, 100
2. Real problem (cantilever beam) with correct parameters
3. Powell's method for optimization
4. Golden section search for line minimization
5. Expansion method available (in opt_alg.cpp)
6. External penalty function for constraints
7. 101 optimizations for w = 0.00, 0.01, ..., 1.00
8. Results in CSV format
9. Validation test matches expected values

## Future Enhancements

Potential improvements:
1. Add visualization tools for Pareto fronts (matplotlib/gnuplot)
2. Implement interior penalty method for comparison
3. Add gradient information to Powell for faster convergence
4. Support for more than 2 objectives (3D Pareto surfaces)
5. Interactive parameter adjustment
6. Excel file generation with formatting
