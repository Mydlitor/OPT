#include "opt_alg.h"
#include <iostream>

using namespace std;

// Simple test function: f(x) = (x1-2)^2 + (x2-3)^2
matrix test_func(matrix x, matrix ud1, matrix ud2) {
    double x1 = m2d(x(0));
    double x2 = m2d(x(1));
    return matrix((x1-2)*(x1-2) + (x2-3)*(x2-3));
}

int main() {
    cout.precision(10);
    
    matrix x0(2, 1);
    x0(0) = 0.0;
    x0(1) = 0.0;
    
    cout << "Testing Powell method on f(x) = (x1-2)^2 + (x2-3)^2" << endl;
    cout << "Starting point: x0 = [" << x0(0) << ", " << x0(1) << "]" << endl;
    
    solution::clear_calls();
    solution opt = Powell(test_func, x0, 1e-6, 10000, NAN, NAN);
    
    cout << "Optimal point: x* = [" << opt.x(0) << ", " << opt.x(1) << "]" << endl;
    cout << "Optimal value: f* = " << opt.y(0) << endl;
    cout << "Function calls: " << solution::f_calls << endl;
    cout << "Expected: x* = [2, 3], f* = 0" << endl;
    
    return 0;
}
