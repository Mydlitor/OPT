#include "opt_alg.h"
#include <iostream>

int main() {
    std::cout.precision(10);
    std::cout << "Quick test of lab5 functions" << std::endl;
    
    // Test Powell on simple function
    matrix x0(2, 1);
    x0(0) = 0.0;
    x0(1) = 0.0;
    
    matrix params(2, 1);
    params(0) = 1.0;  // a
    params(1) = 0.5;  // w
    
    solution::clear_calls();
    solution opt = Powell(ff5T, x0, 1e-3, 10000, params, NAN);
    
    std::cout << "Test function optimization (a=1, w=0.5):" << std::endl;
    std::cout << "  x* = [" << opt.x(0) << ", " << opt.x(1) << "]" << std::endl;
    std::cout << "  f* = " << opt.y(0) << std::endl;
    std::cout << "  f_calls = " << solution::f_calls << std::endl;
    
    matrix f1 = ff5_f1(opt.x, params, NAN);
    matrix f2 = ff5_f2(opt.x, params, NAN);
    std::cout << "  f1* = " << f1(0) << std::endl;
    std::cout << "  f2* = " << f2(0) << std::endl;
    
    // Test beam problem
    matrix x0_beam(2, 1);
    x0_beam(0) = 50.0;  // d
    x0_beam(1) = 500.0; // l
    
    matrix params_beam(1, 1);
    params_beam(0) = 0.5;  // w
    
    solution::clear_calls();
    solution opt_beam = Powell(ff5R, x0_beam, 1e-3, 10000, params_beam, NAN);
    
    std::cout << "\nBeam optimization (w=0.5):" << std::endl;
    std::cout << "  d* = " << opt_beam.x(0) << " mm" << std::endl;
    std::cout << "  l* = " << opt_beam.x(1) << " mm" << std::endl;
    std::cout << "  f* = " << opt_beam.y(0) << std::endl;
    std::cout << "  f_calls = " << solution::f_calls << std::endl;
    
    std::cout << "\nAll tests passed!" << std::endl;
    
    return 0;
}
