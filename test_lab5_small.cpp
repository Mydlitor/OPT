#include "opt_alg.h"
#include <iostream>
#include <fstream>

using namespace std;

int main() {
    cout.precision(10);
    
    cout << "=== SMALL LAB5 TEST ===" << endl;
    
    double epsilon = 1e-3;
    int Nmax = 10000;
    
    // Test with a=1, just 3 values of w
    cout << "\n--- Testing with a = 1 ---" << endl;
    double weights[] = {0.0, 0.5, 1.0};
    
    for (int i = 0; i < 3; i++) {
        double w = weights[i];
        
        matrix x0(2, 1);
        x0(0) = 0.0;
        x0(1) = 0.0;
        
        matrix params(2, 1);
        params(0) = 1.0;  // a
        params(1) = w;     // weight
        
        solution::clear_calls();
        solution opt = Powell(ff5T, x0, epsilon, Nmax, params, NAN);
        
        matrix f1 = ff5_f1(opt.x, params, NAN);
        matrix f2 = ff5_f2(opt.x, params, NAN);
        
        cout << "w=" << w << ": x*=[" << opt.x(0) << "," << opt.x(1) 
             << "], f1*=" << f1(0) << ", f2*=" << f2(0) 
             << ", fcalls=" << solution::f_calls << endl;
    }
    
    // Test beam problem
    cout << "\n--- Testing beam problem ---" << endl;
    for (int i = 0; i < 3; i++) {
        double w = weights[i];
        
        matrix x0(2, 1);
        x0(0) = 50.0;   // d
        x0(1) = 500.0;  // l
        
        matrix params(1, 1);
        params(0) = w;
        
        solution::clear_calls();
        solution opt = Powell(ff5R, x0, epsilon, Nmax, params, NAN);
        
        double d = m2d(opt.x(0));
        double l = m2d(opt.x(1));
        
        const double P = 3000.0;
        const double E = 120e9;
        const double rho = 8920.0;
        
        double d_m = d / 1000.0;
        double l_m = l / 1000.0;
        
        double mass = rho * M_PI * pow(d_m / 2.0, 2) * l_m;
        double u = (64.0 * P * pow(l_m, 3)) / (3.0 * E * M_PI * pow(d_m, 4)) * 1000.0;
        double sigma = (32.0 * P * l_m) / (M_PI * pow(d_m, 3));
        
        cout << "w=" << w << ": d*=" << d << "mm, l*=" << l << "mm, "
             << "m=" << mass << "kg, u=" << u << "mm, sigma=" << sigma/1e6 << "MPa"
             << ", fcalls=" << solution::f_calls << endl;
    }
    
    return 0;
}
