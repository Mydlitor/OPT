#include <iostream>
#include <cmath>

using namespace std;

int main() {
    cout.precision(10);
    
    // Test different combinations
    double test_cases[][2] = {
        {25.0, 200.0},
        {25.0, 1000.0},
        {100.0, 200.0},
        {100.0, 1000.0},
        {50.0, 500.0}
    };
    
    const double P = 3000.0;        // N
    const double E = 120e9;         // Pa
    const double rho = 8920.0;      // kg/m^3
    
    for (int i = 0; i < 5; i++) {
        double d = test_cases[i][0];
        double l = test_cases[i][1];
        
        // Convert to meters
        double d_m = d / 1000.0;
        double l_m = l / 1000.0;
        
        // Calculate
        double mass = rho * M_PI * pow(d_m / 2.0, 2) * l_m;
        double u = (64.0 * P * pow(l_m, 3)) / (3.0 * E * M_PI * pow(d_m, 4));
        u = u * 1000.0;  // convert to mm
        double sigma = (32.0 * P * l_m) / (M_PI * pow(d_m, 3));
        
        cout << "d=" << d << "mm, l=" << l << "mm:" << endl;
        cout << "  m=" << mass << " kg, u=" << u << " mm, sigma=" << sigma/1e6 << " MPa" << endl;
    }
    
    return 0;
}
