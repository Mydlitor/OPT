#include <iostream>
#include <cmath>

using namespace std;

int main() {
    cout.precision(10);
    
    const double P = 3000.0;        // N
    const double E = 120e9;         // Pa
    const double rho = 8920.0;      // kg/m^3
    
    // Try to find d and l that give the expected results
    // Expected: m ≈ 2.19 kg, u ≈ 36.22 mm, σ ≈ 651.9 MPa
    
    // Test with d in 20-30mm range and l in 400-600mm range
    for (double d = 20.0; d <= 30.0; d += 1.0) {
        for (double l = 400.0; l <= 600.0; l += 20.0) {
            double d_m = d / 1000.0;
            double l_m = l / 1000.0;
            
            double mass = rho * M_PI * pow(d_m / 2.0, 2) * l_m;
            double u = (64.0 * P * pow(l_m, 3)) / (3.0 * E * M_PI * pow(d_m, 4)) * 1000.0;
            double sigma = (32.0 * P * l_m) / (M_PI * pow(d_m, 3));
            
            // Check if close to expected
            if (fabs(mass - 2.19) < 0.1 && fabs(u - 36.22) < 2.0 && fabs(sigma/1e6 - 651.9) < 50.0) {
                cout << "MATCH: d=" << d << "mm, l=" << l << "mm:" << endl;
                cout << "  m=" << mass << " kg, u=" << u << " mm, sigma=" << sigma/1e6 << " MPa" << endl;
            }
        }
    }
    
    return 0;
}
