#include <iostream>
#include <cmath>

using namespace std;

int main() {
    cout.precision(10);
    
    // Test case: l=1000mm, d=25mm
    double d = 25.0;  // mm
    double l = 1000.0;  // mm
    
    const double P = 3000.0;        // N
    const double E = 120e9;         // Pa
    const double rho = 8920.0;      // kg/m^3
    
    // Convert to meters
    double d_m = d / 1000.0;
    double l_m = l / 1000.0;
    
    // Calculate
    double mass = rho * M_PI * pow(d_m / 2.0, 2) * l_m;
    double u = (64.0 * P * pow(l_m, 3)) / (3.0 * E * M_PI * pow(d_m, 4));
    u = u * 1000.0;  // convert to mm
    double sigma = (32.0 * P * l_m) / (M_PI * pow(d_m, 3));
    
    cout << "Validation test (l=" << l << "mm, d=" << d << "mm):" << endl;
    cout << "  m = " << mass << " kg (expected ~2.19 kg)" << endl;
    cout << "  u = " << u << " mm (expected ~36.22 mm)" << endl;
    cout << "  sigma = " << sigma / 1e6 << " MPa (expected ~651.9 MPa)" << endl;
    
    return 0;
}
