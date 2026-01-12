#include"user_funs.h"

matrix ff0T(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla przypadku testowego
{
	matrix y;												// y zawiera warto�� funkcji celu
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);		// ud1 zawiera wsp�rz�dne szukanego optimum
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla problemu rzeczywistego
{
	matrix y;												// y zawiera warto�� funkcji celu
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz�tkowe
		MT = matrix(2, new double[2] { m2d(x), 0.5 });		// MT zawiera moment si�y dzia�aj�cy na wahad�o oraz czas dzia�ania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);	// rozwi�zujemy r�wnanie r�niczkowe
	int n = get_len(Y[0]);									// d�ugo�� rozwi�zania
	double teta_max = Y[1](0, 0);							// szukamy maksymalnego wychylenia wahad�a
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));							// warto�� funkcji celu (ud1 to za�o�one maksymalne wychylenie)
	Y[0].~matrix();											// usuwamy z pami�ci rozwi�zanie RR
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);										// definiujemy wektor pochodnych szukanych funkcji
	double m = 1, l = 0.5, b = 0.5, g = 9.81;				// definiujemy parametry modelu
	double I = m * pow(l, 2);
	dY(0) = Y(1);																// pochodna z po�o�enia to pr�dko��
	dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;	// pochodna z pr�dko�ci to przyspieszenie
	return dY;
}

matrix ff1T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	double exponent = -1 * pow((0.1 * m2d(x) - 2 * M_PI), 2);
	y = -cos(0.1 * m2d(x)) * exp(exponent) + 0.002 * pow((0.1 * m2d(x)), 2);
	return y;
}


matrix df1(double x, matrix Y, matrix ud1, matrix ud2)
{
    // Y = [Va, Vb, Tb]
    matrix dY(3, 1);
    double a = 0.98, b = 0.63, g = 9.81;
    double Pa = 2.0, Pb = 1.0;
    double Ta0 = 95.0;                 // temperatura w zbiorniku A
    double TBin = 20.0;                // temperatura dopływu z zewnątrz
    double Fin = 0.01;                // dopływ z zewnątrz [m^3/s]
    double Db = 0.00365665;           // [m^2]

    double Da = m2d(ud1)*0.0001;             // [m^2] pole otworu między A a B


    double As_S = (2 * g * Y(0)) / Pa;
	if (As_S < 0) 
        As_S = 0.0;
    double Fa_out = a * b * Da * sqrt(As_S);   // wypływ z A do B
    
    As_S = (2 * g * Y(1)) / Pb;
	if (As_S < 0)
		As_S = 0.0;

    double Fb_out = a * b * Db * sqrt(As_S);    // wypływ z B na zewnątrz

    // Równania różniczkowe:
    dY(0) = -Fa_out;
    dY(1) = Fa_out - Fb_out + Fin;
    dY(2) = (Fin/Y(1))*(TBin - Y(2)) + (Fa_out/Y(1))*(Ta0 - Y(2));

    return dY;
}


matrix ff1R(matrix x, matrix ud1, matrix ud2)
{
    // x - wektor argumentów (x(0) = D_A w cm^2)
    // ud1, ud2 - nieużywane (opcjonalne)
    
    matrix y(1,1);
    matrix Y0 = matrix(3, 1);
    Y0(0) = 5.0;    // Va0 [m^3]
    Y0(1) = 1.0;    // Vb0 [m^3]
    Y0(2) = 20.0;   // Tb0 [°C]
    
    // Rozwiązanie ODE
    matrix* Y = solve_ode(df1, 0.0, 1.0, 2000.0, Y0, x, NAN);

    // Poprawka: długość (liczba kroków) pobierana z wektora czasu Y[0]
    int n = get_len(Y[0]); // liczba wierszy (kroków czasowych)
    
    // Szukamy maksymalnej temperatury w zbiorniku B -> trzecia zmienna stanu to kolumna 2
    double Tmax = Y[1](0, 2);
    for (int i = 1; i < n; ++i)
        if (Y[1](i, 2) > Tmax)
            Tmax = Y[1](i, 2);
    
    //debug
    //cout << "tmax: " << Tmax << endl;

    // Funkcja celu: różnica względem 50°C (ud1 zawiera cel temperatury)
    y(0,0) = fabs(Tmax - m2d(ud1));
    
    // Zwalniamy pamięć macierzy zgodnie z innymi funkcjami w projekcie
    Y[0].~matrix();
    Y[1].~matrix();
    
    return y;
}

matrix ff2T(matrix x, matrix ud1, matrix ud2) {
    double x1 = m2d(x(0));
    double x2 = m2d(x(1));

    return matrix(x1*x1 + x2*x2 - cos(2.5*M_PI*x1) - cos(2.5*M_PI*x2) + 2);
}

matrix df2(double t, matrix Y, matrix ud1, matrix ud2)
{
    matrix dY(2, 1);
    
    double m = 1.0;     // masa ramienia
    double mc = 5.0;    // masa ciężarka
    double l = 2.0;     // długość ramienia
    double b = 0.25;    // współczynnik tarcia

    // Moment bezwładności
    double I = m * l * l / 3.0 + mc * l * l;
    
    // Parametry regulatora
    double k1 = m2d(ud1(0));
    double k2 = m2d(ud1(1));
    
    double alpha_ref = m2d(ud2(0));
    double omega_ref = m2d(ud2(1));

    double M = k1 * (alpha_ref - Y(0)) + k2 * (omega_ref - Y(1));

    dY(0) = Y(1);
    dY(1) = (M - b * Y(1)) / I;
    
    return dY;
}

matrix ff2R(matrix x, matrix ud1, matrix ud2)
{
    matrix y;

    double k1 = m2d(x(0));
    double k2 = m2d(x(1));

    double alpha_ref = M_PI;
    double omega_ref = 0.0;

    matrix ref(2, 1);
    ref(0) = alpha_ref;
    ref(1) = omega_ref;

    matrix Y0(2, 1);
    Y0(0) = 0.0;
    Y0(1) = 0.0;

    matrix* Y = solve_ode(df2, 0.0, 0.1, 100.0, Y0, x, ref);
    
    int n = get_len(Y[0]);
    
    double Q = 0.0;
    double dt = 0.1;

    for (int i = 0; i < n; i++)
    {
        double alpha = Y[1](i, 0);
        double omega = Y[1](i, 1);

        double M = k1 * (alpha_ref - alpha) + k2 * (omega_ref - omega);

		cout << Y[0](i, 0) << "," << alpha << "," << omega << endl;

        Q += (10.0 * pow(alpha_ref - alpha, 2) + pow(omega_ref - omega, 2) + pow(M, 2)) * dt;
    }
    
    y = Q;

    return y;
}

matrix ff3T(matrix x, matrix ud1, matrix ud2) {
    double x1 = m2d(x(0));
    double x2 = m2d(x(1));
    double common = M_PI * sqrt(pow(x1/M_PI, 2) + pow(x2/M_PI, 2));
    double t1 = sin(common);

    return matrix(t1 / common);
}

matrix ff3T_zewn(matrix x, matrix ud1, matrix ud2) {
    double x1 = m2d(x(0));
    double x2 = m2d(x(1));
    double a = m2d(ud1(0));
    double c = m2d(ud1(1));

    double g1 = -x1 + 1.0;
    double g2 = -x2 + 1.0;
    double g3 = sqrt(x1*x1 + x2*x2) - a;

    double S = pow(fmax(0.0, g1), 2) + pow(fmax(0.0, g2), 2) + pow(fmax(0.0, g3), 2);

    double common = M_PI * sqrt(pow(x1/M_PI, 2) + pow(x2/M_PI, 2));
    double f = sin(common) / common;

    return matrix(f + c * S);
}

matrix ff3T_wewn(matrix x, matrix ud1, matrix ud2) {
    double x1 = m2d(x(0));
    double x2 = m2d(x(1));
    double a = m2d(ud1(0));
    double c = m2d(ud1(1));

    double g1 = -x1 + 1.0;
    double g2 = -x2 + 1.0;
    double g3 = sqrt(x1*x1 + x2*x2) - a;

    if (g1 >= 0 || g2 >= 0 || g3 >= 0) {
        return matrix(1e10);  
    }

    double S = -1.0/g1 - 1.0/g2 - 1.0/g3;

    double common = M_PI * sqrt(pow(x1/M_PI, 2) + pow(x2/M_PI, 2));
    double f = sin(common) / common;

    return matrix(f + c * S);
}

matrix df3(double t, matrix Y, matrix ud1, matrix ud2) {
    matrix dY(4, 1);
    
    double m = 0.6;     
    double r = 0.12;    
    double g = 9.81;    
    double C = 0.47;    
    double rho = 1.2;   
    double S = M_PI * r * r; 
    
    double omega = m2d(ud1(1));  
    
    double vx = Y(2);
    double vy = Y(3);
    
    double Dx = 0.5 * C * rho * S * vx * fabs(vx);
    double Dy = 0.5 * C * rho * S * vy * fabs(vy);
    
    double FMx = rho * vy * omega * M_PI * pow(r, 3);
    double FMy = rho * vx * omega * M_PI * pow(r, 3);
    
    dY(0) = vx;                             
    dY(1) = vy;                             
    dY(2) = (-Dx - FMx) / m;                
    dY(3) = (-Dy - FMy - m * g) / m;        
    
    return dY;
}

matrix ff3R(matrix x, matrix ud1, matrix ud2) {
    double v0x = m2d(x(0));
    double omega = m2d(x(1));
    double c = m2d(ud1(1)); 
    

    double g_v0x_min = -10.0 - v0x;
    double g_v0x_max = v0x - 10.0;
    double g_omega_min = -10.0 - omega;
    double g_omega_max = omega - 10.0;
    
    double S_bounds = pow(fmax(0.0, g_v0x_min), 2) + pow(fmax(0.0, g_v0x_max), 2) +
                      pow(fmax(0.0, g_omega_min), 2) + pow(fmax(0.0, g_omega_max), 2);
    
    matrix Y0(4, 1);
    Y0(0) = 0.0;    
    Y0(1) = 100.0; 
    Y0(2) = v0x;   
    Y0(3) = 0.0;    
    
    matrix params(2, 1);
    params(0) = v0x;
    params(1) = omega;
    
    matrix* Y = solve_ode(df3, 0.0, 0.01, 7.0, Y0, params, NAN);
    
    int n = get_len(Y[0]);
    
    double x_end = 0.0;
    double x_at_y50 = 0.0;
    bool found_y50 = false;
    
    for (int i = 0; i < n; ++i) {
        double yi = Y[1](i, 1);
        double xi = Y[1](i, 0);
        
        if (!found_y50 && yi <= 50.0) {
            x_at_y50 = xi;
            found_y50 = true;
        }
        
        if (yi <= 0.0) {
            x_end = xi;
            break;
        }
        
        if (i == n - 1) {
            x_end = xi;
        }
    }

    double g1 = 3.0 - x_at_y50;
    double g2 = x_at_y50 - 7.0;
    
    double S = pow(fmax(0.0, g1), 2) + pow(fmax(0.0, g2), 2) + S_bounds;
    
    double f = -x_end + c * S;
    
    Y[0].~matrix();
    Y[1].~matrix();
    
    return matrix(f);
}

// ====== LAB 4 FUNCTIONS ======

// Test function: f(x1,x2) = (1/6)x1^6 - 1.05*x1^4 + 2*x1^2 + x2^2 + x1*x2
matrix ff4T(matrix x, matrix ud1, matrix ud2) {
    double x1 = m2d(x(0));
    double x2 = m2d(x(1));
    
    double f = (1.0/6.0) * pow(x1, 6) - 1.05 * pow(x1, 4) + 2.0 * pow(x1, 2) + pow(x2, 2) + x1 * x2;
    
    return matrix(f);
}

// Gradient of test function
matrix gf4T(matrix x, matrix ud1, matrix ud2) {
    double x1 = m2d(x(0));
    double x2 = m2d(x(1));
    
    matrix grad(2, 1);
    // df/dx1 = x1^5 - 4.2*x1^3 + 4*x1 + x2
    grad(0) = pow(x1, 5) - 4.2 * pow(x1, 3) + 4.0 * x1 + x2;
    // df/dx2 = 2*x2 + x1
    grad(1) = 2.0 * x2 + x1;
    
    return grad;
}

// Hessian of test function
matrix Hf4T(matrix x, matrix ud1, matrix ud2) {
    double x1 = m2d(x(0));
    double x2 = m2d(x(1));
    
    matrix H(2, 2);
    // d2f/dx1dx1 = 5*x1^4 - 12.6*x1^2 + 4
    H(0, 0) = 5.0 * pow(x1, 4) - 12.6 * pow(x1, 2) + 4.0;
    // d2f/dx1dx2 = 1
    H(0, 1) = 1.0;
    // d2f/dx2dx1 = 1
    H(1, 0) = 1.0;
    // d2f/dx2dx2 = 2
    H(1, 1) = 2.0;
    
    return H;
}

// Logistic regression cost function J(theta)
// ud1 should contain X data (m x 3 matrix with bias column)
// ud2 should contain Y data (m x 1 matrix)
matrix ff4R_cost(matrix theta, matrix ud1, matrix ud2) {
    int* size = get_size(ud1);
    int m = size[0]; // number of samples
    delete[] size;
    
    const double EPSILON = 1e-15; // Small value to avoid log(0)
    
    double J = 0.0;
    for (int i = 0; i < m; i++) {
        // Compute h = 1 / (1 + exp(-theta^T * x^(i)))
        double z = 0.0;
        for (int j = 0; j < 3; j++) {
            z += m2d(theta(j)) * m2d(ud1(i, j));
        }
        double h = 1.0 / (1.0 + exp(-z));
        
        double y = m2d(ud2(i));
        
        // J += y * log(h) + (1-y) * log(1-h)
        h = fmax(EPSILON, fmin(1.0 - EPSILON, h));
        J += y * log(h) + (1.0 - y) * log(1.0 - h);
    }
    
    J = -J / m;
    
    return matrix(J);
}

// Gradient of logistic regression cost function
// ud1 should contain X data (m x 3 matrix with bias column)
// ud2 should contain Y data (m x 1 matrix)
matrix gf4R_grad(matrix theta, matrix ud1, matrix ud2) {
    int* size = get_size(ud1);
    int m = size[0]; // number of samples
    delete[] size;
    
    matrix grad(3, 1);
    grad(0) = 0.0;
    grad(1) = 0.0;
    grad(2) = 0.0;
    
    for (int i = 0; i < m; i++) {
        // Compute h = 1 / (1 + exp(-theta^T * x^(i)))
        double z = 0.0;
        for (int j = 0; j < 3; j++) {
            z += m2d(theta(j)) * m2d(ud1(i, j));
        }
        double h = 1.0 / (1.0 + exp(-z));
        
        double y = m2d(ud2(i));
        double diff = h - y;
        
        // grad(j) += (h - y) * x_j^(i)
        for (int j = 0; j < 3; j++) {
            grad(j) = grad(j) + diff * m2d(ud1(i, j));
        }
    }
    
    // Average over all samples
    for (int j = 0; j < 3; j++) {
        grad(j) = grad(j) / m;
    }
    
    return grad;
}

// ====== LAB 5 FUNCTIONS ======

// First criterion for test problem: f1(x) = a*((x1-3)^2 + (x2-3)^2)
// ud1(0) = parameter a
matrix ff5_f1(matrix x, matrix ud1, matrix ud2) {
    double x1 = m2d(x(0));
    double x2 = m2d(x(1));
    double a = m2d(ud1(0));
    
    double f1 = a * (pow(x1 - 3.0, 2) + pow(x2 - 3.0, 2));
    
    return matrix(f1);
}

// Second criterion for test problem: f2(x) = (1/a)*((x1+3)^2 + (x2+3)^2)
// ud1(0) = parameter a
matrix ff5_f2(matrix x, matrix ud1, matrix ud2) {
    double x1 = m2d(x(0));
    double x2 = m2d(x(1));
    double a = m2d(ud1(0));
    
    double f2 = (1.0 / a) * (pow(x1 + 3.0, 2) + pow(x2 + 3.0, 2));
    
    return matrix(f2);
}

// Weighted objective for test problem: f(x) = w*f1(x) + (1-w)*f2(x)
// ud1(0) = parameter a
// ud1(1) = weight w
matrix ff5T(matrix x, matrix ud1, matrix ud2) {
    double w = m2d(ud1(1));
    
    matrix f1 = ff5_f1(x, ud1, ud2);
    matrix f2 = ff5_f2(x, ud1, ud2);
    
    double f = w * m2d(f1) + (1.0 - w) * m2d(f2);
    
    return matrix(f);
}

// Weighted objective for real problem: cantilever beam
// x(0) = diameter d [mm]
// x(1) = length l [mm]
// ud1(0) = weight w
// Returns weighted objective with penalties for constraint violations
matrix ff5R(matrix x, matrix ud1, matrix ud2) {
    double d_orig = m2d(x(0));  // diameter in mm
    double l_orig = m2d(x(1));  // length in mm
    double w = m2d(ud1(0)); // weight
    
    // Constants
    const double P = 3000.0;        // force in N (3 kN)
    const double E = 120e9;         // Young's modulus in Pa (120 GPa)
    const double rho = 8920.0;      // density in kg/m^3
    const double u_max = 2.5;       // max deflection in mm
    const double sigma_max = 300e6; // max stress in Pa (300 MPa)
    
    const double penalty_coef = 1e5;  // Reduced from 1e10 to avoid overflow
    
    // Penalty for bounds: d in [0.01, 1000], l in [0.2, 1000]
    double f = 0.0;
    
    // Apply penalties for out-of-bounds, but clamp for calculation
    double d = d_orig;
    double l = l_orig;
    
    if (d_orig < 0.01) {
        f += penalty_coef * pow(0.01 - d_orig, 2);
        d = 0.01;
    }
    if (d_orig > 1000.0) {
        f += penalty_coef * pow(d_orig - 1000.0, 2);
        d = 1000.0;
    }
    if (l_orig < 0.2) {
        f += penalty_coef * pow(0.2 - l_orig, 2);
        l = 0.2;
    }
    if (l_orig > 1000.0) {
        f += penalty_coef * pow(l_orig - 1000.0, 2);
        l = 1000.0;
    }
    
    // Convert to SI units (meters)
    double d_m = d / 1000.0;
    double l_m = l / 1000.0;
    
    // Calculate mass (kg)
    double mass = rho * M_PI * pow(d_m / 2.0, 2) * l_m;
    
    // Calculate deflection (convert to mm)
    double u = (64.0 * P * pow(l_m, 3)) / (3.0 * E * M_PI * pow(d_m, 4));
    u = u * 1000.0;  // convert to mm
    
    // Calculate stress (Pa)
    double sigma = (32.0 * P * l_m) / (M_PI * pow(d_m, 3));
    
    // Criteria
    double f1 = mass;       // minimize mass
    double f2 = u;          // minimize deflection
    
    // Weighted objective
    f += w * f1 + (1.0 - w) * f2;
    
    // Add penalties for constraint violations
    const double constraint_penalty = 1e4;  // Reduced from 1e6
    
    // Penalty for deflection constraint: u <= u_max
    if (u > u_max) {
        f += constraint_penalty * pow(u - u_max, 2);
    }
    
    // Penalty for stress constraint: sigma <= sigma_max
    if (sigma > sigma_max) {
        f += constraint_penalty * pow((sigma - sigma_max) / 1e6, 2);  // normalize to MPa scale
    }
    
    return matrix(f);
}