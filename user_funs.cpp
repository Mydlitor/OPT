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

// Zewnętrzna funkcja kary dla problemu testowego K3
// ud1(0) = a (parametr ograniczenia)
// ud1(1) = c (współczynnik kary)
matrix ff3T_zewn(matrix x, matrix ud1, matrix ud2) {
    double x1 = m2d(x(0));
    double x2 = m2d(x(1));
    double a = m2d(ud1(0));
    double c = m2d(ud1(1));

    // Ograniczenia:
    // g1 = -x1 + 1 <= 0  (x1 >= 1)
    // g2 = -x2 + 1 <= 0  (x2 >= 1)
    // g3 = sqrt(x1^2 + x2^2) - a <= 0
    double g1 = -x1 + 1.0;
    double g2 = -x2 + 1.0;
    double g3 = sqrt(x1*x1 + x2*x2) - a;

    // Funkcja kary zewnętrznej: S = sum(max(0, gi)^2)
    double S = pow(fmax(0.0, g1), 2) + pow(fmax(0.0, g2), 2) + pow(fmax(0.0, g3), 2);

    // Funkcja celu (ff3T)
    double common = M_PI * sqrt(pow(x1/M_PI, 2) + pow(x2/M_PI, 2));
    double f = sin(common) / common;

    return matrix(f + c * S);
}

// Wewnętrzna funkcja kary dla problemu testowego K3
// ud1(0) = a (parametr ograniczenia)
// ud1(1) = c (współczynnik kary)
matrix ff3T_wewn(matrix x, matrix ud1, matrix ud2) {
    double x1 = m2d(x(0));
    double x2 = m2d(x(1));
    double a = m2d(ud1(0));
    double c = m2d(ud1(1));

    // Ograniczenia (muszą być < 0 dla punktu wewnątrz):
    double g1 = -x1 + 1.0;
    double g2 = -x2 + 1.0;
    double g3 = sqrt(x1*x1 + x2*x2) - a;

    // Sprawdzenie czy punkt jest wewnątrz obszaru dopuszczalnego
    if (g1 >= 0 || g2 >= 0 || g3 >= 0) {
        return matrix(1e10);  // Bardzo duża wartość dla punktu poza obszarem
    }

    // Funkcja kary wewnętrznej: S = -sum(1/gi)
    double S = -1.0/g1 - 1.0/g2 - 1.0/g3;

    // Funkcja celu (ff3T)
    double common = M_PI * sqrt(pow(x1/M_PI, 2) + pow(x2/M_PI, 2));
    double f = sin(common) / common;

    return matrix(f + c * S);
}

// Równania ruchu piłki dla problemu rzeczywistego K3
// Y = [x, y, vx, vy]
// ud1 = [v0x, omega]
matrix df3(double t, matrix Y, matrix ud1, matrix ud2) {
    matrix dY(4, 1);
    
    double m = 0.6;       // masa piłki [kg]
    double r = 0.12;      // promień piłki [m]
    double g = 9.81;      // przyspieszenie grawitacyjne [m/s^2]
    double C = 0.47;      // współczynnik oporu
    double rho = 1.2;     // gęstość powietrza [kg/m^3]
    double S = M_PI * r * r;  // pole przekroju
    
    double omega = m2d(ud1(1));  // rotacja [rad/s]
    
    double vx = Y(2);
    double vy = Y(3);
    
    // Siła oporu powietrza
    double Dx = 0.5 * C * rho * S * vx * fabs(vx);
    double Dy = 0.5 * C * rho * S * vy * fabs(vy);
    
    // Siła Magnusa
    double FMx = rho * vy * omega * M_PI * pow(r, 3);
    double FMy = rho * vx * omega * M_PI * pow(r, 3);
    
    // Równania ruchu
    dY(0) = vx;                                    // dx/dt = vx
    dY(1) = vy;                                    // dy/dt = vy
    dY(2) = (-Dx - FMx) / m;                       // dvx/dt
    dY(3) = (-Dy - FMy - m * g) / m;               // dvy/dt = (-Dy - FMy)/m - g
    
    return dY;
}

// Funkcja celu dla problemu rzeczywistego K3 (spadająca piłka)
// x = [v0x, omega]
// ud1(1) = c (współczynnik kary) - zgodnie z konwencją pen()
matrix ff3R(matrix x, matrix ud1, matrix ud2) {
    double v0x = m2d(x(0));
    double omega = m2d(x(1));
    double c = m2d(ud1(1));  // współczynnik kary jest w ud1(1)
    
    // Ograniczenia na zmienne: v0x in [-10, 10], omega in [-10, 10]
    // g_v0x_min = -10 - v0x <= 0  (v0x >= -10)
    // g_v0x_max = v0x - 10 <= 0   (v0x <= 10)
    // g_omega_min = -10 - omega <= 0  (omega >= -10)
    // g_omega_max = omega - 10 <= 0   (omega <= 10)
    double g_v0x_min = -10.0 - v0x;
    double g_v0x_max = v0x - 10.0;
    double g_omega_min = -10.0 - omega;
    double g_omega_max = omega - 10.0;
    
    double S_bounds = pow(fmax(0.0, g_v0x_min), 2) + pow(fmax(0.0, g_v0x_max), 2) +
                      pow(fmax(0.0, g_omega_min), 2) + pow(fmax(0.0, g_omega_max), 2);
    
    // Warunki początkowe: [x0, y0, vx0, vy0]
    matrix Y0(4, 1);
    Y0(0) = 0.0;      // x0 = 0
    Y0(1) = 100.0;    // y0 = 100m
    Y0(2) = v0x;      // vx0 = v0x
    Y0(3) = 0.0;      // vy0 = 0
    
    matrix params(2, 1);
    params(0) = v0x;
    params(1) = omega;
    
    // Symulacja: t0=0, dt=0.01, tend=7
    matrix* Y = solve_ode(df3, 0.0, 0.01, 7.0, Y0, params, NAN);
    
    int n = get_len(Y[0]);
    
    // Szukamy x_end (położenie x gdy y <= 0)
    double x_end = 0.0;
    double x_at_y50 = 0.0;
    bool found_y50 = false;
    
    for (int i = 0; i < n; ++i) {
        double yi = Y[1](i, 1);
        double xi = Y[1](i, 0);
        
        // Szukamy x gdy y jest blisko 50m
        if (!found_y50 && yi <= 50.0) {
            x_at_y50 = xi;
            found_y50 = true;
        }
        
        // Szukamy x_end gdy y <= 0
        if (yi <= 0.0) {
            x_end = xi;
            break;
        }
        
        // Jeśli doszliśmy do końca symulacji
        if (i == n - 1) {
            x_end = xi;
        }
    }
    
    // Ograniczenie: dla y=50, x musi być w [3, 7]
    // g1 = 3 - x_at_y50 <= 0  (x >= 3)
    // g2 = x_at_y50 - 7 <= 0  (x <= 7)
    double g1 = 3.0 - x_at_y50;
    double g2 = x_at_y50 - 7.0;
    
    // Funkcja kary zewnętrznej (ograniczenie na przejście przez punkt + ograniczenia brzegowe)
    double S = pow(fmax(0.0, g1), 2) + pow(fmax(0.0, g2), 2) + S_bounds;
    
    // Funkcja celu: maksymalizacja x_end => minimalizacja -x_end
    double f = -x_end + c * S;
    
    Y[0].~matrix();
    Y[1].~matrix();
    
    return matrix(f);
}