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