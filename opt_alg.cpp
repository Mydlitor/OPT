#include "opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	// Zmienne wejściowe:
	// ff - wskaźnik do funkcji celu
	// N - liczba zmiennych funkcji celu
	// lb, ub - dolne i górne ograniczenie
	// epslion - zakłądana dokładność rozwiązania
	// Nmax - maksymalna liczba wywołań funkcji celu
	// ud1, ud2 - user data
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);									// losujemy macierz Nx1 stosując rozkład jednostajny na przedziale [0,1]
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);// przeskalowywujemy rozwiązanie do przedziału [lb, ub]
			Xopt.fit_fun(ff, ud1, ud2);							// obliczmy wartość funkcji celu
			if (Xopt.y < epsilon)								// sprawdzmy 1. kryterium stopu
			{
				Xopt.flag = 1;									// flaga = 1 ozancza znalezienie rozwiązanie z zadaną dokładnością
				break;
			}
			if (solution::f_calls > Nmax)						// sprawdzmy 2. kryterium stopu
			{
				Xopt.flag = 0;									// flaga = 0 ozancza przekroczenie maksymalne liczby wywołań funkcji celu
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2]{ 0, 0 };

		int i = 0;

		solution X0(x0);
		solution X1(x0 + d);

		X0.fit_fun(ff, ud1, ud2);
		X1.fit_fun(ff, ud1, ud2);

		if (X1.y == X0.y)
		{
			p[0] = X0.x(0);
			p[1] = X1.x(0);
			return p;
		}

		if (X1.y > X0.y)
		{
			d = -d;
			X1.x = x0 + d;
			X1.fit_fun(ff, ud1, ud2);

			if (X1.y >= X0.y)
			{
				p[0] = X1.x(0);
				p[1] = X0.x(0) - d;
				return p;
			}
		}

		solution Xi_prev = X0;
		solution Xi = X1;
		solution Xi_next;

		do
		{
			if (solution::f_calls >= Nmax)
			{
				throw string("Przekroczono maksymalna liczbe wywolan funkcji celu");
			}

			i++;
			Xi_next.x = x0 + pow(alpha, i) * d;
			Xi_next.fit_fun(ff, ud1, ud2);

			if (Xi.y <= Xi_next.y)
			{
				break;
			}

			Xi_prev = Xi;
			Xi = Xi_next;

		} while (true);

		if (d > 0)
		{
			p[0] = Xi_prev.x(0);
			p[1] = Xi_next.x(0);
		}
		else
		{
			p[0] = Xi_next.x(0);
			p[1] = Xi_prev.x(0);
		}

		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;

		double phi = (1.0 + sqrt(5.0)) / 2.0;
		double r_phi = 1.0 / phi;

		int k = 0;
		double phi_k = 1.0;

		while (phi_k <= (b - a) / epsilon)
		{
			phi_k *= phi;
			k++;
		}

		solution A(a), B(b), C, D;

		C.x = B.x - r_phi * (B.x - A.x);
		D.x = A.x + B.x - C.x;

		C.fit_fun(ff, ud1, ud2);
		D.fit_fun(ff, ud1, ud2);

		for (int i = 0; i <= k - 3; i++)
		{
			if (C.y < D.y)
			{
				B = D;
			}
			else
			{
				A = C;
			}

			C.x = B.x - r_phi * (B.x - A.x);
			D.x = A.x + B.x - C.x;

			C.fit_fun(ff, ud1, ud2);
			D.fit_fun(ff, ud1, ud2);

			//debug - wielkosc przedzialu w zaleznosci od iteracji
			//cout << fabs(C.x(0) - D.x(0)) << "\n";
		}

		Xopt = C;

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}
}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		//double* interval = expansion(ff, a, 2, 1.01, Nmax, ud1, ud2);
		double interval[2] = {a,b};
		solution Xopt;
		solution::clear_calls();

		double aa = interval[0];
		double bb = interval[1];
		double cc = (interval[1] + interval[0]) / 2;
		double dd = 0, dd_prev = 0;
		int iter = 0;

		do {
			// Oblicz wartości funkcji celu
			Xopt.x = aa; Xopt.fit_fun(ff, ud1, ud2); double f_a = m2d(Xopt.y);
			Xopt.x = bb; Xopt.fit_fun(ff, ud1, ud2); double f_b = m2d(Xopt.y);
			Xopt.x = cc; Xopt.fit_fun(ff, ud1, ud2); double f_c = m2d(Xopt.y);

			// Interpolacja kwadratowa
			double l = f_a * (pow(bb, 2) - pow(cc, 2)) + f_b * (pow(cc, 2) - pow(aa, 2)) + f_c * (pow(aa, 2) - pow(bb, 2));
			double m = f_a * (bb - cc) + f_b * (cc - aa) + f_c * (aa - bb);
			if (m <= 0) throw std::string("M <= 0");

			dd_prev = dd;
			dd = 0.5 * l / m;

			// Zabezpieczenie przed wyjściem poza przedział
			if (dd < aa || dd > bb) {
				//std::cout << "aa: " << aa << ", bb: " << bb << ", cc: " << cc << ", dd: " << dd << std::endl;
				throw std::string("dd poza przedziałem!");
			}

			Xopt.x = dd; Xopt.fit_fun(ff, ud1, ud2); double f_d = m2d(Xopt.y);

			// Warunki wyboru nowego przedziału
			if ((aa < dd) && (dd < cc)) {
				if (f_d < f_c) {
					bb = cc;
					cc = dd;
				}
				else {
					aa = dd;
				}
			}
			else if ((cc < dd) && (dd < bb)) {
				if (f_d < f_c) {
					aa = cc;
					cc = dd;
				}
				else {
					bb = dd;
				}
			}
			else {
				// Zabezpieczenie przed zdegenerowaniem przedziału
				if (fabs(bb - aa) < epsilon || fabs(dd - dd_prev) < gamma) break;
				throw std::string("Error 4321");
			}

			iter++;
			if (solution::f_calls > Nmax) throw std::string("Too many iterations");

			// Warunek zbieżności
			if (fabs(bb - aa) < epsilon || fabs(dd - dd_prev) < gamma) break;

			//debug - wielkosc przedzialu w zaleznosci od iteracji
			//cout << fabs(aa-bb) << "\n";

		} while (true);

		Xopt.x = dd;
		Xopt.fit_fun(ff, ud1, ud2);
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

matrix sol2mat(solution sol) {
	matrix ret(2,1);
	ret(0) = m2d(sol.x);
	ret(1) = m2d(sol.y);
	return ret;
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		matrix x = x0;
		matrix X_b_prev = x;
		matrix X_b = x;

		solution Xopt;
		
		double f_x, f_xb; //f(x), f(x_b)

		while(s >= epsilon) {
			X_b = x;
			x = sol2mat(HJ_trial(ff, X_b, s)); //Convert solution.x, y to matrix (x, y), function ff takes matrix as input.
			Xopt.x = x; Xopt.fit_fun(ff); f_x = m2d(Xopt.y);
			Xopt.x = X_b; Xopt.fit_fun(ff); f_xb = m2d(Xopt.y);
			if (f_x < f_xb) {
				while (f_x < f_xb) {
					X_b_prev = X_b;
					X_b = x;
					x = 2 * X_b - X_b_prev;
					x = sol2mat(HJ_trial(ff, x, s));
					
					Xopt.x = x; Xopt.fit_fun(ff); f_x = m2d(Xopt.y);
					Xopt.x = X_b; Xopt.fit_fun(ff); f_xb = m2d(Xopt.y);

					if (solution::f_calls > Nmax) {
						throw std::string("Dupas zbitas");
					}
				}
				x = X_b;
			} else {
				s = alpha * s;
			}
			if (solution::f_calls > Nmax) {
				throw std::string("Dupas zbitas 2");
			}
		}
		Xopt.x = X_b;
		Xopt.fit_fun(ff);
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), matrix XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		matrix ej = ident_mat(2); //spontanicznie przypadkowo maciez jednostkowa pasuje tutaj

		solution Xopt;
		double f_x, f_forward, f_backward;
		
		for (int i = 0; i < 2; i++) {
			Xopt.x = XB; 
			Xopt.fit_fun(ff);
			f_x = m2d(Xopt.y);

			Xopt.x = XB + s * ej[i]; 
			Xopt.fit_fun(ff);
			 f_forward = m2d(Xopt.y);


			Xopt.x = XB - s * ej[i]; Xopt.fit_fun(ff); f_backward = m2d(Xopt.y);
			
			if (f_forward < f_x) {
				XB = XB + s*ej[i];
			} else if (f_backward < f_x) {
				XB = XB - s*ej[i];
			}
		}
		Xopt.x = XB(0);
		Xopt.y = XB(1);
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
