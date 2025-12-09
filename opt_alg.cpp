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
	int iteracja = 0;
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
					cout << iteracja << "," << Xopt.x(0) << "," << Xopt.x(1) << endl;
					iteracja++;
					X_b_prev = X_b;
					X_b = x;
					x = 2 * X_b - X_b_prev;
					x = sol2mat(HJ_trial(ff, x, s));
					
					Xopt.x = x; Xopt.fit_fun(ff); f_x = m2d(Xopt.y);
					Xopt.x = X_b; Xopt.fit_fun(ff); f_xb = m2d(Xopt.y);

					if (solution::f_calls >= Nmax)
					{
						throw string("Przekroczono maksymalna liczbe wywolan funkcji celu");
					}
				}
				x = X_b;
			} else {
				s = alpha * s;
				cout << iteracja << "," << Xopt.x(0) << "," << Xopt.x(1) << endl;
				iteracja++;
			}
			if (solution::f_calls >= Nmax)
			{
				throw string("Przekroczono maksymalna liczbe wywolan funkcji celu");
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

solution Rosen(matrix (*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	int iteracja = 0;
	try
	{
		solution Xopt;

		int n = get_len(x0);
		int i = 0;

		matrix D(n, n);
		for (int j = 0; j < n; j++)
		{
			for (int k = 0; k < n; k++)
			{
				D(k, j) = (j == k) ? 1.0 : 0.0;
			}
		}

		matrix s = s0;
		matrix lambda(n, 1);
		matrix p(n, 1);

		for (int j = 0; j < n; j++)
		{
			lambda(j) = 0.0;
			p(j) = 0.0;
		}

		solution xB(x0);
		xB.fit_fun(ff, ud1, ud2);

		while (true)
		{
			cout << iteracja << "," << xB.x(0) << "," << xB.x(1) << endl;
			iteracja++;
			for (int j = 0; j < n; j++)
			{
				matrix d_j(n, 1);
				for (int k = 0; k < n; k++)
				{
					d_j(k) = D(k, j);
				}

				matrix x_new = xB.x + s(j) * d_j;
				solution xNew(x_new);
				xNew.fit_fun(ff, ud1, ud2);

				if (xNew.y(0) < xB.y(0))
				{
					xB.x = xNew.x;
					xB.y = xNew.y;
					lambda(j) = lambda(j) + s(j);
					s(j) = alpha * s(j);
				}
				else
				{
					s(j) = -beta * s(j);
					p(j) = p(j) + 1;
				}
			}

			i++;

			bool all_lambda_nonzero = true;
			bool all_p_nonzero = true;

			for (int j = 0; j < n; j++)
			{
				if (lambda(j) == 0.0)
					all_lambda_nonzero = false;
				if (p(j) == 0.0)
					all_p_nonzero = false;
			}

			if (all_lambda_nonzero && all_p_nonzero)
			{
				matrix Q(n, n);
				for (int row = 0; row < n; row++)
				{
					for (int col = 0; col < n; col++)
					{
						if (row >= col)
						{
							Q(row, col) = lambda(col);
						}
						else
						{
							Q(row, col) = 0.0;
						}
					}
				}

				matrix Q_star = D * Q;

				for (int j = 0; j < n; j++)
				{
					matrix v_j(n, 1);
					for (int k = 0; k < n; k++)
					{
						v_j(k) = Q_star(k, j);
					}

					for (int k = 0; k < j; k++)
					{
						matrix d_k(n, 1);
						for (int l = 0; l < n; l++)
						{
							d_k(l) = D(l, k);
						}

						double dot_product = 0.0;
						for (int l = 0; l < n; l++)
						{
							dot_product += v_j(l) * d_k(l);
						}

						for (int l = 0; l < n; l++)
						{
							v_j(l) -= dot_product * d_k(l);
						}
					}

					double norm_val = 0.0;
					for (int k = 0; k < n; k++)
					{
						norm_val += v_j(k) * v_j(k);
					}
					norm_val = sqrt(norm_val);

					if (norm_val > 1e-10)
					{
						for (int k = 0; k < n; k++)
						{
							D(k, j) = v_j(k) / norm_val;
						}
					}
				}

				for (int j = 0; j < n; j++)
				{
					lambda(j) = 0.0;
					p(j) = 0.0;
					s(j) = s0(j);
				}
			}

			if (solution::f_calls > Nmax)
			{
				Xopt = xB;
				Xopt.flag = 0;
				return Xopt;
			}

			double max_s = 0.0;
			for (int j = 0; j < n; j++)
			{
				double abs_s = fabs(s(j));
				if (abs_s > max_s)
					max_s = abs_s;
			}

			if (max_s < epsilon)
			{
				break;
			}
		}

		Xopt = xB;
		Xopt.flag = 1;

		return Xopt;
	}
	catch (string ex_info)
	{
		throw("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		
		// Parametry dla sym_NM
		double alpha = 1.0;
		double beta = 0.5;
		double gamma = 2.0;
		double delta = 0.5;
		double s = 0.5;
		double eps_NM = 1e-2;
		
		matrix x_prev = x0;
		matrix x_curr = x0;
		
		matrix params(2, 1);
		params(0) = ud1(0); 
		params(1) = c;       
		
		int iter = 0;
		
		do {
			x_prev = x_curr;
			
			params(1) = c;
			
			Xopt = sym_NM(ff, x_curr, s, alpha, beta, gamma, delta, eps_NM, Nmax, params, ud2);
			
			x_curr = Xopt.x;
			
			c = dc * c;
			
			iter++;
			
			if (solution::f_calls > Nmax) {
				Xopt.flag = 0;
				return Xopt;
			}
			
			double dist = 0.0;
			int len = get_len(x_curr);
			for (int i = 0; i < len; ++i) {
				dist += pow(x_curr(i) - x_prev(i), 2);
			}
			dist = sqrt(dist);
			
			if (dist < epsilon) {
				break;
			}
			
		} while (true);
		
		Xopt.flag = 1;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

bool odleglosc_wieksza_od_epsilon(matrix* p, int n, int min_idx, double epsilon) {
	for (int i = 0; i <= n; i++) {
		if (i == min_idx) continue;
		double dist = 0.0;
		for (int j = 0; j < n; j++) {
			dist += pow(p[min_idx](j) - p[i](j), 2);
		}
		if (sqrt(dist) >= epsilon) return true;
	}
	return false;
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		int n = get_len(x0);
		
		matrix* p = new matrix[n + 1];
		p[0] = x0;
		for (int i = 1; i <= n; i++) {
			p[i] = x0 + s * ident_mat(n)[i-1];
		}
		
		double* f = new double[n + 1]; 
		
		do {
			for (int i = 0; i <= n; i++) {
				Xopt.x = p[i]; 
				Xopt.fit_fun(ff, ud1, ud2); 
				f[i] = m2d(Xopt.y);
			}
			
			int min_idx = 0, max_idx = 0;
			for (int i = 1; i <= n; i++) {
				if (f[i] < f[min_idx]) min_idx = i;
				if (f[i] > f[max_idx]) max_idx = i;
			}
			
			matrix p_bar(n, 1, 0.0);
			for (int i = 0; i <= n; i++) {
				if (i != max_idx) {
					p_bar = p_bar + p[i];
				}
			}
			p_bar = p_bar / n;
			
			matrix p_odb = p_bar + alpha * (p_bar - p[max_idx]);
			Xopt.x = p_odb; Xopt.fit_fun(ff, ud1, ud2);
			double f_odb = m2d(Xopt.y);
			
			if (f_odb < f[min_idx]) {
				matrix p_e = p_bar + gamma * (p_odb - p_bar);
				Xopt.x = p_e; Xopt.fit_fun(ff, ud1, ud2);
				double f_e = m2d(Xopt.y);
				
				if (f_e < f_odb) {
					p[max_idx] = p_e;
				} else {
					p[max_idx] = p_odb;
				}
			} else {
				if (f[min_idx] <= f_odb && f_odb < f[max_idx]) {
					p[max_idx] = p_odb;
				} else {
					matrix p_z = p_bar + beta * (p[max_idx] - p_bar);
					Xopt.x = p_z; Xopt.fit_fun(ff, ud1, ud2);
					double f_z = m2d(Xopt.y);
					
					if (f_z >= f[max_idx]) {
						for (int i = 0; i <= n; i++) {
							if (i != min_idx) {
								p[i] = delta * (p[i] + p[min_idx]);
							}
						}
					} else {
						p[max_idx] = p_z;
					}
				}
			}
			
			if (solution::f_calls > Nmax) {
				Xopt.flag = 0;
				break;
			}
			
			for (int i = 0; i <= n; i++) {
				Xopt.x = p[i]; 
				Xopt.fit_fun(ff, ud1, ud2); 
				f[i] = m2d(Xopt.y);
			}
			min_idx = 0;
			for (int i = 1; i <= n; i++) {
				if (f[i] < f[min_idx]) min_idx = i;
			}
			
			if (!odleglosc_wieksza_od_epsilon(p, n, min_idx, epsilon)) {
				break;
			}
			
		} while (true);
		
		int best_idx = 0;
		for (int i = 1; i <= n; i++) {
			Xopt.x = p[i]; Xopt.fit_fun(ff, ud1, ud2);
			double fi = m2d(Xopt.y);
			Xopt.x = p[best_idx]; Xopt.fit_fun(ff, ud1, ud2);
			double fb = m2d(Xopt.y);
			if (fi < fb) best_idx = i;
		}
		
		Xopt.x = p[best_idx];
		Xopt.fit_fun(ff, ud1, ud2);
		Xopt.flag = 1;
		
		delete[] p;
		delete[] f;
		
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

// Helper function for golden section line search along a direction
// Returns optimal step size alpha for x_new = x_current + alpha * direction
double golden_line_search(matrix(*ff)(matrix, matrix, matrix), matrix x_current, matrix direction, 
                          double a, double b, double tol, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double alpha_gold = (sqrt(5.0) - 1.0) / 2.0;
		
		double c = b - alpha_gold * (b - a);
		double d_val = a + alpha_gold * (b - a);
		
		matrix x_c = x_current + c * direction;
		matrix x_d = x_current + d_val * direction;
		double fc = m2d(ff(x_c, ud1, ud2));
		double fd = m2d(ff(x_d, ud1, ud2));
		
		int ls_iters = 0;
		while (b - a > tol && solution::f_calls < Nmax && ls_iters < 100) {
			if (fc < fd) {
				b = d_val;
				d_val = c;
				fd = fc;
				c = b - alpha_gold * (b - a);
				x_c = x_current + c * direction;
				fc = m2d(ff(x_c, ud1, ud2));
			} else {
				a = c;
				c = d_val;
				fc = fd;
				d_val = a + alpha_gold * (b - a);
				x_d = x_current + d_val * direction;
				fd = m2d(ff(x_d, ud1, ud2));
			}
			ls_iters++;
		}
		
		return (a + b) / 2.0;
	}
	catch (string ex_info)
	{
		throw ("double golden_line_search(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		Xopt.x = x0;
		Xopt.fit_fun(ff, ud1, ud2);
		
		solution XB = Xopt;
		
		int iteration = 0;
		const int MAX_ITERATIONS = 10000; // Maximum iterations to prevent infinite loops
		
		while (true) {
			// Compute gradient using solution class method
			matrix grad = XB.grad(gf, ud1, ud2);
			
			// Check gradient norm for convergence
			if (norm(grad) < epsilon) {
				Xopt.flag = 1;
				break;
			}
			
			if (solution::f_calls > Nmax || iteration >= MAX_ITERATIONS) {
				Xopt.flag = 0;
				break;
			}
			
			// Direction is negative gradient
			matrix d = -grad;
			
			// If h0 > 0, use fixed step size
			// If h0 == 0, use line search (golden section)
			double step_size;
			if (h0 > 0) {
				step_size = h0;
			} else {
				// Use golden section line search
				step_size = golden_line_search(ff, XB.x, d, 0.0, 2.0, 1e-6, Nmax, ud1, ud2);
			}
			
			// Update position
			Xopt.x = XB.x + step_size * d;
			Xopt.fit_fun(ff, ud1, ud2);
			
			// Update best solution
			XB = Xopt;
			
			iteration++;
		}
		
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
		Xopt.x = x0;
		Xopt.fit_fun(ff, ud1, ud2);
		
		solution XB = Xopt;
		matrix grad = XB.grad(gf, ud1, ud2);
		matrix d = -grad; // Initial direction
		matrix grad_prev;
		
		int iteration = 0;
		const int MAX_ITERATIONS = 10000; // Maximum iterations to prevent infinite loops
		
		while (true) {
			// Check gradient norm for convergence
			if (norm(grad) < epsilon) {
				Xopt.flag = 1;
				break;
			}
			
			if (solution::f_calls > Nmax || iteration >= MAX_ITERATIONS) {
				Xopt.flag = 0;
				break;
			}
			
			// If h0 > 0, use fixed step size
			// If h0 == 0, use line search (golden section)
			double step_size;
			if (h0 > 0) {
				step_size = h0;
			} else {
				// Use golden section line search
				step_size = golden_line_search(ff, XB.x, d, 0.0, 2.0, 1e-6, Nmax, ud1, ud2);
			}
			
			// Update position
			Xopt.x = XB.x + step_size * d;
			Xopt.fit_fun(ff, ud1, ud2);
			
			// Compute new gradient using solution class method
			grad_prev = grad;
			grad = Xopt.grad(gf, ud1, ud2);
			
			// Fletcher-Reeves formula: beta = ||grad_new||^2 / ||grad_old||^2
			double grad_norm_sq = norm(grad) * norm(grad);
			double grad_prev_norm_sq = norm(grad_prev) * norm(grad_prev);
			
			double beta = 0.0;
			if (grad_prev_norm_sq > 1e-15) {
				beta = grad_norm_sq / grad_prev_norm_sq;
			}
			
			// Update direction: d = -grad + beta * d_prev
			d = -grad + beta * d;
			
			// Update best solution
			XB = Xopt;
			
			iteration++;
		}
		
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
		Xopt.x = x0;
		Xopt.fit_fun(ff, ud1, ud2);
		
		solution XB = Xopt;
		
		int iteration = 0;
		const int MAX_ITERATIONS = 10000; // Maximum iterations to prevent infinite loops
		
		while (true) {
			// Compute gradient using solution class method
			matrix grad = XB.grad(gf, ud1, ud2);
			
			// Check gradient norm for convergence
			if (norm(grad) < epsilon) {
				Xopt.flag = 1;
				break;
			}
			
			if (solution::f_calls > Nmax || iteration >= MAX_ITERATIONS) {
				Xopt.flag = 0;
				break;
			}
			
			// Compute Hessian using solution class method
			matrix H = XB.hess(Hf, ud1, ud2);
			
			// Newton direction: d = -H^{-1} * grad
			matrix H_inv = inv(H);
			matrix d = -(H_inv * grad);
			
			// If h0 > 0, use fixed step size
			// If h0 == 0, use line search (golden section)
			double step_size;
			if (h0 > 0) {
				step_size = h0;
			} else {
				// Use golden section line search
				step_size = golden_line_search(ff, XB.x, d, 0.0, 2.0, 1e-6, Nmax, ud1, ud2);
			}
			
			// Update position
			Xopt.x = XB.x + step_size * d;
			Xopt.fit_fun(ff, ud1, ud2);
			
			// Update best solution
			XB = Xopt;
			
			iteration++;
		}
		
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
		
		double alpha = (sqrt(5.0) - 1.0) / 2.0; // golden ratio constant
		
		double c = b - alpha * (b - a);
		double d = a + alpha * (b - a);
		
		solution fc(c);
		solution fd(d);
		fc.fit_fun(ff, ud1, ud2);
		fd.fit_fun(ff, ud1, ud2);
		
		while (true) {
			if (solution::f_calls > Nmax) {
				Xopt.flag = 0;
				break;
			}
			
			if (b - a < epsilon) {
				Xopt.x = (a + b) / 2.0;
				Xopt.fit_fun(ff, ud1, ud2);
				Xopt.flag = 1;
				break;
			}
			
			if (fc.y < fd.y) {
				b = d;
				d = c;
				fd = fc;
				c = b - alpha * (b - a);
				fc.x = c;
				fc.fit_fun(ff, ud1, ud2);
			} else {
				a = c;
				c = d;
				fc = fd;
				d = a + alpha * (b - a);
				fd.x = d;
				fd.fit_fun(ff, ud1, ud2);
			}
		}
		
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
