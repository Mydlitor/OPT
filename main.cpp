/*********************************************
Kod stanowi uzupe�nienie materia��w do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udost�pniony na licencji CC BY-SA 3.0
Autor: dr in�. �ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G�rniczo-Hutnicza
Data ostatniej modyfikacji: 30.09.2025
*********************************************/

#include "opt_alg.h"

void lab0();
void lab1(int aN);
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main(int argc, char *argv[])
{
	try
	{
		// Current lab: Lab4 - Gradient-based optimization
		// Change this line to switch between different labs
		// To run lab5: uncomment the line below and comment out lab4()
		// lab5();
		lab6();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl
			 << endl;
	}
	return 0;
}

void lab0()
{
	// Funkcja testowa
	double epsilon = 1e-2;			  // dok�adno��
	int Nmax = 10000;				  // maksymalna liczba wywo�a� funkcji celu
	matrix lb(2, 1, -5), ub(2, 1, 5), // dolne oraz g�rne ograniczenie
		a(2, 1);					  // dok�adne rozwi�zanie optymalne
	solution opt;					  // rozwi�zanie optymalne znalezione przez algorytm
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a); // wywo�anie procedury optymalizacji
	cout << opt << endl
		 << endl;			 // wypisanie wyniku
	solution::clear_calls(); // wyzerowanie licznik�w

	// Wahadlo
	Nmax = 1000;										// dok�adno��
	epsilon = 1e-2;										// maksymalna liczba wywo�a� funkcji celu
	lb = 0, ub = 5;										// dolne oraz g�rne ograniczenie
	double teta_opt = 1;								// maksymalne wychylenie wahad�a
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt); // wywo�anie procedury optymalizacji
	cout << opt << endl
		 << endl;			 // wypisanie wyniku
	solution::clear_calls(); // wyzerowanie licznik�w

	// Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1),							 // Y0 zawiera warunki pocz�tkowe
		MT = matrix(2, new double[2]{m2d(opt.x), 0.5});	 // MT zawiera moment si�y dzia�aj�cy na wahad�o oraz czas dzia�ania
	matrix *Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT); // rozwi�zujemy r�wnanie r�niczkowe
	ofstream Sout("symulacja_lab0.csv");				 // definiujemy strumie� do pliku .csv
	Sout << hcat(Y[0], Y[1]);							 // zapisyjemy wyniki w pliku
	Sout.close();										 // zamykamy strumie�
	Y[0].~matrix();										 // usuwamy z pami�ci rozwi�zanie RR
	Y[1].~matrix();
}

void lab1(int aN)
{
#pragma region zadanie
	// DANE
	if (aN < 0) aN = 0; // wybor alpha
	if (aN > 2) aN = 2;
	solution opt;
	int Nmax = 1000;
	matrix cel = matrix(50.0);
	int x0 = 0;

	// ekspansja
	int ld = 0, ud = 100;
	double alpha[] = {1.1, 1.5, 2.0};
	double d = 1.5;
	double *p; // uzywane do obu jako przedzial

	// lagrange
	double epsilon = 1e-3; // dla fib i lag
	double gamma = 1e-6;

	cout.precision(10);

	int n = 1;

	for (int i = 0; i < n; i++)
	{
		// LOSOWANIE PUNKTU POCZATKOWEGO
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<> dis(ld, ud);
		x0 = dis(gen);

		// EKSPANSJA
		solution::clear_calls();
		p = expansion(ff1R, x0, d, alpha[aN], Nmax, cel);
		cout << "Ekspansja: ";
		cout << fixed << x0 << "," << p[0] << "," << p[1] << "," << solution::f_calls << ",\n";

		// FIBONACCI
		solution::clear_calls();
		opt = fib(ff1R, p[0], p[1], epsilon, cel);
		cout << "Fibonacci: ";
		cout << fixed << opt.x(0) << "," << opt.y(0) << "," << solution::f_calls << "," << "\n";

		// LAGRANGE
		solution::clear_calls();
		try
		{
			opt = lag(ff1R, ld, ud, epsilon, gamma, Nmax, cel, NAN);
		}
		catch (string ex_info)
		{
			cout << "error,error,error,\n";
			continue;
		}
		cout << "Lagrange: ";
		cout << fixed << opt.x(0) << "," << opt.y(0) << "," << solution::f_calls << "," << "\n";
	}
	delete p;
#pragma endregion

#pragma region symulacja
	double xFab = 20.1374350875;
	double xLag = 20.1380339112;

	matrix Y0 = matrix(3, 1);
	Y0(0) = 5.0;
	Y0(1) = 1.0;
	Y0(2) = 20.0;

	matrix *Y = solve_ode(df1, 0.0, 1.0, 2000.0, Y0, matrix(xLag), NAN);

	cout << hcat(Y[0], Y[1]) << endl;
#pragma endregion
}

void lab2()
{
	/*// matrix xy_start(2, 1, 100.0);

	// double step_size = 1.0;
	// double alfa = 0.8;
	// double epsilon = 1e-3;

	// int max_n = 1000;

	// solution opt;

	// opt = HJ(ff2T, xy_start, step_size, alfa, epsilon, max_n);

	// std::cout << "X: " << opt.x(0) << " Y: " << opt.x(1) << ", f(x,y): " << opt.y << "\n";

	// TESTOWA FUNKCJA CELU

	double alphaHJ = 0.5;
	double alphaRos = 2.0;
	double beta = 0.5;
	double epsilon = 1e-4;
	int Nmax = 10000;

	double step_sizes[3] = {0.01, 0.36, 0.876};
	int n = 100;

	solution optHJ, optRos;

	ofstream Sout("tabela3_lab2.csv");

	// cout << "Lp.,x1,x2,x1*,x2*,y*,f_calls" << endl;
	Sout << "Lp.,x1,x2,x1*,x2*,y*,f_calls,,x1*,x2*,y*,f_calls" << endl;

	double s = 0.12;
		matrix k0 = rand_mat(2);

		for (int j = 0; j < 2; ++j)
		{
			k0(j) = 20.0 * k0(j);
		}

		k0(0) = 3.00612;
		k0(1) = 10.8401;

		matrix s0(2, 1);
		s0(0) = s;
		s0(1) = s;

		// cout << (i + 1) << ",";
		//Sout << (i + 1) << "," << x0(0) << "," << x0(1);

		// METODA HOOKE'A JEEVESA
		cout << "HJ method:\n";
		solution::clear_calls();
		optHJ = HJ(ff2R, k0, s, alphaHJ, epsilon, Nmax);

		// cout << x0(0) << "," << x0(1) << "," << opt.x(0) << "," << opt.x(1) << "," << opt.y(0) << "," << opt.f_calls << endl;
		//cout << optHJ.x(0) << "," << optHJ.x(1) << "," << optHJ.y(0) << "," << optHJ.f_calls;

		// METODA ROSENBROCKA
		cout << "Rosenbrock method:\n";
		solution::clear_calls();
		optRos = Rosen(ff2R, k0, s0, alphaRos, beta, epsilon, Nmax);
		//cout << "," << optRos.x(0) << "," << optRos.x(1) << "," << optRos.y(0) << "," << optRos.f_calls << endl;

	Sout.close();


	// PROBLEM RZECZYWISTY - TEST
	/*matrix k(2, 1);
	k(0) = 5.0;
	k(1) = 5.0;

	matrix Q = ff2R(k);

	cout << "Q(k1,k2) = " << Q(0) << " (powinna wynosic: okolo 775.229)\n";*/

	// HJ
	/*matrix k0(2, 1);
	k0(0) = 3.00612;
	k0(1) = 10.8401;

	ff2R(k0);*/

	// Rosen
	matrix k0(2, 1);
	k0(0) = 3.00586;
	k0(1) = 10.8388;

	ff2R(k0);
}

void lab3()
{
	// Parametry algorytmu Nelder-Mead
	double epsilon = 1e-2;
	int Nmax = 10000;

	int n = 100; 

	// Wartości parametru a
	double a_values[3] = {4.0, 4.4934, 5.0};

	solution optX;

	cout.precision(6);

	cout << "=== ZEWNETRZNA FUNKCJA KARY ===" << endl;

	for (int a_idx = 0; a_idx < 3; ++a_idx)
	{
		double a = a_values[a_idx];
		cout << "\na = " << a << endl;
		cout << "Lp,x1,x2,x1*,x2*,y*,r*,f_calls" << endl;

		for (int i = 0; i < n; ++i)
		{
			matrix x0(2, 1);
			do
			{
				x0 = rand_mat(2);
				x0(0) = x0(0) * (a - 1.0) + 1.0;
				x0(1) = x0(1) * (a - 1.0) + 1.0;
			} while (sqrt(x0(0) * x0(0) + x0(1) * x0(1)) > a); 

			matrix params(2, 1);
			params(0) = a;
			params(1) = 1.0; 

			solution::clear_calls();

			optX = pen(ff3T_zewn, x0, 0.5, 1.5, epsilon, Nmax, params, NAN);

			matrix y_bez_kary = ff3T(optX.x); 

			double r = sqrt(optX.x(0) * optX.x(0) + optX.x(1) * optX.x(1));

			cout << fixed << (i + 1) << "," << x0(0) << "," << x0(1) << ","
				 << optX.x(0) << "," << optX.x(1) << "," << r << ","
				 << y_bez_kary(0) << "," << solution::f_calls << endl;
		}
	}

	cout << "\n=== WEWNETRZNA FUNKCJA KARY ===" << endl;

	for (int a_idx = 0; a_idx < 3; ++a_idx)
	{
		double a = a_values[a_idx];
		cout << "\na = " << a << endl;
		cout << "Lp,x1,x2,x1*,x2*,y*,r*,f_calls" << endl;

		for (int i = 0; i < n; ++i)
		{
			matrix x0(2, 1);
			do
			{
				x0 = rand_mat(2);
				x0(0) = x0(0) * (a - 1.0 - 0.1) + 1.0 + 0.05;
				x0(1) = x0(1) * (a - 1.0 - 0.1) + 1.0 + 0.05;
			} while (sqrt(x0(0) * x0(0) + x0(1) * x0(1)) >= a - 0.01); 

			matrix params(2, 1);
			params(0) = a;
			params(1) = 10.0;

			solution::clear_calls();

			optX = pen(ff3T_wewn, x0, 5.0, 0.5, epsilon, Nmax, params, NAN);

			matrix y_bez_kary = ff3T(optX.x);

			double r = sqrt(optX.x(0) * optX.x(0) + optX.x(1) * optX.x(1));

			cout << fixed << (i + 1) << "," << x0(0) << "," << x0(1) << ","
				 << optX.x(0) << "," << optX.x(1) << "," << r << ","
				 << y_bez_kary(0) << "," << solution::f_calls << endl;
		}
	}

	cout << "\n=== PROBLEM RZECZYWISTY ===" << endl;
	cout << "v0x,omega,x_end,f_calls" << endl;

	matrix x0_R(2, 1);
	x0_R(0) = 0.0;
	x0_R(1) = 0.0;

	matrix params_R(2, 1);
	params_R(0) = 1.0;
	params_R(1) = 1.0;

	solution::clear_calls();

	optX = pen(ff3R, x0_R, 1.0, 2.0, epsilon, Nmax, params_R, NAN);

	cout << fixed << optX.x(0) << "," << optX.x(1) << "," << -optX.y(0) << "," << solution::f_calls << endl;

	matrix Y0_sim(4, 1);
	Y0_sim(0) = 0.0;   
	Y0_sim(1) = 100.0; 
	Y0_sim(2) = optX.x(0);  
	Y0_sim(3) = 0.0;    

	matrix params_sim(2, 1);
	params_sim(0) = optX.x(0);  
	params_sim(1) = optX.x(1);  

	matrix *Y = solve_ode(df3, 0.0, 0.01, 7.0, Y0_sim, params_sim, NAN);  
	ofstream Sout("symulacja_lab3.csv");
	Sout << hcat(Y[0], Y[1]);
	Sout.close();
	Y[0].~matrix();
	Y[1].~matrix();
}

void lab4()
{
	cout.precision(10);
	
	// ============================================================
	// PART 1: TEST FUNCTION EXPERIMENTS
	// ============================================================
	cout << "=== TEST FUNCTION EXPERIMENTS ===" << endl;
	
	// Set random seed for reproducibility
	srand(42);
	
	double epsilon = 1e-3;
	int Nmax = 10000;
	int n_experiments = 100;
	
	// Generate one set of random starting points to be used by all methods
	// This ensures fair comparison between methods
	vector<matrix> starting_points;
	for (int i = 0; i < n_experiments; i++) {
		matrix x0(2, 1);
		x0(0) = (rand() / (double)RAND_MAX) * 4.0 - 2.0;
		x0(1) = (rand() / (double)RAND_MAX) * 4.0 - 2.0;
		starting_points.push_back(x0);
	}
	
	// Define step sizes for each method
	double sd_steps[] = {0.05, 0.25};
	double cg_steps[] = {0.05, 0.25};
	double newton_steps[] = {0.01, 0.0001};
	
	// Arrays to store results
	ofstream table1("table1_lab4.csv");
	table1 << "Method,StepSize,AvgFinalF,AvgIterations,AvgFCalls,Convergences" << endl;
	
	// Steepest Descent experiments
	for (int s_idx = 0; s_idx < 2; s_idx++) {
		double step = sd_steps[s_idx];
		double sum_f = 0.0, sum_iters = 0.0, sum_fcalls = 0.0;
		int convergences = 0;
		
		for (int i = 0; i < n_experiments; i++) {
			// Use the same starting point for all methods
			matrix x0 = starting_points[i];
			
			solution::clear_calls();
			solution opt = SD(ff4T, gf4T, x0, step, epsilon, Nmax);
			
			sum_f += m2d(opt.y);
			sum_fcalls += solution::f_calls;
			if (opt.flag == 1) convergences++;
		}
		
		table1 << "SD," << step << "," << (sum_f / n_experiments) << ","
		       << "N/A," << (sum_fcalls / n_experiments) << "," << convergences << endl;
		
		cout << "SD (h=" << step << "): avg_f=" << (sum_f / n_experiments) 
		     << ", avg_fcalls=" << (sum_fcalls / n_experiments) 
		     << ", convergences=" << convergences << "/" << n_experiments << endl;
	}
	
	// Conjugate Gradients experiments
	for (int s_idx = 0; s_idx < 2; s_idx++) {
		double step = cg_steps[s_idx];
		double sum_f = 0.0, sum_iters = 0.0, sum_fcalls = 0.0;
		int convergences = 0;
		
		for (int i = 0; i < n_experiments; i++) {
			// Use the same starting point for all methods
			matrix x0 = starting_points[i];
			
			solution::clear_calls();
			solution opt = CG(ff4T, gf4T, x0, step, epsilon, Nmax);
			
			sum_f += m2d(opt.y);
			sum_fcalls += solution::f_calls;
			if (opt.flag == 1) convergences++;
		}
		
		table1 << "CG," << step << "," << (sum_f / n_experiments) << ","
		       << "N/A," << (sum_fcalls / n_experiments) << "," << convergences << endl;
		
		cout << "CG (h=" << step << "): avg_f=" << (sum_f / n_experiments) 
		     << ", avg_fcalls=" << (sum_fcalls / n_experiments) 
		     << ", convergences=" << convergences << "/" << n_experiments << endl;
	}
	
	// Newton's method experiments
	for (int s_idx = 0; s_idx < 2; s_idx++) {
		double step = newton_steps[s_idx];
		double sum_f = 0.0, sum_iters = 0.0, sum_fcalls = 0.0;
		int convergences = 0;
		
		for (int i = 0; i < n_experiments; i++) {
			// Use the same starting point for all methods
			matrix x0 = starting_points[i];
			
			solution::clear_calls();
			solution opt = Newton(ff4T, gf4T, Hf4T, x0, step, epsilon, Nmax);
			
			sum_f += m2d(opt.y);
			sum_fcalls += solution::f_calls;
			if (opt.flag == 1) convergences++;
		}
		
		table1 << "Newton," << step << "," << (sum_f / n_experiments) << ","
		       << "N/A," << (sum_fcalls / n_experiments) << "," << convergences << endl;
		
		cout << "Newton (h=" << step << "): avg_f=" << (sum_f / n_experiments) 
		     << ", avg_fcalls=" << (sum_fcalls / n_experiments) 
		     << ", convergences=" << convergences << "/" << n_experiments << endl;
	}
	
	table1.close();
	
	// ============================================================
	// PART 2: LOGISTIC REGRESSION
	// ============================================================
	cout << "\n=== LOGISTIC REGRESSION ===" << endl;
	
	// Load data
	ifstream xfile("XData.txt");
	ifstream yfile("YData.txt");
	
	if (!xfile.is_open() || !yfile.is_open()) {
		cerr << "Error: Could not open data files" << endl;
		return;
	}
	
	// Count rows
	int m = 0;
	string line;
	while (getline(yfile, line)) m++;
	yfile.close();
	yfile.open("YData.txt");
	
	// Create data matrices
	matrix X(m, 3); // 3 columns: bias (1), x1, x2
	matrix Y(m, 1);
	
	// Read data
	for (int i = 0; i < m; i++) {
		double x1, x2;
		xfile >> x1 >> x2;
		X(i, 0) = 1.0; // bias term
		X(i, 1) = x1;
		X(i, 2) = x2;
		
		double y;
		yfile >> y;
		Y(i, 0) = y;
	}
	
	xfile.close();
	yfile.close();
	
	// Validation at theta = [0, 0, 0]
	matrix theta0(3, 1);
	theta0(0) = 0.0;
	theta0(1) = 0.0;
	theta0(2) = 0.0;
	
	matrix J0 = ff4R_cost(theta0, X, Y);
	matrix grad0 = gf4R_grad(theta0, X, Y);
	
	cout << "Validation at theta=[0,0,0]:" << endl;
	cout << "  J(theta) = " << J0(0) << endl;
	cout << "  grad J = [" << grad0(0) << ", " << grad0(1) << ", " << grad0(2) << "]" << endl;
	
	// Logistic regression optimization experiments
	ofstream table2("table2_lab4.csv");
	table2 << "Method,StepSize,Theta0,Theta1,Theta2,J,Accuracy,FCalls" << endl;
	
	// Fixed step optimizations
	double lr_sd_steps[] = {0.05, 0.25};
	double lr_cg_steps[] = {0.05, 0.25};
	double lr_newton_steps[] = {0.01, 0.0001};
	
	solution best_solution;
	double best_accuracy = 0.0;
	
	// Helper function to compute accuracy
	auto compute_accuracy = [&](matrix theta) -> double {
		int correct = 0;
		for (int i = 0; i < m; i++) {
			double z = 0.0;
			for (int j = 0; j < 3; j++) {
				z += m2d(theta(j)) * m2d(X(i, j));
			}
			double h = 1.0 / (1.0 + exp(-z));
			int predicted = (h >= 0.5) ? 1 : 0;
			if (predicted == (int)m2d(Y(i))) correct++;
		}
		return (double)correct / m;
	};
	
	// SD fixed step
	for (int s_idx = 0; s_idx < 2; s_idx++) {
		double step = lr_sd_steps[s_idx];
		
		solution::clear_calls();
		solution opt = SD(ff4R_cost, gf4R_grad, theta0, step, epsilon, Nmax, X, Y);
		double accuracy = compute_accuracy(opt.x);
		
		table2 << "SD_fixed," << step << "," << opt.x(0) << "," << opt.x(1) << "," << opt.x(2) 
		       << "," << opt.y(0) << "," << accuracy << "," << solution::f_calls << endl;
		
		cout << "SD (h=" << step << "): theta=[" << opt.x(0) << "," << opt.x(1) << "," << opt.x(2) 
		     << "], J=" << opt.y(0) << ", accuracy=" << accuracy << endl;
		
		if (accuracy > best_accuracy) {
			best_accuracy = accuracy;
			best_solution = opt;
		}
	}
	
	// SD with line search
	solution::clear_calls();
	solution opt_sd_ls = SD(ff4R_cost, gf4R_grad, theta0, 0.0, epsilon, Nmax, X, Y);
	double acc_sd_ls = compute_accuracy(opt_sd_ls.x);
	
	table2 << "SD_linesearch,0," << opt_sd_ls.x(0) << "," << opt_sd_ls.x(1) << "," << opt_sd_ls.x(2) 
	       << "," << opt_sd_ls.y(0) << "," << acc_sd_ls << "," << solution::f_calls << endl;
	
	cout << "SD (line search): theta=[" << opt_sd_ls.x(0) << "," << opt_sd_ls.x(1) << "," << opt_sd_ls.x(2) 
	     << "], J=" << opt_sd_ls.y(0) << ", accuracy=" << acc_sd_ls << endl;
	
	if (acc_sd_ls > best_accuracy) {
		best_accuracy = acc_sd_ls;
		best_solution = opt_sd_ls;
	}
	
	// CG fixed step
	for (int s_idx = 0; s_idx < 2; s_idx++) {
		double step = lr_cg_steps[s_idx];
		
		solution::clear_calls();
		solution opt = CG(ff4R_cost, gf4R_grad, theta0, step, epsilon, Nmax, X, Y);
		double accuracy = compute_accuracy(opt.x);
		
		table2 << "CG_fixed," << step << "," << opt.x(0) << "," << opt.x(1) << "," << opt.x(2) 
		       << "," << opt.y(0) << "," << accuracy << "," << solution::f_calls << endl;
		
		cout << "CG (h=" << step << "): theta=[" << opt.x(0) << "," << opt.x(1) << "," << opt.x(2) 
		     << "], J=" << opt.y(0) << ", accuracy=" << accuracy << endl;
		
		if (accuracy > best_accuracy) {
			best_accuracy = accuracy;
			best_solution = opt;
		}
	}
	
	// CG with line search - SKIPPED due to numerical issues
	cout << "CG (line search): Skipped (numerical instability)" << endl;
	table2 << "CG_linesearch,0,N/A,N/A,N/A,N/A,N/A,N/A" << endl;
	
	table2.close();
	
	// ============================================================
	// TABLE: CONJUGATE GRADIENT METHOD RESULTS
	// ============================================================
	cout << "\n=== TABELA: METODA GRADIENTÓW SPRZĘŻONYCH ===" << endl;
	cout << "Długość kroku | θ0*        | θ1*        | θ2*        | J(θ*)      | P(θ*)      | g_calls" << endl;
	cout << "----------------------------------------------------------------------------------------------" << endl;
	
	double cg_table_steps[] = {0.01, 0.001, 0.0001};
	
	for (int i = 0; i < 3; i++) {
		double step = cg_table_steps[i];
		
		solution::clear_calls();
		solution opt = CG(ff4R_cost, gf4R_grad, theta0, step, epsilon, Nmax, X, Y);
		double accuracy = compute_accuracy(opt.x);
		
		cout << left << setw(14) << step 
		     << "| " << setw(11) << opt.x(0) 
		     << "| " << setw(11) << opt.x(1) 
		     << "| " << setw(11) << opt.x(2) 
		     << "| " << setw(11) << opt.y(0) 
		     << "| " << setw(11) << accuracy 
		     << "| " << solution::f_calls << endl;
	}
	
	cout << endl;
	
	// Save decision boundary data
	cout << "\nBest model: theta=[" << best_solution.x(0) << "," << best_solution.x(1) 
	     << "," << best_solution.x(2) << "], accuracy=" << best_accuracy << endl;
	
	// Generate decision boundary plot data
	ofstream boundary_file("decision_boundary_lab4.csv");
	boundary_file << "x1,x2,y,predicted" << endl;
	
	for (int i = 0; i < m; i++) {
		double x1 = m2d(X(i, 1));
		double x2 = m2d(X(i, 2));
		double y = m2d(Y(i));
		
		double z = m2d(best_solution.x(0)) + m2d(best_solution.x(1)) * x1 + m2d(best_solution.x(2)) * x2;
		double h = 1.0 / (1.0 + exp(-z));
		int predicted = (h >= 0.5) ? 1 : 0;
		
		boundary_file << x1 << "," << x2 << "," << y << "," << predicted << endl;
	}
	
	boundary_file.close();
	
	cout << "\nResults saved to:" << endl;
	cout << "  - table1_lab4.csv (test function experiments)" << endl;
	cout << "  - table2_lab4.csv (logistic regression results)" << endl;
	cout << "  - decision_boundary_lab4.csv (decision boundary data)" << endl;
	
	// ============================================================
	// PART 3: GENERATE ITERATION HISTORY FOR CONTOUR PLOTS
	// ============================================================
	cout << "\n=== GENERATING ITERATION HISTORY FOR CONTOUR PLOTS ===" << endl;
	
	// Use the first starting point for consistency across all methods
	matrix x0_plot = starting_points[0];
	cout << "Using starting point: x0 = [" << x0_plot(0) << ", " << x0_plot(1) << "]" << endl;
	
	// Helper function to save history to CSV
	auto save_history = [](const std::vector<matrix>& history, const string& filename) {
		ofstream file(filename);
		file << "iteration,x1,x2" << endl;
		for (size_t i = 0; i < history.size(); i++) {
			file << i << "," << history[i](0) << "," << history[i](1) << endl;
		}
		file.close();
		cout << "  Saved: " << filename << " (" << history.size() << " iterations)" << endl;
	};
	
	// Steepest Descent with different step sizes
	{
		// SD with h = 0.05
		std::vector<matrix> history;
		solution::clear_calls();
		SD(ff4T, gf4T, x0_plot, 0.05, epsilon, Nmax, NAN, NAN, &history);
		save_history(history, "history_SD_0.05.csv");
		
		// SD with h = 0.25
		history.clear();
		solution::clear_calls();
		SD(ff4T, gf4T, x0_plot, 0.25, epsilon, Nmax, NAN, NAN, &history);
		save_history(history, "history_SD_0.25.csv");
		
		// SD with variable step (line search)
		history.clear();
		solution::clear_calls();
		SD(ff4T, gf4T, x0_plot, 0.0, epsilon, Nmax, NAN, NAN, &history);
		save_history(history, "history_SD_variable.csv");
	}
	
	// Conjugate Gradients with different step sizes
	{
		// CG with h = 0.05
		std::vector<matrix> history;
		solution::clear_calls();
		CG(ff4T, gf4T, x0_plot, 0.05, epsilon, Nmax, NAN, NAN, &history);
		save_history(history, "history_CG_0.05.csv");
		
		// CG with h = 0.25
		history.clear();
		solution::clear_calls();
		CG(ff4T, gf4T, x0_plot, 0.25, epsilon, Nmax, NAN, NAN, &history);
		save_history(history, "history_CG_0.25.csv");
		
		// CG with variable step (line search)
		history.clear();
		solution::clear_calls();
		CG(ff4T, gf4T, x0_plot, 0.0, epsilon, Nmax, NAN, NAN, &history);
		save_history(history, "history_CG_variable.csv");
	}
	
	// Newton's method with different step sizes
	{
		// Newton with h = 0.05 (for plot 1)
		std::vector<matrix> history;
		solution::clear_calls();
		Newton(ff4T, gf4T, Hf4T, x0_plot, 0.05, epsilon, Nmax, NAN, NAN, &history);
		save_history(history, "history_Newton_0.05.csv");
		
		// Newton with h = 0.25 (for plot 2)
		history.clear();
		solution::clear_calls();
		Newton(ff4T, gf4T, Hf4T, x0_plot, 0.25, epsilon, Nmax, NAN, NAN, &history);
		save_history(history, "history_Newton_0.25.csv");
		
		// Newton with h = 0.01 (for plot 6)
		history.clear();
		solution::clear_calls();
		Newton(ff4T, gf4T, Hf4T, x0_plot, 0.01, epsilon, Nmax, NAN, NAN, &history);
		save_history(history, "history_Newton_0.01.csv");
		
		// Newton with h = 0.0001 (for plot 6)
		history.clear();
		solution::clear_calls();
		Newton(ff4T, gf4T, Hf4T, x0_plot, 0.0001, epsilon, Nmax, NAN, NAN, &history);
		save_history(history, "history_Newton_0.0001.csv");
		
		// Newton with variable step (line search)
		history.clear();
		solution::clear_calls();
		Newton(ff4T, gf4T, Hf4T, x0_plot, 0.0, epsilon, Nmax, NAN, NAN, &history);
		save_history(history, "history_Newton_variable.csv");
	}
	
	cout << "\nAll history files generated successfully!" << endl;
	cout << "Run 'python3 plot_contours.py' to generate contour plots." << endl;
}

void lab5()
{
	cout.precision(10);
	
	// ============================================================
	// PART 1: TEST PROBLEM - MULTI-CRITERIA OPTIMIZATION
	// ============================================================
	cout << "=== LAB5: MULTI-CRITERIA OPTIMIZATION ===" << endl;
	cout << "=== PART 1: TEST PROBLEM ===" << endl;
	
	// Set random seed for reproducibility (can be changed to time(NULL) for true randomness)
	srand(42);
	
	double epsilon = 1e-3;
	int Nmax = 10000;
	
	// Parameter values for a
	double a_values[] = {1.0, 10.0, 100.0};
	
	// For each value of parameter a
	for (int a_idx = 0; a_idx < 3; a_idx++) {
		double a = a_values[a_idx];
		
		cout << "\n--- Testing with a = " << a << " ---" << endl;
		
		// Create output file
		stringstream filename;
		filename << "lab5_test_a" << (int)a << ".csv";
		ofstream outfile(filename.str());
		outfile << "w,x1_start,x2_start,x1_opt,x2_opt,f1_opt,f2_opt,f_opt,f_calls" << endl;
		outfile.precision(10);
		
		// 101 optimizations for w = 0.00, 0.01, 0.02, ..., 1.00
		for (int w_idx = 0; w_idx <= 100; w_idx++) {
			double w = w_idx / 100.0;
			
			// Random starting point in [-10, 10] x [-10, 10]
			matrix x0(2, 1);
			x0(0) = (rand() / (double)RAND_MAX) * 20.0 - 10.0;
			x0(1) = (rand() / (double)RAND_MAX) * 20.0 - 10.0;
			
			// Setup parameters: ud1(0) = a, ud1(1) = w
			matrix params(2, 1);
			params(0) = a;
			params(1) = w;
			
			// Run Powell optimization
			solution::clear_calls();
			solution opt = Powell(ff5T, x0, epsilon, Nmax, params, NAN);
			
			// Calculate individual criteria at optimal point
			matrix f1_opt = ff5_f1(opt.x, params, NAN);
			matrix f2_opt = ff5_f2(opt.x, params, NAN);
			
			// Write results
			outfile << w << ","
			        << x0(0) << "," << x0(1) << ","
			        << opt.x(0) << "," << opt.x(1) << ","
			        << f1_opt(0) << "," << f2_opt(0) << ","
			        << opt.y(0) << "," << solution::f_calls << endl;
			
			// Print progress every 10 iterations
			if (w_idx % 10 == 0) {
				cout << "  w=" << w << ": x*=[" << opt.x(0) << "," << opt.x(1) 
				     << "], f1*=" << f1_opt(0) << ", f2*=" << f2_opt(0) << endl;
			}
		}
		
		outfile.close();
		cout << "Results saved to " << filename.str() << endl;
	}
	
	// ============================================================
	// PART 2: REAL PROBLEM - CANTILEVER BEAM
	// ============================================================
	cout << "\n=== PART 2: REAL PROBLEM - CANTILEVER BEAM ===" << endl;
	
	// Validation test - verify formulas work
	cout << "\nValidation test (l=500mm, d=30mm):" << endl;
	{
		matrix x_test(2, 1);
		x_test(0) = 30.0;   // d = 30 mm
		x_test(1) = 500.0;  // l = 500 mm
		
		// Constants for validation
		const double P = 3000.0;
		const double E = 120e9;
		const double rho = 8920.0;
		
		double d_m = 30.0 / 1000.0;
		double l_m = 500.0 / 1000.0;
		
		double mass = rho * M_PI * pow(d_m / 2.0, 2) * l_m;
		double u = (64.0 * P * pow(l_m, 3)) / (3.0 * E * M_PI * pow(d_m, 4)) * 1000.0;
		double sigma = (32.0 * P * l_m) / (M_PI * pow(d_m, 3));
		
		cout << "  m = " << mass << " kg" << endl;
		cout << "  u = " << u << " mm" << endl;
		cout << "  sigma = " << sigma / 1e6 << " MPa" << endl;
		cout << "  Constraints: u_max=2.5mm, sigma_max=300MPa" << endl;
		if (u <= 2.5) cout << "  Deflection OK" << endl;
		else cout << "  Deflection VIOLATED" << endl;
		if (sigma <= 300e6) cout << "  Stress OK" << endl;
		else cout << "  Stress VIOLATED" << endl;
	}
	
	// Create output file for beam optimization
	ofstream beam_file("lab5_beam.csv");
	beam_file << "w,d_start,l_start,d_opt,l_opt,mass,deflection,stress,f_opt,f_calls" << endl;
	beam_file.precision(10);
	
	cout << "\nRunning 101 optimizations for beam problem..." << endl;
	
	// 101 optimizations for w = 0.00, 0.01, 0.02, ..., 1.00
	for (int w_idx = 0; w_idx <= 100; w_idx++) {
		double w = w_idx / 100.0;
		
		// Random starting point: d in [0.01, 1000], l in [0.2, 1000]
		matrix x0(2, 1);
		x0(0) = (rand() / (double)RAND_MAX) * 999.99 + 0.01;   // d
		x0(1) = (rand() / (double)RAND_MAX) * 999.8 + 0.2;     // l
		
		// Setup parameters: ud1(0) = w
		matrix params(1, 1);
		params(0) = w;
		
		// Run Powell optimization
		solution::clear_calls();
		solution opt = Powell(ff5R, x0, epsilon, Nmax, params, NAN);
		
		// Calculate mass, deflection, and stress at optimal point
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
		
		// Write results
		beam_file << w << ","
		          << x0(0) << "," << x0(1) << ","
		          << d << "," << l << ","
		          << mass << "," << u << "," << sigma / 1e6 << ","
		          << opt.y(0) << "," << solution::f_calls << endl;
		
		// Print progress every 10 iterations
		if (w_idx % 10 == 0) {
			cout << "  w=" << w << ": d*=" << d << "mm, l*=" << l 
			     << "mm, m=" << mass << "kg, u=" << u << "mm, sigma=" 
			     << sigma / 1e6 << "MPa" << endl;
		}
	}
	
	beam_file.close();
	cout << "\nResults saved to lab5_beam.csv" << endl;
	
	cout << "\n=== LAB5 COMPLETED ===" << endl;
	cout << "Generated files:" << endl;
	cout << "  - lab5_test_a1.csv (test problem with a=1)" << endl;
	cout << "  - lab5_test_a10.csv (test problem with a=10)" << endl;
	cout << "  - lab5_test_a100.csv (test problem with a=100)" << endl;
	cout << "  - lab5_beam.csv (cantilever beam problem)" << endl;
}

void lab6()
{
	int N = 2;
	int mi = 20;
	int lambda = 40;
	double epsilon = 1e-5;
	int Nmax = 10000;

	matrix lb(N, 1);
	lb(0) = -5;
	lb(1) = -5;

	matrix ub(N, 1);
	ub(0) = 5;
	ub(1) = 5;

	double sigma_values[] = { 0.01, 0.1, 1, 10, 100 };

	ofstream Sout("lab6_tabela1.csv");

	Sout << "Lp.,x1(0),x2(0),x1*,x2*,y*,f_calls,,x1*,x2*,y*,f_calls,,x1*,x2*,y*,f_calls,,x1*,x2*,y*,f_calls,,x1*,x2*,y*,f_calls" << endl;

	for (int i = 0; i < 100; i++) {
		matrix x0 = rand_mat(2);
		x0(0) = x0(0) * 10.0 - 5.0;
		x0(1) = x0(1) * 10.0 - 5.0;

		Sout << i << "," << x0(0) << "," << x0(1) << ",";

		for (double sigma : sigma_values) {
			matrix sigma0(1, 1);
			sigma0(0) = sigma;
			
			solution opt = EA(ff6T, N, lb, ub, mi, lambda, sigma0, epsilon, Nmax, NAN, NAN);

			Sout << opt.x(0) << "," << opt.x(1) << "," << opt.y(0) << "," << solution::f_calls << ",,";

			solution::clear_calls();
		}

		Sout << endl;
	}
}