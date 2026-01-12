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
		lab5();
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
	aN = std::clamp(aN, 0, 2); // wybor alpha
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
	// // TESTOWA FUNKCJA CELU
	// double epsilon = 1e-3;
    // int Nmax = 10000;
    
    // double steps[] = {0.05, 0.25, 0.0};
    
    // ofstream Sout("lab4_tabela1.csv");
 
    // Sout << "Lp.,x1(0),x2(0),x1*,x2*,y*,f_calls,g_calls,,x1*,x2*,y*,f_calls,g_calls,,x1*,x2*,y*,f_calls,g_calls,H_calls" << endl;
 
    // for (double s : steps) {
    //     for (int i = 0; i < 100; i++) {
    //         matrix x0 = rand_mat(2);
    //         x0(0) = x0(0) * 4.0 - 2.0;
    //         x0(1) = x0(1) * 4.0 - 2.0;
 
    //         Sout << (i + 1) << "," << x0(0) << "," << x0(1) << ",";
            
    //         solution::clear_calls();
    //         solution optSD = SD(ff4T, gf4T, x0, s, epsilon, Nmax);
    //         Sout << optSD.x(0) << "," << optSD.x(1) << "," << optSD.y(0) << ","
    //             << solution::f_calls << "," << solution::g_calls << ",,";
 
    //         solution::clear_calls();
    //         solution optCG = CG(ff4T, gf4T, x0, s, epsilon, Nmax);
    //         Sout << optCG.x(0) << "," << optCG.x(1) << "," << optCG.y(0) << ","
    //         	<< solution::f_calls << "," << solution::g_calls << ",,";
            
    //         solution::clear_calls();
    //         solution optNewton = Newton(ff4T, gf4T, Hf4T, x0, s, epsilon, Nmax);
    //         Sout << optNewton.x(0) << "," << optNewton.x(1) << "," << optNewton.y(0) << ","
    //         	<< solution::f_calls << "," << solution::g_calls << "," << solution::H_calls << "\n";
    //     }
    // }
	// Sout.close();
	
	// // PROBLEM RZECZYWISTY
	// Sout.open("table2_lab4.csv");
	// Sout << "Method,StepSize,Theta0,Theta1,Theta2,J,Accuracy,FCalls" << endl;
	
	// ifstream xfile("XData.txt");
	// ifstream yfile("YData.txt");
	
	// if (!xfile.is_open() || !yfile.is_open()) {
	// 	cerr << "Error: Could not open data files" << endl;
	// 	return;
	// }
	
	// int m = 0;
	// string line;
	// while (getline(yfile, line)) m++;
	// yfile.close();
	// yfile.open("YData.txt");
	
	// matrix X(m, 3);
	// matrix Y(m, 1);
	
	// // Read data
	// for (int i = 0; i < m; i++) {
	// 	double x1, x2;
	// 	xfile >> x1 >> x2;
	// 	X(i, 0) = 1.0;
	// 	X(i, 1) = x1;
	// 	X(i, 2) = x2;
		
	// 	double y;
	// 	yfile >> y;
	// 	Y(i, 0) = y;
	// }
	
	// xfile.close();
	// yfile.close();
	
	// matrix theta0(3, 1);
	// theta0(0) = 0.0;
	// theta0(1) = 0.0;
	// theta0(2) = 0.0;
	
	// matrix J0 = ff4R_cost(theta0, X, Y);
	// matrix grad0 = gf4R_grad(theta0, X, Y);
	
	// // Fixed step optimizations
	// double lr_sd_steps[] = {0.05, 0.25};
	// double lr_cg_steps[] = {0.05, 0.25};
	// double lr_newton_steps[] = {0.01, 0.0001};
	
	// solution best_solution;
	// double best_accuracy = 0.0;
	
	// auto compute_accuracy = [&](matrix theta) -> double {
	// 	int correct = 0;
	// 	for (int i = 0; i < m; i++) {
	// 		double z = 0.0;
	// 		for (int j = 0; j < 3; j++) {
	// 			z += m2d(theta(j)) * m2d(X(i, j));
	// 		}
	// 		double h = 1.0 / (1.0 + exp(-z));
	// 		int predicted = (h >= 0.5) ? 1 : 0;
	// 		if (predicted == (int)m2d(Y(i))) correct++;
	// 	}
	// 	return (double)correct / m;
	// };
	
	// for (int s_idx = 0; s_idx < 2; s_idx++) {
	// 	double step = lr_sd_steps[s_idx];
		
	// 	solution::clear_calls();
	// 	solution opt = SD(ff4R_cost, gf4R_grad, theta0, step, epsilon, Nmax, X, Y);
	// 	double accuracy = compute_accuracy(opt.x);
		
	// 	Sout << "SD_fixed," << step << "," << opt.x(0) << "," << opt.x(1) << "," << opt.x(2) 
	// 	     << "," << opt.y(0) << "," << accuracy << "," << solution::f_calls << endl;
		
	// 	if (accuracy > best_accuracy) {
	// 		best_accuracy = accuracy;
	// 		best_solution = opt;
	// 	}
	// }
	
	// solution::clear_calls();
	// solution opt_sd_ls = SD(ff4R_cost, gf4R_grad, theta0, 0.0, epsilon, Nmax, X, Y);
	// double acc_sd_ls = compute_accuracy(opt_sd_ls.x);
	
	// Sout << "SD_linesearch,0," << opt_sd_ls.x(0) << "," << opt_sd_ls.x(1) << "," << opt_sd_ls.x(2) 
	//      << "," << opt_sd_ls.y(0) << "," << acc_sd_ls << "," << solution::f_calls << endl;
	
	// if (acc_sd_ls > best_accuracy) {
	// 	best_accuracy = acc_sd_ls;
	// 	best_solution = opt_sd_ls;
	// }
	
	// for (int s_idx = 0; s_idx < 2; s_idx++) {
	// 	double step = lr_cg_steps[s_idx];
		
	// 	solution::clear_calls();
	// 	solution opt = CG(ff4R_cost, gf4R_grad, theta0, step, epsilon, Nmax, X, Y);
	// 	double accuracy = compute_accuracy(opt.x);
		
	// 	Sout << "CG_fixed," << step << "," << opt.x(0) << "," << opt.x(1) << "," << opt.x(2) 
	// 	     << "," << opt.y(0) << "," << accuracy << "," << solution::f_calls << endl;
		
	// 	if (accuracy > best_accuracy) {
	// 		best_accuracy = accuracy;
	// 		best_solution = opt;
	// 	}
	// }
	
	// Sout << "CG_linesearch,0,N/A,N/A,N/A,N/A,N/A,N/A" << endl;
	// Sout.close();
	
	// ofstream boundary_file("decision_boundary_lab4.csv");
	// boundary_file << "x1,x2,y,predicted" << endl;
	
	// for (int i = 0; i < m; i++) {
	// 	double x1 = m2d(X(i, 1));
	// 	double x2 = m2d(X(i, 2));
	// 	double y = m2d(Y(i));
		
	// 	double z = m2d(best_solution.x(0)) + m2d(best_solution.x(1)) * x1 + m2d(best_solution.x(2)) * x2;
	// 	double h = 1.0 / (1.0 + exp(-z));
	// 	int predicted = (h >= 0.5) ? 1 : 0;
		
	// 	boundary_file << x1 << "," << x2 << "," << y << "," << predicted << endl;
	// }
	
	// boundary_file.close();

	// // Select one starting point for trajectory visualization
	// matrix x0_plot(2, 1);
	// x0_plot(0) = -1.66484;  // x1
	// x0_plot(1) = -1.28852;   // x2
	
	// double epsilon_plot = 1e-3;
	// int Nmax_plot = 10000;
	
	// // SD with h=0.05
	// g_track_trajectory = true;
	// g_trajectory_file.open("trajectory_SD_h005.csv");
	// g_trajectory_file << "x1,x2,f" << endl;
	// solution::clear_calls();
	// SD(ff4T, gf4T, x0_plot, 0.05, epsilon_plot, Nmax_plot);
	// g_trajectory_file.close();
	// g_track_trajectory = false;
	
	// // SD with h=0.25
	// g_track_trajectory = true;
	// g_trajectory_file.open("trajectory_SD_h025.csv");
	// g_trajectory_file << "x1,x2,f" << endl;
	// solution::clear_calls();
	// SD(ff4T, gf4T, x0_plot, 0.25, epsilon_plot, Nmax_plot);
	// g_trajectory_file.close();
	// g_track_trajectory = false;
	
	// // SD with line search
	// g_track_trajectory = true;
	// g_trajectory_file.open("trajectory_SD_linesearch.csv");
	// g_trajectory_file << "x1,x2,f" << endl;
	// solution::clear_calls();
	// SD(ff4T, gf4T, x0_plot, 0.0, epsilon_plot, Nmax_plot);
	// g_trajectory_file.close();
	// g_track_trajectory = false;
	
	// // CG with h=0.05
	// g_track_trajectory = true;
	// g_trajectory_file.open("trajectory_CG_h005.csv");
	// g_trajectory_file << "x1,x2,f" << endl;
	// solution::clear_calls();
	// CG(ff4T, gf4T, x0_plot, 0.05, epsilon_plot, Nmax_plot);
	// g_trajectory_file.close();
	// g_track_trajectory = false;
	
	// // CG with h=0.25
	// g_track_trajectory = true;
	// g_trajectory_file.open("trajectory_CG_h025.csv");
	// g_trajectory_file << "x1,x2,f" << endl;
	// solution::clear_calls();
	// CG(ff4T, gf4T, x0_plot, 0.25, epsilon_plot, Nmax_plot);
	// g_trajectory_file.close();
	// g_track_trajectory = false;
	
	// // CG with line search
	// g_track_trajectory = true;
	// g_trajectory_file.open("trajectory_CG_linesearch.csv");
	// g_trajectory_file << "x1,x2,f" << endl;
	// solution::clear_calls();
	// CG(ff4T, gf4T, x0_plot, 0.0, epsilon_plot, Nmax_plot);
	// g_trajectory_file.close();
	// g_track_trajectory = false;
	
	// // Newton with h=0.05
	// g_track_trajectory = true;
	// g_trajectory_file.open("trajectory_Newton_h005.csv");
	// g_trajectory_file << "x1,x2,f" << endl;
	// solution::clear_calls();
	// Newton(ff4T, gf4T, Hf4T, x0_plot, 0.05, epsilon_plot, Nmax_plot);
	// g_trajectory_file.close();
	// g_track_trajectory = false;
	
	// // Newton with h=0.25
	// g_track_trajectory = true;
	// g_trajectory_file.open("trajectory_Newton_h025.csv");
	// g_trajectory_file << "x1,x2,f" << endl;
	// solution::clear_calls();
	// Newton(ff4T, gf4T, Hf4T, x0_plot, 0.25, epsilon_plot, Nmax_plot);
	// g_trajectory_file.close();
	// g_track_trajectory = false;
	
	// // Newton with line search
	// g_track_trajectory = true;
	// g_trajectory_file.open("trajectory_Newton_linesearch.csv");
	// g_trajectory_file << "x1,x2,f" << endl;
	// solution::clear_calls();
	// Newton(ff4T, gf4T, Hf4T, x0_plot, 0.0, epsilon_plot, Nmax_plot);
	// g_trajectory_file.close();
	// g_track_trajectory = false;
}

void lab5()
{
    // TESTOWA FUNKCJA CELU
    double epsilon = 1e-3;
    int Nmax = 10000;
   
    double a_values[] = {1.0, 10.0, 100.0};
   
    ofstream Sout("lab5_tabela1.csv");
 
    Sout << "w,x1(0),x2(0),x1*,x2*,f1*,f2*,f_calls,x1*,x2*,f1*,f2*,f_calls,x1*,x2*,f1*,f2*,f_calls" << endl;
 
    for (int i = 0; i <= 100; i++) {
        double w = i / 100.0;
        
        matrix x0 = rand_mat(2);
        x0(0) = x0(0) * 20.0 - 10.0;
        x0(1) = x0(1) * 20.0 - 10.0;
 
        Sout << w << "," << x0(0) << "," << x0(1);
        
        for (int j = 0; j < 3; j++) {
            double a = a_values[j];
            
            matrix params(2, 1);
            params(0) = a;
            params(1) = w;
            
            solution::clear_calls();
            solution opt = Powell(ff5T, x0, epsilon, Nmax, matrix(NAN), params);
            
            double f1 = a * (pow(opt.x(0) - 3.0, 2) + pow(opt.x(1) - 3.0, 2));
            double f2 = (1.0 / a) * (pow(opt.x(0) + 3.0, 2) + pow(opt.x(1) + 3.0, 2));

            Sout << "," << opt.x(0) << "," << opt.x(1) << ","
                 << f1 << "," << f2 << ","
                 << solution::f_calls;
        }
        
        Sout << endl;
    }
    
    Sout.close();
}

void lab6()
{
}
