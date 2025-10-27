/*********************************************
Kod stanowi uzupe�nienie materia��w do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udost�pniony na licencji CC BY-SA 3.0
Autor: dr in�. �ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G�rniczo-Hutnicza
Data ostatniej modyfikacji: 30.09.2025
*********************************************/

#include"opt_alg.h"

void lab0();
void lab1(int aN);
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main(int argc, char* argv[])
{
	try
	{
		lab1(atoi(argv[1]));
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	return 0;
}

void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;									// dok�adno��
	int Nmax = 10000;										// maksymalna liczba wywo�a� funkcji celu
	matrix lb(2, 1, -5), ub(2, 1, 5),						// dolne oraz g�rne ograniczenie
		a(2, 1);											// dok�adne rozwi�zanie optymalne
	solution opt;											// rozwi�zanie optymalne znalezione przez algorytm
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);			// wywo�anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie licznik�w

	//Wahadlo
	Nmax = 1000;											// dok�adno��
	epsilon = 1e-2;											// maksymalna liczba wywo�a� funkcji celu
	lb = 0, ub = 5;											// dolne oraz g�rne ograniczenie
	double teta_opt = 1;									// maksymalne wychylenie wahad�a
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);		// wywo�anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie licznik�w

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz�tkowe
		MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });	// MT zawiera moment si�y dzia�aj�cy na wahad�o oraz czas dzia�ania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);	// rozwi�zujemy r�wnanie r�niczkowe
	ofstream Sout("symulacja_lab0.csv");					// definiujemy strumie� do pliku .csv
	Sout << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
	Sout.close();											// zamykamy strumie�
	Y[0].~matrix();											// usuwamy z pami�ci rozwi�zanie RR
	Y[1].~matrix();
}

void lab1(int aN)
{
#pragma region zadanie
	//DANE
	double alpha[] = { 1.1, 1.5, 2.0 };

	aN = std::clamp(aN, 0, 2);

	solution opt;
	int Nmax = 1000;
	double d = 1.5;

	double epsilon = 1e-3;
	double gamma = 1e-6;

	matrix cel = matrix(50.0);

	int x0 = 0;

	int ld = 0, ud = 100;

	double* p;

	cout.precision(10);

	int n = 1;

	for (int i = 0; i < n; i++)
	{
		//LOSOWANIE PUNKTU POCZATKOWEGO
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<> dis(ld, ud);
		x0 = dis(gen);

		//EKSPANSJA
		solution::clear_calls();
		p = expansion(ff1R, x0, d, alpha[aN], Nmax, cel);
		cout << "Ekspansja: ";
		cout << fixed << x0 << "," << p[0] << "," << p[1] << "," << solution::f_calls << ",\n";

		//FIBONACCI
		solution::clear_calls();
		opt = fib(ff1R, p[0], p[1], epsilon, cel);
		cout << "Fibonacci: ";
		cout << fixed << opt.x(0) << "," << opt.y(0) << "," << solution::f_calls << "," << "\n";

		//LAGRANGE 
		solution::clear_calls();
		try {
			opt = lag(ff1R, ld, ud, epsilon, gamma, Nmax, matrix(50.0), NAN);
		}
		catch (string ex_info) {
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

	matrix* Y = solve_ode(df1, 0.0, 1.0, 2000.0, Y0, matrix(xLag), NAN);

	cout << hcat(Y[0], Y[1]) << endl;
#pragma endregion
}

void lab2()
{

}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}
