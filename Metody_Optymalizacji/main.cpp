/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

#include <time.h>
#include"opt_alg.h"
#include <fstream>
#include <iomanip>
#include <random>

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main() {

	// algorytmy s¹ w opt_alg.cpp
	
	try {
		lab1();
	} catch (string EX_INFO) {
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	
	system("pause");
	return 0;
}

void lab0() {
	//Funkcja testowa
	double epsilon = 1e-2;
	int Nmax = 10000;
	matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);
	solution opt;
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Wahadlo
	Nmax = 1000;
	epsilon = 1e-2;
	lb = 0;
	ub = 5;
	double teta_opt = 1;
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(opt.x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
	ofstream Sout("symulacja_lab0.csv");
	Sout << hcat(Y[0], Y[1]);
	Sout.close();
	Y[0].~matrix();
	Y[1].~matrix();
}

void lab1() {

	//Generowanie liczb losowych
	double dolna_granica_pkt = -100;
	double gorna_granica_pkt = 100;
	std::uniform_real_distribution<double> unif(dolna_granica_pkt, gorna_granica_pkt);
	std::default_random_engine re;

	//Zapisywanie do pliku
	fstream file;
	//file << std::fixed << setprecision(7);

	//file.open("ekspansja.txt");
	//
	//if (file.good()) {

	//	for (int j = 1; j < 4; j++) {

	//		//Wspolczynnik ekspansji
	//		double wsp = j;

	//		//Optymalizacja 1
	//		for (int i = 0; i < 100; i++) {
	//			//Losowy punkt
	//			double pkt = unif(re);

	//			//Eskpansja
	//			double* p = expansion(f1, pkt, 1, wsp, 10000);

	//			file << pkt << ";" << p[0] << ";" << p[1] << ";" << solution::f_calls << ";";

	//			//Fib
	//			solution min_fib = fib(f1, p[0], p[1], 0.01);

	//			file << min_fib.x(0) << ";" << m2d(f1(min_fib.x(0))) << ";" << solution::f_calls << ";" << "NULL" << ";";

	//			//Lag
	//			solution min_lag = lag(f1, p[0], p[1], 0.01, 0.01, 10000);

	//			file << min_lag.x(0) << ";" << m2d(f1(min_lag.x(0))) << ";" << solution::f_calls << ";" << "NULL\n";
	//		}

	//	}

	//}

	////Bez ekspansji

	//file.close();
	//file.open("bezekspansji.txt");

	////Fib
	//solution min_fib = fib(f1, dolna_granica_pkt, gorna_granica_pkt, 0.01);
	//file << min_fib.x(0) << ";" << m2d(f1(min_fib.x(0))) << ";" << solution::f_calls << ";" << "NULL" << ";";

	////Lag
	//solution min_lag = lag(f1, dolna_granica_pkt, gorna_granica_pkt, 0.01, 0.0001, 10000);
	//file << min_lag.x(0) << ";" << m2d(f1(min_lag.x(0))) << ";" << solution::f_calls << ";" << "NULL\n";

	// przedzial w zadaniu 1 - 100 cm^2 -> 0,0001 - 0,01 m^2, nw jaki ma byæ krok bo zawsze koniec przedzia³u to pocz¹tek + krok
	// zmiana tych podanych w instrukcji parametrów w df1 nie wp³ywa w ogóle na wynik

	random_device R;
	double d = 1, alpha = 2, epsilon = 1e-4, gamma = 1e-4;
	int Nmax = 1000;

	double Da = 100.0 * R() / R.max();
	cout << "Da: " << Da << endl;
	Da = Da * pow(10, -4);


	double* p = expansion(ff1R, Da, d, alpha, Nmax);

	cout << "Przedzial: " << p[0] << " - " << p[1] << endl;

	solution min_fib = fib(ff1R, p[0], p[1], epsilon);
	solution min_lag = lag(ff1R, p[0], p[1], epsilon, gamma, Nmax);

	//min_fib.fit_fun(ff1R);
	//min_lag.fit_fun(ff1R);

	cout << "Wynik z metody Fibbonacciego " << min_fib << endl;
	cout << "Wynik z metody Lagrange'a: " << min_lag << endl;

	file.close();
}

void lab2() {

}

void lab3() {

}

void lab4() {

}

void lab5() {

}

void lab6() {

}
