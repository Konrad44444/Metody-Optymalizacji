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
		lab2();
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
	int Nmax = 1000;
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

	//random_device R;
	//double d = 1, alpha = 2, epsilon = 1e-6, gamma = 1e-12;
	//int Nmax = 1000;
	//solution min_fib;
	//solution min_lag;

	//	double* p = expansion(ff1R, 0.00192143, d, alpha, Nmax);

	//	cout << "Przedzial: " << p[0] << " - " << p[1] << endl;

	//	min_fib = fib(ff1R, p[0], p[1], epsilon);
	//	//min_lag = lag(ff1R, p[0], p[1], epsilon, gamma, Nmax);

	//	min_fib.fit_fun(ff1R);
	//	//min_lag.fit_fun(ff1R);

	//	cout << "Wynik z metody Fibbonacciego " << min_fib << endl;
	//	//cout << "Wynik z metody Lagrange'a: " << min_lag << endl;
	matrix Y0 = matrix(3, new double[3] { 5, 1, 10 });
	matrix ud1;
	matrix* Y = solve_ode(df1, 0, 1, 1000, Y0, ud1, 0.00201729);

	file.open("lag.txt");

	if (file.good()) {
		file << Y[1];
	}
	file.close();

}

void lab2() {
	//random_device R;

	//fstream file1;
	//fstream file2;

	//// symulacja 100 razy z 3 ró¿nymi krokami

	///*file1.open("hj.txt");
	//file2.open("rosen.txt");

	//double s[] = { 0.001, 0.01, 0.005 };
	//
	//for (int i = 0; i < 3; i++) {

	//	double step = s[i];


	//	for (int j = 0; j < 100; j++) {

	//		double x1 = 2.0 * R() / R.max() - 1.0;
	//		double x2 = 2.0 * R() / R.max() - 1.0;
	//
	//		double epsilon = 0.0001;
	//		int Nmax = 10000;
	//
	//		double alphaHJ = 0.5;
	//		double alphaR = 2, beta = 0.5;

	//		matrix x = matrix(2, new double[2] {x1, x2});
	//		matrix s0 = matrix(2, new double[2] {step, step});

	//		solution hj = HJ(f2, x, step, alphaHJ, epsilon, Nmax);
	//		hj.fit_fun(f2);
	//		file1 << x1 << ";" << x2 << ";" << hj.x(0) << ";" << hj.x(1) << ";" << hj.y(0) << ";" << solution::f_calls << ";" << (abs(hj.y(0)) < 0.001 ? "TAK" : "NIE") << "\n";

	//		solution r = Rosen(f2, x, s0, alphaR, beta, epsilon, Nmax);
	//		r.fit_fun(f2);
	//		file2 << r.x(0) << ";" << r.x(1) << ";" << r.y(0) << ";" << solution::f_calls << ";" << (abs(r.y(0)) < 0.001 ? "TAK" : "NIE") << "\n";

	//	}

	//}

	//file1.close();
	//file2.close();*/

	//// do wykresu

	///*double x1 = 0.2661680;
	//double x2 = 0.1869780;

	//double alphaHJ = 0.5;
	//double alphaR = 2, beta = 0.5;
	//double epsilon = 0.0001;
	//int Nmax = 10000;
	//double step = 0.001;

	//matrix x = matrix(2, new double[2] {x1, x2});
	//matrix s0 = matrix(2, new double[2] {step, step});


	//solution hj = HJ(f2, x, step, alphaHJ, epsilon, Nmax);
	//solution r = Rosen(f2, x, s0, alphaR, beta, epsilon, Nmax);*/

	//double k1 = 10.0 * R() / R.max();
	//double k2 = 10.0 * R() / R.max();

	//cout << "k1: " << k1 << ", k2: " << k2 << "\n";

	//double alphaHJ = 0.75;
	//double alphaR = 2.0, beta = 0.7;
	//double epsilon = 0.001;
	//int Nmax = 1000;
	//double step = 0.1;

	//matrix x0 = matrix(2, new double[2]{ k1, k2 });
	//matrix s0 = matrix(2, new double[2]{ step, step });

	////2,8; 3,8
	//
	//solution hj = HJ(ff2R, x0, step, alphaHJ, epsilon, Nmax);
	//hj.fit_fun(ff2R);
	//cout << "HJ:\n" << hj << "\n"; 
	//
	//solution r = Rosen(ff2R, x0, s0, alphaR, beta, epsilon, Nmax);
	//r.fit_fun(ff2R);
	//cout << "Rosen:\n" << r << "\n";

	matrix Y0(2, 1), Yref(2, new double[2] { 3.14, 0 }), X(2, new double[2] { 2.81863, 3.84227 });
	matrix* Y = solve_ode(df2, 0, 0.1, 1000, Y0, Yref, X);

	fstream file;
	file.open("Rosen_sim.txt");

	if (file.good()) {
		file << Y[1];
	}
	file.close();
}

void lab3() {



}

void lab4() {

}

void lab5() {

}

void lab6() {

}
