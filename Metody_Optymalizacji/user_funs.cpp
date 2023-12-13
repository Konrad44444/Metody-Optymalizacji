#define _USE_MATH_DEFINES
#include"user_funs.h"
#include<math.h>
#include<vector>

matrix ff0T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);
	int n = get_len(Y[0]);
	double teta_max = Y[1](0, 0);
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));
	Y[0].~matrix();
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);
	double m = 1, l = 0.5, b = 0.5, g = 9.81;
	double I = m*pow(l, 2);
	dY(0) = Y(1);
	dY(1) = ((t <= ud2(1))*ud2(0) - m*g*l*sin(Y(0)) - b*Y(1)) / I;
	return dY;
}

matrix f1(matrix x, matrix ud1, matrix ud2)
{
	return -cos(0.1 * m2d(x)) * exp(-pow((0.1 * m2d(x) - 2 * M_PI), 2)) + 0.002 * pow((0.1 * m2d(x)), 2);
}

matrix df1(double t, matrix ud1, matrix ud2, matrix ud3) {
	double A = 0.98, B = 0.63, Pa = 0.7, Pb = 1.0, Va = 5.0, Vb = 1.0, Ta = 90, Tb = 10, Db = 36.5665e-4, G = 9.81, Fin = 0.01, Tin = 10.0;

	
	double FaOUT;
	double FbOUT;
	FaOUT= A * B * m2d(ud3) * sqrt(2 * G * (ud1(0) / Pa)) ;
	FbOUT=  A * B * Db * sqrt(2 * G * (ud1(1) / Pb));

	matrix dY = matrix(3, 1);
	dY(0) = -1 * FaOUT;
	dY(1) = FaOUT + Fin - FbOUT;
	dY(2) = Fin / ud1(1) * (Tin - ud1(2)) + FaOUT / ud1(1) * (Ta - ud1(2));

	return dY;
}

matrix ff1R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(3, new double[3]{ 5, 1, 10 });
	matrix* Y = solve_ode(df1, 0, 1, 1000, Y0, ud1, x);

	int n = get_len(Y[0]);
	double max = Y[1](0, 2);
	for (int i = 1; i < n; i++) {
		if (max < Y[1][2](i))
			max = Y[1][2](i);
	}
	y(0) = abs(max - 50);
	return y;

}

matrix f2(matrix x, matrix ud1, matrix ud2)
{
	return pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * M_PI * x(0)) - cos(2.5 * M_PI * x(1)) + 2;
}

matrix df2(double x, matrix Y, matrix ud1, matrix ud2) {
	double mc = 9.5, mr = 1, l = 0.6, b = 0.5;
	double I = (mr * l * l) / 3 + (mc * l * l);
	matrix dY(2, 1);
	dY(0) = Y(1);
	dY(1) = (ud2(0) * (ud1(0) - Y(0)) + ud2(1) * (ud1(1) - Y(1)) - b * Y(1)) / I;
	return dY;
}

matrix ff2R(matrix x, matrix ud1, matrix ud2) {

	matrix y;
	matrix Y0(2, 1), Yref(2, new double[2]{ 3.14,0 });
	matrix* Y = solve_ode(df2, 0, 0.1, 1000, Y0, Yref, x);
	int n = get_len(Y[0]);

	double temp = 0;
	for (int i = 0; i < n; i++) {

		temp = temp + 10 * pow(Yref(0) - Y[1](i, 0), 2) + pow(Yref(1) - Y[1](i, 1), 2) +
			pow(x(0) * (Yref(0) - Y[1](i, 0)) + x(1) * (Yref(1) - Y[1](i, 1)), 2);

	}
	temp = temp * 0.1;

	y = temp;

	return y;

}

matrix f3_zewn(matrix x, matrix ud1, matrix ud2) {

	matrix y;

	y = (sin(M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2))) / (M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2))));

	if (-x(0) + 1 > 0) {
		y = y + ud2 * pow(-x(0) + 1, 2);
	}

	if (-x(1) + 1 > 0) {
		y = y + ud2 * pow(-x(1) + 1, 2);
	}

	if (norm(x) - ud1 > 0) {
		y = y + ud2 * pow(norm(x) - ud1, 2);
	}

	return y;

}

matrix f3_wewn(matrix x, matrix ud1, matrix ud2) {
	// ud2 - c
	// ud1 - a
	matrix y;

	y = (sin(M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2))) / (M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2))));

	if (-x(0) + 1 > 0) {
		y = 1e10;
	} else  {
		y = y - ud2 / (-x(0) + 1);
	}

	if (-x(1) + 1 > 0) {
		y = 1e10;
	} else {
		y = y - ud2 / (-x(1) + 1);
	}

	if (norm(x) - ud1 > 0) {
		y = 1e10;
	} else {
		y = y - ud2 / (norm(x) - ud1);
	}

	return y;
}

matrix df3(double x, matrix Y, matrix ud1, matrix ud2) {
	double m = 0.6, r = 0.12, y0 = 100;
	double g = 9.81, C = 0.47, ro = 1.2, S = M_PI * r * r, omega = ud1(0);

	double dX = 0.5 * C * ro * S * Y(1) * Y(1);
	double dY = 0.5 * C * ro * S * Y(3) * Y(3);
	double FMx = M_PI * ro * Y(3) * omega * r * r * r;
	double FMy = M_PI * ro * Y(1) * omega * r * r * r;

	matrix DY(4, 1);
	DY(0) = Y(1);
	DY(1) = (-dX - FMx) / m;
	DY(2) = Y(3);
	DY(3) = (-m * g - dY - FMy) / m;

	return DY;
}

matrix ff3R(matrix x, matrix ud1, matrix ud2) {
	double T0 = 0, dT = 0.01, Tend = 7;

	matrix Y0(4, new double[4] {0, x(0), 100, 0});
	matrix X1 = x(1);
	matrix* Y = solve_ode(df3, T0, dT, Tend, Y0, X1);

	int n = get_len(Y[0]);
	int i50 = 0;
	int i0 = 0;

	for (int i = 0; i < n; i++) {
		if (abs(Y[1](i, 2) - 50) < abs(Y[1](i50, 2) - 50)) {
			i50 = i;
		}

		if (abs(Y[1](i, 2)) < abs(Y[1](i0, 2))) {
			i0 = i;
		}
	}

	matrix y = -Y[1](i0, 0);

	if (abs(x(0)) - 10 > 0) {
		y = y + ud2(0) * pow(abs(x(0)) - 10, 2);
	}

	if (abs(x(1)) - 23 > 0) {
		y = y + ud2(0) * pow(abs(x(1)) - 23, 2);
	}

	if (abs(Y[1](i50, 0) - 5) - 0.95 > 0) {
		y = y + ud2(0) * pow(abs(Y[1](i50, 0) - 5) - 0.95, 2);
	}

	return y;
}

matrix f4(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	if (isnan(ud2(0, 0))) {
		y = pow(x(0) + 2 * x(1) - 7, 2) + pow(2 * x(0) + x(1) - 5, 2);
	}
	else {
		y = f4(ud2[0] + x * ud2[1], ud1);
	}
	return y;
}

matrix f4_grad(matrix x, matrix ud1, matrix ud2)
{
	matrix g(2, 1);
	g(0) = 10 * x(0) + 8 * x(1) - 34;
	g(1) = 8 * x(0) + 10 * x(1) - 38;
	return g;
}

matrix f4_hess(matrix x, matrix ud1, matrix ud2)
{
	matrix h(2, 2);
	h(0, 0) = 10;
	h(0, 1) = 8;
	h(1, 0) = 8;
	h(1, 1) = 10;
	return h;
}

matrix ff4R(matrix teta, matrix x, matrix y) {

	matrix J;

	int* n = get_size(y);
	double suma = 0.0;

	for (int i = 0; i < n[1]; i++) {
		suma += y[i](0) * log(h0(teta, x[i])) + (1 - y[i](0)) * log(1 - h0(teta, x[i]));
	}

	suma *= -1;
	suma /= n[1];

	J(0) = suma;

	return J;
}

matrix ff4R_grad(matrix teta, matrix x, matrix y) {

	matrix J(3, 1);

	int* n = get_size(y);

	for (int j = 0; j < get_len(teta); j++) {
		double suma = 0.0;

		for (int i = 0; i < n[1]; i++) {
			suma += (h0(teta, x[i]) - y[i](0)) * x[i](j);
		}
	
		suma /= n[1];

		J(j) = suma;
	}

	return J;
}

double h0(matrix teta, matrix x) {	
	return 1 / (1 + exp(m2d(-1 * trans(teta) * x)));
}


matrix f5(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	
	y = pow(x[0](0), 2) + pow(x[0](1), 2) - cos(2.5 * M_PI * x[0](0)) - cos(2.5 * M_PI * x[0](1)) + 2;

	return y;
}

void bubbleSort(double y[], int index[], int n)
{
	int i, j;
	bool swapped;
	for (i = 0; i < n - 1; i++) {
		swapped = false;
		for (j = 0; j < n - i - 1; j++) {
			if (y[j] > y[j + 1]) {
				swap(y[j], y[j + 1]);
				swap(index[j], index[j + 1]);
				swapped = true;
			}
		}

		if (swapped == false)
			break;
	}
}