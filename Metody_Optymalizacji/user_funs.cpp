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

	y(0) = 0;
	for (int i = 0; i < n; i++) {

		y(0) = y(0) + 10 * pow(Yref(0) - Y[1](i, 0), 2) + pow(Yref(1) - Y[1](i, 1), 2) +
			pow(x(0) * (Yref(0) - Y[1](i, 0)) + x(1) * (Yref(1) - Y[1](i, 1)), 2);
		y(0) = y(0) * 0.1;

	}

	return y;

}
