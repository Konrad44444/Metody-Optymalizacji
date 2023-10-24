#include"opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2) {
	try	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
			Xopt.fit_fun(ff, ud1, ud2);
			if (Xopt.y < epsilon)
			{
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = 0;
				break;
			}
		}
		return Xopt;
	} catch (string ex_info) {
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2) {

	try {
		double* p = new double[2]{ 0,0 };

		solution::clear_calls();
		
		solution X0(x0);
		solution X1(x0 + d);

		X0.fit_fun(ff, ud1, ud2);
		X1.fit_fun(ff, ud1, ud2);

		if (X0.y(0) == X1.y(0))
		{
			return p = new double[2]{ X0.x(0), X1.x(0) };
		}

		if (X1.y(0) > X0.y(0)) {
		
			d *= -1;
			X1.x(0) = X0.x(0) + d;
			X1.fit_fun(ff, ud1, ud2);
			if (X1.y(0) >= X0.y(0))
			{
				return p = new double[2]{ X1.x(0), X0.x(0)-d };
			}
		}

		solution X2;	

		int i = 1;
		while(true) {

			X2.x(0) = x0 + pow(alpha, i) * d;
			X2.fit_fun(ff, ud1, ud2);
			X1.fit_fun(ff, ud1, ud2);

			if (solution::f_calls > Nmax || X1.y(0) <= X2.y(0))
				break;

			X0.x = X1.x;
			X1.x = X2.x;
			i++;

		}

		if(d > 0)
			return p = new double[2]{ X0.x(0), X2.x(0)};

		return p = new double[2]{ X2.x(0), X0.x(0) };

	} catch (string ex_info) {
		throw ("double* expansion(...):\n" + ex_info);
	}
}

int fib_seq(int n) {
	return (1 / sqrt(5)) * (pow(((1 + sqrt(5)) / (2)), n) - pow(((1 - sqrt(5)) / (2)), n));
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2) {
	
	solution::clear_calls();

	try {
		int k = 0;
		int fi = 0;

		while (true) {
			fi = fib_seq(k);
			if (fi > (b - a) / epsilon) {
				break;
			}
			else {
				k++;
			}
		}

		solution A(a);
		solution B(b);
		solution C(B.x(0) - ((double)fib_seq(k - 1) / (double)fib_seq(k) * (B.x(0) - A.x(0))));
		solution D(A.x(0) + B.x(0) - C.x(0));

		int i = 0;
		while (i <= k - 3) {

			C.fit_fun(ff, ud1, ud2);
			D.fit_fun(ff, ud1, ud2);

			if (C.y(0) < D.y(0)) {
				B.x = D.x(0);
			}
			else {
				A.x = C.x(0);
			}

			C.x(0) = B.x(0) - ((double)fib_seq(k - i - 2) / (double)fib_seq(k - i - 1)) * (B.x(0) - A.x(0));
			D.x(0) = A.x(0) + B.x(0) - C.x(0);

			i++;
		}

		return C;

	}
	catch (string ex_info) {
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2) {
	
	solution::clear_calls();

	try {
		solution A(a);
		solution B(b);
		solution C(((a + b) / 2));
		solution D(DBL_MAX);
		solution prevD(a);

		A.fit_fun(ff, ud1, ud2);
		B.fit_fun(ff, ud1, ud2);
		C.fit_fun(ff, ud1, ud2);

		do {
			double licznik = A.y(0) * (B.x(0) * B.x(0) - C.x(0) * C.x(0)) + B.y(0) * (C.x(0) * C.x(0) - A.x(0) * A.x(0)) + C.y(0) * (A.x(0) * A.x(0) - B.x(0) * B.x(0));
			double mianownik = A.y(0) * (B.x(0) - C.x(0)) + B.y(0) * (C.x(0) - A.x(0)) + C.y(0) * (A.x(0) - B.x(0));

			if (mianownik <= 0) break;

			prevD.x(0) = D.x(0);
			D.x(0) = 0.5 * (licznik / mianownik);
			D.fit_fun(ff, ud1, ud2);

			if ((A.x(0) < D.x(0)) && (D.x(0) < C.x(0))) {

				if (D.y(0) < C.y(0)) {
					C.x(0) = D.x(0);
					B.x(0) = C.x(0);
					B.fit_fun(ff, ud1, ud2);
					C.fit_fun(ff, ud1, ud2);
				}
				else {
					A.x(0) = D.x(0);
					A.fit_fun(ff, ud1, ud2);
				}

			} else {

				if ((C.x(0) < D.x(0)) && (D.x(0) < B.x(0))) {

					if (D.y(0) < C.y(0)) {
						A.x(0) = C.x(0);
						C.x(0) = D.x(0);
						A.fit_fun(ff, ud1, ud2);
						C.fit_fun(ff, ud1, ud2);
					}
					else {
						B.x(0) = D.x(0);
						B.fit_fun(ff, ud1, ud2);
					}

				}
				else {
					break;
				}

			}

			if (solution::f_calls > Nmax) break;

		} while ((B.x(0) - A.x(0) >= epsilon) || abs((D.x(0) - prevD.x(0))) >= gamma);

		return D.x(0);

	}
	catch (string ex_info) {
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2) {
	
	try {
		solution Xopt;
		solution XB(x0);
		solution::clear_calls();
		do {

			solution X = HJ_trial(ff, XB, s);
			XB.fit_fun(ff);
			X.fit_fun(ff);

			if (X.y(0) < XB.y(0)) {

				do {
					solution _XB(XB.x(0));
					XB.x(0) = X.x(0);
					X.x(0) = 2 * XB.x(0) - _XB.x(0);
					X = HJ_trial(ff, XB, s);
					
					if (solution::f_calls > Nmax) break;

					X.fit_fun(ff);
					XB.fit_fun(ff);
				} while (X.y(0) >= XB.y(0));

				X.x(0) = XB.x(0);

			}
			else {

				s = alpha * s;

			}

			if (solution::f_calls > Nmax) break;

		} while (s < epsilon);

		XB.fit_fun(ff);
		Xopt = XB;
		return Xopt;
	
	} catch (string ex_info) {
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2) {
	
	try {

		int j = 1;

		for (int j = 1; j < 2; j++) {

			solution X1(XB.x(0) + s * exp(j));
			solution X2(XB.x(0) - s * exp(j));
			X1.fit_fun(ff);
			X2.fit_fun(ff);

			if (X1.y(0) < XB.y(0)) {
				XB = X1;
			}
			else if (X2.y(0) < XB.y(0)) {
				XB = X2;
			}

		}

		return XB;
	}
	catch (string ex_info) {
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2) {
	
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	
	} catch (string ex_info) {
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2) {
	
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	} catch (string ex_info) {
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2) {
	
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	} catch (string ex_info) {
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2) {
	
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	
	} catch (string ex_info) {
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2) {
	
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	
	} catch (string ex_info) {
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2) {
	
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	
	} catch (string ex_info) {
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2) {
	
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	
	} catch (string ex_info) {
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2) {
	
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	} catch (string ex_info) {
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2) {
	
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	
	} catch (string ex_info) {
		throw ("solution EA(...):\n" + ex_info);
	}
}
