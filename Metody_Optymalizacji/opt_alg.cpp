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
		//Tu wpisz kod funkcji

		return p;

	} catch (string ex_info) {
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2) {
	
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	
	} catch (string ex_info) {
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2) {
	
	try {
		solution Xopt;
		
		solution A = solution(a);
		solution B = solution(b);
		solution C = solution(((b - a) / 2));
		solution D = solution(0);
		solution prevD = solution(D.x);

		while ((B.x - A.x < epsilon) || (D - prevD) < gamma) {
			double licznik = A.fit_fun(ff) * (B.x * b.x - C.x * C.x) + B.fit_fun(ff) * (C.x * C.x - A.x * A.x) + C.fit_fun(ff) * (A.x * A.x - B.x * B.x);
			double mianownik = A.fit_fun(ff) * (B.x - C.x) + B.fit_fun(ff) * (C.x - A.x) + C.fit_fun(ff) * (A.x - B.x);

			if (mianownik <= 0) throw("Dzielenie przez 0!");

			prevD.x = D.x;
			D.x = 0, 5 * licznik / mianownik;

			if ((A.x < D) && (D < C.x)) {
				
				if (D.fit_fun(ff) < C.fit_fun(ff)) {
					// A.x = A.x;
					C.x = D;
					B.x = C.x;
				} else {
					A.x = D;
					// C.x = C.x;
					// B.x = B.x;
				}

			} else {

				if ((C.x < D.x) && (D.x < B.x)) {

					if (D.fit_fun(ff) < C.fit_fun()) {
						A.x = C.x;
						C.x = D.x;
						//B.x = B.x;
					} else {
						//A.x = A.x;
						//C.x = C.x;
						B.x = D.x;
					}

				} else {
					throw("Error in if condition!");
				}

			}
			
			if (solution::f_calls > Nmax) throw("f_calls > Nmax");

		}

		Xopt = solution(D.x);

		return Xopt;
	
	} catch (string ex_info) {
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2) {
	
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	
	} catch (string ex_info) {
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2) {
	
	try {
		//Tu wpisz kod funkcji

		return XB;
	} catch (string ex_info) {
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
