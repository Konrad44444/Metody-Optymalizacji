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

		/*fstream file;
		file.open("hj.txt");*/

		solution::clear_calls();

		solution XB, _XB, X;
		XB.x = x0;
		XB.fit_fun(ff);

		while(true) {

			X = HJ_trial(ff, XB, s);

			if (X.y < XB.y) {

				while(true) {

					//file << X.x(0) << ";" << X.x(1) << "\n";

					_XB = XB;
					XB = X;
					X.x = 2.0 * XB.x - _XB.x;
					X = HJ_trial(ff, XB, s);
					X.fit_fun(ff);

					if (solution::f_calls > Nmax) break;
					if (X.y >= XB.y) break;

				};
			}
			else {
				s *= alpha;
			}

			if (solution::f_calls < Nmax) break;
			if (s < epsilon) break;

		}

		//file.close();

		return XB;
	
	} catch (string ex_info) {
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2) {
	
	try {
		int* n = get_size(XB.x);
		matrix D(n[0], n[0]);
		for (int i = 0; i < n[0]; i++) {
			D(i, i) = 1;
		}
		solution X;

		for (int j = 0; j < n[0]; j++) {

			X.x = (XB.x + s * D[j]);
			X.fit_fun(ff);
			if (X.y < XB.y) {
				XB = X;
			}
			else {
				X.x = (XB.x - s * D[j]);
				X.fit_fun(ff);
				if (X.y < XB.y) {
					XB = X;
				}
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
		//fstream file;
		//file.open("rosen.txt");

		solution::clear_calls();

		solution Xopt;
		//Tu wpisz kod funkcji

		int* n = get_size(x0);
		matrix l(n[0], 1), p(n[0], 1), s(s0); // l -> lambda
		matrix D(n[0], n[0]); // macierz kierunku
		for (int i = 0; i < n[0]; i++) {
			D(i, i) = 1;
		}

		solution XB, XBt;
		XB.x = x0;
		XB.fit_fun(ff);

		double max_s;

		do {

			for (int i = 0; i < n[0]; i++) {

				//file << XB.x(1) << ";" << XB.x(1) << "\n";

				XBt.x = XB.x + s(i) * D[i];
				XBt.fit_fun(ff);

				if (XBt.y(0) < XB.y(0)) {
					XB = XBt;
					l(i) += s(i);
					s(i) *= alpha;
				} else {
					s(i) *= (-beta);
					p(i)++;
				}
			}

			bool change = true;
			for (int i = 0; i < n[0]; i++) {
				if (l(i) != 0 || p(i) != 0) {
					change = false;
					break;
				}
			}

			// zmiana bazy kierunków
			if (change) {
				matrix Q(n[0], n[0]), v(n[0], 1);
				for (int i = 0; i < n[0]; ++i) {
					for (int j = 0; j <= i; ++j) {
						Q(i, j) = l(i);
					}
				}

				Q = D * Q;
				v = Q[0] / norm(Q[0]);
				D.set_col(v, 0);

				for (int i = 1; i < n[0]; i++) {
					matrix temp(n[0], 1);

					for (int j = 0; j < i; j++) {
						temp = temp + trans(Q[i]) * D[j] * D[j];
					}

					v = Q[i] - temp / norm(Q[i] - temp);
					D.set_col(v, i);
				}

				l = matrix(n[0], 1);
				p = matrix(n[0], 1);
				s = s0;
			}

			max_s = abs(s(0));

			for (int i = 1; i < n[0]; ++i) {
				if (max_s < abs(s(i))) max_s = abs(s(i));
			}

			if (solution::f_calls > Nmax) break;

		} while (max_s > epsilon);

		//file.close();

		Xopt.x = XB.x;

		return Xopt;
	
	} catch (string ex_info) {
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2) {
	
	try {
		double alpha = 1, beta = 0.5, gamma = 2, delta = 0.5, s = 0.5, C = c;
		solution X(x0), X1;
		while (true)
		{
			X1 = sym_NM(ff, X.x, s, alpha, beta, gamma, delta, epsilon, Nmax, ud1, C);
			if (norm(X.x - X1.x) < epsilon || solution::f_calls > Nmax)
				return X1;
			C *= dc;
			X = X1;
		}
		return X;
	} catch (string ex_info) {
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2) {
	
	try {

		int n = get_len(x0);
		matrix D = ident_mat(n);
		int N = n + 1;//liczba wierzcjo�k�w symlpeksu
		solution* S = new solution[N];//sympleks
		S[0].x = x0;
		S[0].fit_fun(ff, ud1, ud2);
		for (int i = 1; i < N; ++i)
		{
			S[i].x = S[0].x + s * D[i - 1];
			S[i].fit_fun(ff, ud1, ud2);
		}
		solution PR, PE, PN;//reflection,expansion,zaw�enie
		matrix pc;
		int i_min, i_max;
		while (true)
		{
			i_min = i_max = 0;
			for (int i = 1; i < N; ++i)
			{
				if (S[i_min].y > S[i].y)
					i_min = i;
				if (S[i_max].y < S[i].y)
					i_max = i;
			}

			pc = matrix(n, 1);
			for (int i = 0; i < N; ++i)
				if (i != i_max)
					pc = pc + S[i].x;
			pc = pc / (N - 1);
			PR.x = pc + alpha * (pc - S[i_max].x);
			PR.fit_fun(ff, ud1, ud2);
			if (S[i_min].y <= PR.y && PR.y < S[i_max].y)
				S[i_max] = PR;
			else if (PR.y < S[i_min].y)
			{
				PE.x = pc + gamma * (PR.y - pc);
				PE.fit_fun(ff, ud1, ud2);
				if (PR.y <= PE.y)
					S[i_max] = PR;
				else
					S[i_max] = PE;
			}
			else
			{
				PN.x = pc + beta * (S[i_max].x - pc);
				PN.fit_fun(ff, ud1, ud2);
				if (PN.y < S[i_max].y)
					S[i_max] = PN;
				else
				{
					for (int i = 0; i < N; ++i)
						if (i != i_min)
						{
							S[i].x = delta * (S[i].x + S[i_min].x);
							S[i].fit_fun(ff, ud1, ud2);
						}
				}
			}
			double max_s = norm(S[0].x - S[i_min].x);
			for (int i = 1; i < N; ++i)
				if (max_s < norm(S[i].x - S[i_min].x))
					max_s = norm(S[i].x - S[i_min].x);
			if (max_s < epsilon || solution::f_calls > Nmax)
				return S[i_min];
		}

	} catch (string ex_info) {
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

// najszybszy spadek
solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2) {
	
	try {
		//Tu wpisz kod funkcji
		int n = get_len(x0);
		
		solution X, X1;
		X.x = x0;

		matrix d(n, 1), P(n, 2), limits = ud1;
		solution h;
		double b;

		while (true) {
			X.grad(gf, ud1, ud2);
			d = -X.g;
			P.set_col(X.x, 0);
			P.set_col(d, 1);

			if (h0 < 0) {
				b = compute_b(X.x, d, limits);
				h = golden(ff, 0, b, epsilon, Nmax, ud1, P);
				X1.x = X.x + h.x * d;
			}
			else {
				X1.x = X.x + h0 * d;
			}
			
			if (norm(X1.x - X.x) < epsilon || solution::f_calls + solution::g_calls > Nmax) {
				X1.fit_fun(ff, ud1, ud2);
				return X1;
			}

			X = X1;
			
		}
	
	} catch (string ex_info) {
		throw ("solution SD(...):\n" + ex_info);
	}
}

// gradienty sprzężone
solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2) {
	
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		int n = get_len(x0);
		solution X, X1;
		X.x = x0;
		matrix d(n, 1), P(n, 2), limits = ud1;
		solution h;
		double b, beta;
		X.grad(gf, ud1, ud2);
		d = -X.g;

		while (true) {
			P.set_col(X.x, 0);
			P.set_col(d, 1);

			if (h0 < 0) {
				b = compute_b(X.x, d, limits);
				h = golden(ff, 0, b, epsilon, Nmax, ud1, P);
				X1.x = X.x + h.x * d;
			}
			else {
				X1.x = X.x + h0 * d;
			}

			if (norm(X1.x - X.x) < epsilon || solution::f_calls + solution::g_calls > Nmax || norm(X1.g) < epsilon) {
				X1.fit_fun(ff, ud1, ud2);
				return X1;
			}

			X1.grad(gf, ud1, ud2);
			beta = pow(norm(X1.g), 2) / pow(norm(X.g), 2);
			d = -X1.g + beta * d;
			X = X1;
		}

		return Xopt;
	
	} catch (string ex_info) {
		throw ("solution CG(...):\n" + ex_info);
	}
}

// metoda Newtona
solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2) {
	
	try {
		solution Xopt;
		//Tu wpisz kod funkcji
		int n = get_len(x0);
		solution X, X1;
		X.x = x0;
		matrix d(n, 1), P(n, 2), limits = ud1;
		solution h;
		double b;

		while (true) {
			X.grad(gf, ud1, ud2);
			X.hess(Hf, ud1, ud2);
			d = -inv(X.H) * X.g;
			P.set_col(X.x, 0);
			P.set_col(d, 1);

			if (h0 < 0) {
				b = compute_b(X.x, d, limits);
				h = golden(ff, 0, b, epsilon, Nmax, ud1, P);
				X1.x = X.x + h.x * d;
			}
			else {
				X1.x = X.x + d;
			}

			if (norm(X1.x - X.x) < epsilon || (solution::f_calls + solution::g_calls) > Nmax ||	norm(X1.g) < epsilon ||	det(X1.H) == 0)	{
				X1.fit_fun(ff, ud1, ud2);
				return X1;
			}

			X = X1;
		}

		return Xopt;
	
	} catch (string ex_info) {
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2) {
	
	try {
		solution Xopt;
		//Tu wpisz kod funkcji
		double alpha = (sqrt(5) - 1) / 2;
		solution A, B, C, D;

		matrix AB(2, 1);
		AB(0) = a;
		AB(1) = b;
		A.x = a;
		B.x = b;

		matrix CD(2, 1);
		CD(0) = AB(1) - alpha * (AB(1) - AB(0));
		CD(1) = AB(0) + alpha * (AB(1) - AB(0));

		C.x = CD(0);
		C.fit_fun(ff, ud1, ud2);
		D.x = CD(1);
		D.fit_fun(ff, ud1, ud2);

		while (true) {
			if (C.y < D.y) {
				AB = D.x;
				D = C;
				C.x = B.x - alpha * (B.x - A.x);
				C.fit_fun(ff, ud1, ud2);
			}
			else {
				A = C;
				C = D;
				D.x = A.x + alpha * (B.x - A.x);
				D.fit_fun(ff, ud1, ud2);
			}

			if ((B.x - A.x) < epsilon || solution::f_calls > Nmax) {
				A.x = (A.x + B.x) / 2.0;
				A.fit_fun(ff, ud1, ud2);
				return A;
			}
		}

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

double compute_b(matrix x, matrix d, matrix limits)
{
	int* n = get_size(x);
	double b = 1e9, bi;
	for (int i = 0; i < n[0]; ++i)
	{
		if (d(i) == 0)
			bi = 1e9;
		else if (d(i) > 0)
			bi = (limits(i, 1) - x(i)) / d(i);
		else
			bi = (limits(i, 0) - x(i)) / d(i);
		if (b > bi)
			b = bi;
	}
	return b;
}