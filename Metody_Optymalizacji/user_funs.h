#pragma once

#include"ode_solver.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);
matrix f1(matrix x, matrix = NAN, matrix = NAN);
matrix df1(double, matrix, matrix, matrix);
matrix ff1R(matrix x, matrix ud1, matrix ud2);
matrix f2(matrix x1, matrix x2, matrix ud1);
matrix df2(double, matrix Y, matrix ud1, matrix ud2);
matrix ff2R(matrix x, matrix ud1, matrix ud2);
matrix f3_zewn(matrix x, matrix ud1, matrix ud2);
matrix f3_wewn(matrix x, matrix ud1, matrix ud2);
matrix df3(double, matrix Y, matrix ud1, matrix ud2);
matrix ff3R(matrix x, matrix ud1, matrix ud2);
matrix f4(matrix x1, matrix ud1 = NAN, matrix ud2 = NAN);
matrix f4_grad(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix f4_hess(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);