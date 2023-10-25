#pragma once

#include"ode_solver.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);
matrix f1(matrix x, matrix = NAN, matrix = NAN);
matrix df1(double, matrix, matrix, matrix);
matrix ff1R(matrix x, matrix ud1, matrix ud2);
matrix f2(matrix x1, matrix x2, matrix ud1);