#pragma once

#include"ode_solver.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);
matrix df1(double, matrix, matrix = NAN, matrix = NAN);

matrix ff1R(matrix, matrix = NAN, matrix = NAN);

matrix ff1T(matrix, matrix = NAN, matrix = NAN);

matrix ff3T(matrix, matrix = NAN, matrix = NAN);
matrix df3(double, matrix, matrix = NAN, matrix = NAN);
matrix ff3R(matrix, matrix = NAN, matrix = NAN);

matrix gradient(matrix x, matrix ud1, matrix ud2);
matrix hesjan(matrix x, matrix ud1, matrix ud2);
matrix ff4T(matrix x, matrix ud1, matrix ud2);
matrix ff4R(matrix x, matrix ud1, matrix ud2);
matrix gf(matrix x, matrix ud1, matrix ud2);

matrix ff5T(matrix x, matrix = NAN, matrix = NAN);

matrix ff5T_f1(matrix x, matrix = NAN, matrix = NAN);
matrix ff5T_f2(matrix x, matrix = NAN, matrix = NAN);