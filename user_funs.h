#pragma once

#include"ode_solver.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);
matrix ff1T(matrix, matrix = NAN, matrix = NAN);
matrix ff1R(matrix, matrix = NAN, matrix = NAN);
matrix df1(double, matrix, matrix = NAN, matrix = NAN);
matrix ff2T(matrix, matrix = NAN, matrix = NAN);
matrix ff2R(matrix, matrix = NAN, matrix = NAN);
matrix df2(double, matrix, matrix = NAN, matrix = NAN);
matrix ff3T(matrix, matrix = NAN, matrix = NAN);
matrix ff3T_zewn(matrix, matrix = NAN, matrix = NAN);
matrix ff3T_wewn(matrix, matrix = NAN, matrix = NAN);
matrix ff3R(matrix, matrix = NAN, matrix = NAN);
matrix df3(double, matrix, matrix = NAN, matrix = NAN);
matrix ff4T(matrix, matrix = NAN, matrix = NAN);
matrix gf4T(matrix, matrix = NAN, matrix = NAN);
matrix Hf4T(matrix, matrix = NAN, matrix = NAN);
matrix ff4R_cost(matrix, matrix = NAN, matrix = NAN);
matrix gf4R_grad(matrix, matrix = NAN, matrix = NAN);

// Lab 5 - Multi-criteria optimization
matrix ff5_f1(matrix, matrix = NAN, matrix = NAN);
matrix ff5_f2(matrix, matrix = NAN, matrix = NAN);
matrix ff5T(matrix, matrix = NAN, matrix = NAN);
matrix ff5R(matrix, matrix = NAN, matrix = NAN);