//
//  AlDir.hpp
//  somecode9.1
//
//  Created by Михаил on 10.03.2025.
//

#ifndef AlDir_hpp
#define AlDir_hpp

#include <iostream>
#include <cmath>
#include <algorithm>
#include <fstream>

const double EPSILON = 1e-10;
const int MAX_ITER = 49999;
const double PI = 3.14159265358979323846;

const int N = 99;
const double h = 1.0 / (N + 1);
const double sigma = 2.0 / (h * h);


double** create_grid(int size);

void free_grid(double** grid, int size);

void thomas_algorithm_x(double* a, double* b, double* c, double* d, int n, double** solution, int j);

void thomas_algorithm_y(double* a, double* b, double* c, double* d, int n, double** solution, int i);

void ADI(double** u_old, double** u_half, double** u_new, double** f);

void save_solution(const std::string& filename, double** u_new);




#endif /* AlDir_hpp */
