//
//  AlDir.cpp
//  somecode9.1
//
//  Created by Михаил on 10.03.2025.
//

#include "AlDir.hpp"


double f_func_2_real(double x, double y) {
    return  sin(PI * x) * sin(PI * y);
}


// Функция для выделения памяти под двумерный массив
double** create_grid(int size) {
    double** grid = new double*[size];
    for (int i = 0; i < size; i++) {
        grid[i] = new double[size]();
    }
    return grid;
}

// Функция для освобождения памяти
void free_grid(double** grid, int size) {
    for (int i = 0; i < size; i++) {
        delete[] grid[i];
    }
    delete[] grid;
}

// Метод прогонки для x-направления (столбцы)
void thomas_algorithm_x(double* a, double* b, double* c, double* d, int n, double** solution, int j) {
    double* c_prime = new double[n+1];
    double* d_prime = new double[n+1];

    c_prime[1] = c[1] / b[1];
    d_prime[1] = d[1] / b[1];

    for (int i = 2; i <= n; i++) {
        double denom = b[i] - a[i] * c_prime[i-1];
        c_prime[i] = c[i] / denom;
        d_prime[i] = (d[i] - a[i] * d_prime[i-1]) / denom;
    }

    solution[n][j] = d_prime[n];
    for (int i = n-1; i >= 1; i--) {
        solution[i][j] = d_prime[i] - c_prime[i] * solution[i+1][j];
    }

    delete[] c_prime;
    delete[] d_prime;
}

// Метод прогонки для y-направления (строки)
void thomas_algorithm_y(double* a, double* b, double* c, double* d, int n, double** solution, int i) {
    double* c_prime = new double[n+1];
    double* d_prime = new double[n+1];

    c_prime[1] = c[1] / b[1];
    d_prime[1] = d[1] / b[1];

    for (int j = 2; j <= n; j++) {
        double denom = b[j] - a[j] * c_prime[j-1];
        c_prime[j] = c[j] / denom;
        d_prime[j] = (d[j] - a[j] * d_prime[j-1]) / denom;
    }

    solution[i][n] = d_prime[n];
    for (int j = n-1; j >= 1; j--) {
        solution[i][j] = d_prime[j] - c_prime[j] * solution[i][j+1];
    }

    delete[] c_prime;
    delete[] d_prime;
}


// Функция ADI
void ADI(double** u_old, double** u_half, double** u_new, double** f) {
    for (int iter = 0; iter < MAX_ITER; iter++) {
        // Полушаг по x
        for (int j = 1; j <= N; j++) {
            double* a = new double[N+1];
            double* b = new double[N+1];
            double* c = new double[N+1];
            double* d = new double[N+1];

            for (int i = 1; i <= N; i++) {
                double x_i = i * h;
                double cos_term = cos(PI * x_i) * cos(PI * x_i);

                a[i] = -1.0;
                b[i] = 2.0 + sigma * h * h;
                c[i] = -1.0;

                d[i] = sigma * h * h * u_old[i][j]
                    + cos_term * (u_old[i][j+1] - 2 * u_old[i][j] + u_old[i][j-1])
                    + h * h * f[i][j];
            }
            thomas_algorithm_x(a, b, c, d, N, u_half, j);

            delete[] a;
            delete[] b;
            delete[] c;
            delete[] d;
        }

        // Полушаг по y
        for (int i = 1; i <= N; i++) {
            double* a = new double[N+1];
            double* b = new double[N+1];
            double* c = new double[N+1];
            double* d = new double[N+1];

            for (int j = 1; j <= N; j++) {
                double x_i = i * h;
                double cos_term = cos(PI * x_i) * cos(PI * x_i);

                a[j] = -cos_term;
                b[j] = 2.0 * cos_term + sigma * h * h;
                c[j] = -cos_term;

                d[j] = sigma * h * h * u_half[i][j]
                    + (u_half[i+1][j] - 2 * u_half[i][j] + u_half[i-1][j])
                    + h * h * f[i][j];
            }
            thomas_algorithm_y(a, b, c, d, N, u_new, i);

            delete[] a;
            delete[] b;
            delete[] c;
            delete[] d;
        }

        // Проверка сходимости
        double max_diff = 0.0;
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                max_diff = std::max(max_diff, abs(u_new[i][j] - u_old[i][j]));
            }
        }

        if (max_diff < EPSILON) {
            std::cout << "Сходимость достигнута на итерации: " << iter << std::endl;
            break;
        }

        std::swap(u_old, u_new);
    }
}

// Функция сохранения результатов в файл
void save_solution(const std::string& filename, double** u_new) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Ошибка открытия файла для записи!" << std::endl;
        return;
    }
    double max_diff = 0.0;
    double u_real = 0.0;
    for (int i = 0; i <= N+1; i++) {
        for (int j = 0; j <= N+1; j++) {
            out << i * h << " " << j * h << " " << u_new[i][j] << std::endl;
            u_real = f_func_2_real(i * h, j*h);
            max_diff = std::max(max_diff, abs(u_new[i][j] - u_real));
            
        }
        out << std::endl;
    }
    std::cout<<"Максимальная разность : "<< max_diff<<std::endl;
    out.close();
}
