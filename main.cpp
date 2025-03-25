#include "AlDir.hpp"

using namespace std;



// Функция правой части
double f_func(double x, double y) {
    double cX = cos(PI*x);
    return 2*y*(1-y) + 2*cX*cX * x*(1-x);
}


double f_func_2(double x, double y) {
    return PI * PI * sin(PI * x) * sin(PI * y) * (1 + cos(PI * x) * cos(PI * x));
}



int main() {
    // Выделение памяти под сетки
    double** u_old = create_grid(N+2);
    double** u_half = create_grid(N+2);
    double** u_new = create_grid(N+2);
    double** f = create_grid(N+2);

    // Заполнение правой части
    for (int i = 0; i <= N+1; i++) {
        for (int j = 0; j <= N+1; j++) {
            double x = i * h;
            double y = j * h;
            f[i][j] = f_func_2(x, y);
        }
    }

    clock_t start, end;
 
    start = clock();
    // Запуск метода ADI
    ADI(u_old, u_half, u_new, f);

    end = clock();
    cout<<"Время выполнения метода : " << double(end - start) / double(CLOCKS_PER_SEC) << endl;
    // Сохранение результатов
    save_solution("/Users/mihail/Desktop/somecode9.1/somecode9.1/solution.txt", u_new);

    // Освобождение памяти
    free_grid(u_old, N+2);
    free_grid(u_half, N+2);
    free_grid(u_new, N+2);
    free_grid(f, N+2);

    return 0;
}
