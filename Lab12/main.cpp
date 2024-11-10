#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

// Функция для вычисления правой части ОДУ
double f(double x, double y) {
    return exp(x) / x + y;
}

// Метод Кутты-Мерсона для решения задачи Коши с адаптивным шагом и сравнением двух последовательных приближений
double kuttaMerson(double a, double y_0, double b, int* n, double (*f)(double x, double y), double epsilon) {
    double k1, k2, k3, k4, k5;
    double x = a, y = y_0, h, y_prev;
    bool condition = false;
    *n = 1;

    while (!condition) {
        h = (b - a) / (*n);
        bool stepAccepted = true;
        double prev_k1 = f(x, y);

        for (int i = 0; i < (*n); i++) {
            k1 = prev_k1;
            k2 = f(x + h / 3, y + h * k1 / 3);
            k3 = f(x + h / 3, y + h * k1 / 6 + h * k2 / 6);
            k4 = f(x + h / 2, y + h * k1 / 8 + 3 * h * k2 / 8);
            k5 = f(x + h, y + h * k1 / 2 - 3 * h * k3 / 2 + 2 * h * k4);

            y_prev = y;
            y += (k1 + 4 * k4 + k5) * h / 6;
            double R = fabs(2 * k1 - 9 * k3 + 8 * k4 - k5) * h / 30;

            if (R <= epsilon) {
                x += h;
                prev_k1 = f(x, y); // Обновляем prev_k1 только если шаг принят

                // Оценка точности через разность двух последовательных значений y
                if (fabs(y - y_prev) < epsilon) {
                    condition = true;
                    break;
                }
            } else {
                stepAccepted = false;
                break;
            }
        }

        if (!stepAccepted) {
            (*n) *= 2; // Уменьшаем шаг
        } else if (fabs(y - y_prev) < epsilon / 10) {
            (*n) = max(1, (*n) / 2); // Увеличиваем шаг, если точность достигнута с запасом
        }
    }

    return fabs(y - y_prev);  // Возвращаем оценку ошибки
}

// Функция для исследования метода Кутты-Мерсона
vector<double> Research(double a, double b, double y_0, int n, double (*f)(double, double), vector<double>& x_vals, vector<double>& y_vals) {
    double h = (b - a) / n;
    double x = a;
    double y = y_0;
    vector<double> error;

    x_vals.push_back(x);
    y_vals.push_back(y);

    for (int i = 0; i < n; i++) {
        double k1 = f(x, y);
        double k2 = f(x + h / 3, y + h * k1 / 3);
        double k3 = f(x + h / 3, y + h * k1 / 6 + h * k2 / 6);
        double k4 = f(x + h / 2, y + h * k1 / 8 + 3 * h * k2 / 8);
        double k5 = f(x + h, y + h * k1 / 2 - 3 * h * k3 / 2 + 2 * h * k4);

        y += (k1 + 4 * k4 + k5) * h / 6;
        x += h;

        x_vals.push_back(x);
        y_vals.push_back(y);

        if (i > 0) {
            error.push_back(fabs(y_vals[i] - y_vals[i - 1])); // Оцениваем ошибку как разность последовательных y
        }
    }
    return error;
}

int main() {
    ofstream file("data.txt");

    int i, n = 1;
    double epsilon = 0.1;
    double error;
    vector<double> error1, error2;
    double a = 1;
    double b = 3;
    double y_0 = exp(1);
    vector<double> x_h1, y_h1, x_h2, y_h2;

    error1 = Research(a, b, y_0, 4, f, x_h1, y_h1);
    error2 = Research(a, b, y_0, 8, f, x_h2, y_h2);

    for (size_t i = 0; i < x_h1.size(); i++) {
        file << x_h1[x_h1.size() - i - 1] << " " << y_h1[x_h1.size() - i - 1] << " " << error1[i] << "\n";
    }
    for (size_t i = 0; i < x_h2.size(); i++) {
        file << x_h2[x_h2.size() - i - 1] << " " << y_h2[x_h2.size() - i - 1] << "\n";
    }

    for (i = 0; i < 9; i++) {
        error = kuttaMerson(a, y_0, b, &n, f, epsilon);
        file << scientific;
        file.precision(15);
        file << error << " " << epsilon << " " << (b - a) / n << "\n";
        epsilon /= 10;
    }

    file.close();
    return 0;
}
