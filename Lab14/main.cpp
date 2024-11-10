#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#pragma warning(disable: 4996)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <algorithm>
#include <fstream>

using namespace std;

const double a = 0, b = 1;

double p(double x) {
    return 1;
}

double q(double x) {
    return 1 + pow(sin(x), 2);
}

double r(double x) {
    return pow(cos(x), 2);
}

double f(double x) {
    return 3 * exp(x);
}

double exact(double x) {
    return exp(x);
}

double* get_grid(double a, double b, int n) {
    double* grid = (double*)malloc((n + 1) * sizeof(double));
    for (int k = 0; k <= n; k++) {
        grid[k] = a + (double)(b - a) / n * k;
    }
    return grid;
}

double* get_refined_grid(double a, double b, int n) {
    double* grid = (double*)malloc((n + 1) * sizeof(double));
    for (int k = 0; k <= n; k++) {
        double t = (double)k / n;
        grid[k] = a + (b - a) * 0.5 * (1 - cos(M_PI * t));
    }
    return grid;
}

double** createArray(int n, int m) {
    double** array = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        array[i] = (double*)malloc(m * sizeof(double));
    }
    return array;
}

void freeArray(double** array, int n) {
    for (int i = 0; i < n; i++) {
        free(array[i]);
    }
    free(array);
}

double** system_equations(double h, int n, double* grid) {
    double** Y = createArray(n + 1, n + 2);
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= n + 1; j++) {
            Y[i][j] = 0;
        }
    }

    Y[0][0] = 1;
    Y[0][n + 1] = 1;

    for (int i = 1; i <= n - 1; i++) {
        Y[i][i - 1] = p(grid[i]) - q(grid[i]) / 2.0 * h;
        Y[i][i] = pow(h, 2) * r(grid[i]) - 2 * p(grid[i]);
        Y[i][i + 1] = p(grid[i]) + q(grid[i]) / 2.0 * h;
        Y[i][n + 1] = pow(h, 2) * f(grid[i]);
    }

    Y[n][n - 1] = -1 / h;
    Y[n][n] = 1 + 1 / h;
    Y[n][n + 1] = 2 * exp(1);

    double coefficient1 = -Y[1][0] / Y[0][0];
    double coefficient2 = -Y[n - 1][n] / Y[n][n];

    for (int i = 0; i <= n + 1; i++) {
        Y[1][i] += Y[0][i] * coefficient1;
        Y[n - 1][i] += Y[n][i] * coefficient2;
    }

    return Y;
}

double* thomas_algorithm(int n, double** matrix) {
    double* x = (double*)malloc((n + 1) * sizeof(double));
    double* lower_diagonal = (double*)malloc((n - 1) * sizeof(double));
    double* main_diagonal = (double*)malloc((n - 1) * sizeof(double));
    double* upper_diagonal = (double*)malloc((n - 1) * sizeof(double));
    double* solution = (double*)malloc((n - 1) * sizeof(double));
    double* d = (double*)malloc((n - 1) * sizeof(double));
    double* l = (double*)malloc((n - 1) * sizeof(double));

    for (int i = 1; i <= n - 1; i++) {
        lower_diagonal[i - 1] = matrix[i][i - 1];
        main_diagonal[i - 1] = matrix[i][i];
        upper_diagonal[i - 1] = matrix[i][i + 1];
        solution[i - 1] = matrix[i][n + 1];
    }

    double di_old = 1, li_old = 1;
    for (int i = 0; i <= n - 2; i++) {
        d[i] = -upper_diagonal[i] / (di_old * lower_diagonal[i] + main_diagonal[i]);
        l[i] = (solution[i] - lower_diagonal[i] * li_old) / (lower_diagonal[i] * di_old + main_diagonal[i]);
        di_old = d[i];
        li_old = l[i];
    }

    x[n - 1] = l[n - 2];
    for (int i = n - 2; i >= 1; i--) {
        x[i] = d[i - 1] * x[i + 1] + l[i - 1];
    }
    x[0] = (matrix[0][n + 1] - matrix[0][n] * x[n]) / matrix[0][0];
    x[n] = (matrix[n][n + 1] - matrix[n][n - 1] * x[n - 1]) / matrix[n][n];

    free(lower_diagonal);
    free(main_diagonal);
    free(upper_diagonal);
    free(solution);
    free(d);
    free(l);

    return x;
}

double** system_equations_with_accuracy(double a, double b, int& n, double eps, double*& grid) {
    double h = (b - a) / n;
    double max_error;

    do {
        double** Y = createArray(n + 1, n + 2);
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= n + 1; j++) {
                Y[i][j] = 0;
            }
        }

        Y[0][0] = 1;
        Y[0][n + 1] = 1;

        for (int i = 1; i <= n - 1; i++) {
            Y[i][i - 1] = p(grid[i]) - q(grid[i]) / 2.0 * h;
            Y[i][i] = pow(h, 2) * r(grid[i]) - 2 * p(grid[i]);
            Y[i][i + 1] = p(grid[i]) + q(grid[i]) / 2.0 * h;
            Y[i][n + 1] = pow(h, 2) * f(grid[i]);
        }

        Y[n][n - 1] = -1 / h;
        Y[n][n] = 1 + 1 / h;
        Y[n][n + 1] = 2 * exp(1);

        double coefficient1 = -Y[1][0] / Y[0][0];
        double coefficient2 = -Y[n - 1][n] / Y[n][n];

        for (int i = 0; i <= n + 1; i++) {
            Y[1][i] += Y[0][i] * coefficient1;
            Y[n - 1][i] += Y[n][i] * coefficient2;
        }

        double* solution = thomas_algorithm(n, Y);
        max_error = 0;
        for (int i = 0; i <= n; i++) {
            double error = fabs(solution[i] - exact(grid[i]));
            if (error > max_error) {
                max_error = error;
            }
        }

        if (max_error < eps) {
            free(solution);
            return Y;
        }
        else {
            n *= 2;
            h /= 2;
            grid = (double*)realloc(grid, (n + 1) * sizeof(double));
            for (int i = 0; i <= n; i++) {
                grid[i] = a + i * h;
            }
            freeArray(Y, n / 2 + 1);
            free(solution);
        }
    } while (true);
}

int main() {
    int n = 5;
    double* grid = get_grid(a, b, n);
    double h = grid[1] - grid[0];
    double** table_coeff = system_equations(h, n, grid);
    double* solution = thomas_algorithm(n, table_coeff);
    ofstream file5("nodes5.txt");

    for (int i = 0; i <= n; i++) {
        file5 << grid[i] << " " << solution[i] << " " << exact(grid[i]) << " " << fabs(solution[i] - exact(grid[i])) << endl;
    }
    file5.close();
    freeArray(table_coeff, n + 1);
    free(solution);

    n = 10;
    grid = get_grid(a, b, n);
    h = grid[1] - grid[0];
    table_coeff = system_equations(h, n, grid);
    solution = thomas_algorithm(n, table_coeff);
    ofstream file10("nodes10.txt");
    for (int i = 0; i <= n; i++) {
        file10 << grid[i] << " " << solution[i] << " " << exact(grid[i]) << " " << fabs(solution[i] - exact(grid[i])) << endl;
    }
    file10.close();
    freeArray(table_coeff, n + 1);
    free(solution);

    ofstream ferror("error.txt");
    double error = 0;
    double eps = 1e-1;
    for (int k = 1; k < 5; k++) {
        error = 0;
        grid = get_grid(a, b, n);
        h = grid[1] - grid[0];
        table_coeff = system_equations_with_accuracy(0, 1, n, eps, grid);
        solution = thomas_algorithm(n, table_coeff);

        for (int i = 0; i <= n; i++) {
            error = fmax(error, fabs(solution[i] - exact(grid[i])));
        }
        ferror << eps << " " << error << endl;
        freeArray(table_coeff, n + 1);
        free(solution);
        free(grid);
        eps /= 10;
    }
    ferror.close();

    n = 10;
    grid = get_refined_grid(a, b, n);
    h = grid[1] - grid[0];
    table_coeff = system_equations(h, n, grid);
    solution = thomas_algorithm(n, table_coeff);
    ofstream file_refined("nodes_refined.txt");

    for (int i = 0; i <= n; i++) {
        file_refined << grid[i] << " " << solution[i] << " " << exact(grid[i]) << " " << fabs(solution[i] - exact(grid[i])) << endl;
    }
    file_refined.close();
    freeArray(table_coeff, n + 1);
    free(solution);
    free(grid);

    return 0;
}
