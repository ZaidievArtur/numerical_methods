#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

double f(double x, double y) {
    return exp(x) / x + y;
}

double exact(double x) {
    return exp(x) * (log(abs(x)) + 1);
}

double kuttaMerson(double a, double y_0, double b, int* n, double (*f)(double x, double y), double epsilon) {
    double k1, k2, k3, k4, k5;
    double x, y, h, error;
    bool condition = false;
    *n = 1;
    y = y_0;
    x = a;

    while (!condition) {
        h = (b - a) / (*n);
        bool stepAccepted = true;

        for (int i = 0; i < (*n); i++) {
            k1 = f(x, y);
            k2 = f(x + h / 3, y + h * k1 / 3);
            k3 = f(x + h / 3, y + h * k1 / 6 + h * k2 / 6);
            k4 = f(x + h / 2, y + h * k1 / 8 + 3 * h * k2 / 8);
            k5 = f(x + h, y + h * k1 / 2 - 3 * h * k3 / 2 + 2 * h * k4);

            double y_new = y + (k1 + 4 * k4 + k5) * h / 6;
            double R = fabs(2 * k1 - 9 * k3 + 8 * k4 - k5) * h / 30;

            if (R <= epsilon) {
                y = y_new;
                x += h;
            }
            else {
                stepAccepted = false;
                break;
            }
        }

        if (stepAccepted) {
            condition = true;
            error = fabs(y - exact(x));
        }
        else {
            (*n) *= 2;
        }
    }

    return error;
}


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
        error.resize(x_vals.size());
        error[i] = fabs(y - exact(x));
    }
    return error;
}

int main() {
    ofstream file("data.txt");

    int i, n = 1;
    double epsilon = 0.1;
    double error;
    double x;
    vector<double> error1;
    vector<double> error2;
    double a = 1;
    double b = 3;
    double y_0 = exp(1);

    vector<double> x_h1, y_h1, x_h2, y_h2;



    error1 = Research(a, b, y_0, 4, f, x_h1, y_h1);
    error2 = Research(a, b, y_0, 8, f, x_h2, y_h2);


    for (size_t i = 0; i < x_h1.size(); i++) {
        file << x_h1[x_h1.size()-i-1] << " " << y_h1[x_h1.size() - i-1] << " " << error1[i] << "\n";
    }
    for (size_t i = 0; i < x_h2.size(); i++) {
        file << x_h2[x_h2.size()-i-1] << " " << y_h2[x_h2.size() - i-1] << " " << error2[i] << "\n";
    }

    for (i = 0; i < 9; i++) {
        /*error*/x = kuttaMerson(a, y_0, b, &n, f, epsilon);
        file << scientific;
        file.precision(15);
        file << /*error*/ x << " " << epsilon << " " << (b - a) / n << "\n";
        epsilon /= 10;

    }

    file.close();
    return 0;
}
