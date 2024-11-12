#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

int countAdamsBashforth = 0;

double f(double x, double y) {
    return exp(x) / x + y;
}

double exact(double x) {
    return exp(x) * (log(abs(x)) + 1);
}

double kuttaMersonSetup(double a, double y0, double h, double y[4], double epsilon) {
    int i = 0;
    y[0] = y0;
    double x = a;

    while (i < 3) {
        double k1 = f(x, y0);
        double k2 = f(x + h / 3, y0 + h * k1 / 3);
        double k3 = f(x + h / 3, y0 + h * k1 / 6 + h * k2 / 6);
        double k4 = f(x + h / 2, y0 + h * k1 / 8 + 3 * h * k2 / 8);
        double k5 = f(x + h, y0 + h * k1 / 2 - 3 * h * k3 / 2 + 2 * h * k4);

        if (fabs(2 * k1 - 9 * k3 + 8 * k4 - k5) * h / 30 > epsilon) {
            h /= 2;
            continue;
        }

        y0 += (k1 + 4 * k4 + k5) * h / 6;
        x += h;
        y[++i] = y0;
    }

    return x;
}

double AdamsBashforth4(double a, double y0, double b, double h, double epsilon) {
    double y[4];
    double fxy[4];
    double x = kuttaMersonSetup(a, y0, h, y, epsilon);

    for (int i = 0; i < 4; i++) {
        fxy[i] = f(x - (3 - i) * h, y[i]);
    }

    double y_temp = y[3];
    while (x + h <= b) {
        y_temp = y[3] + (55 * fxy[3] - 59 * fxy[2] + 37 * fxy[1] - 9 * fxy[0]) * h / 24;
        x += h;

        if (fabs(y_temp - y[3]) / 15 > epsilon) {
            h /= 2;
        } else {
            for (int j = 0; j < 3; j++) {
                y[j] = y[j + 1];
                fxy[j] = fxy[j + 1];
            }
            y[3] = y_temp;
            fxy[3] = f(x, y_temp);
        }
    }

    return y_temp;
}

void saveData(const string& filename, const vector<double>& steps, const vector<int>& calls_adams) {
    ofstream file(filename);
    for (size_t i = 0; i < steps.size(); ++i) {
        file << steps[i] << "," << calls_adams[i] << "\n";
    }

    file.close();
}

int main() {
    double a = 1.0;
    double b = 3.0;
    double y0 = 0.0;
    double epsilon = 1e-3;
    vector<double> steps;
    vector<int> calls_adams;

    for (double h = 0.1; h >= 1e-6; h /= 2) {
        countAdamsBashforth = 0;

        AdamsBashforth4(a, y0, b, h, epsilon);

        steps.push_back(h);
        calls_adams.push_back(countAdamsBashforth);
    }

    saveData("adams_comparison.csv", steps, calls_adams);

    return 0;
}
