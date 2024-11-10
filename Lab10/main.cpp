
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

const double A = -3.0;
const double B = 3.0;

using namespace std;

double f(double x) {
    return x * x * x * x - 4.2 * x * x * x + 3.5 * x * x - 7 * x - 7.4;
}

double Trapezoid(int n, double previous, double (*f)(double x)) {
    double result = 0.0;

    double h = (B - A) / static_cast<double>(n);

    /*int i = 0;
    double xk = A + i * h;*/
    if (n == 1) {
        result = ((f(A) + f(B)) / 2);
        for (int i = 1; i < n; i++) {
            result += f(A + i * h)*h;

        }
        

    }
    else {
         double result = (previous / 2.0)/ h;
        for (int i = 1; i < n; i += 2) {

            result += f(A + i * h)*h;

        }
       
        
    }
    result *= h;

  
    return result;
}

double Integral(double epsilon, int* N, int* n, double (*f)(double x)) {
    (*n)= 1;
    double next = Trapezoid(*n, 0, f);
    double prev;
    N = 0;

    while (true) {
        prev = next;
        *n *= 2;
        next = Trapezoid(*n, prev, f);
        N++;
        if (fabs(next - prev) / 3 < epsilon) {
            break;
        }
    }
    return next;
}

int main() {
    ofstream file("data.txt");

    int N, n;
    double Iexact = 115.8;
    double epsilon = 0.1;
    double result;

    file << fixed << setprecision(15);
    for (int i = 0; i < 9; i++) {
        result = Integral(epsilon, &N, &n, f);
        file << fabs(Iexact - result) << " " << N << " " << log2(static_cast<double>(B - A) / n) << "\n";
        epsilon /= 10;
    }

    file.close();
    return 0;
}
