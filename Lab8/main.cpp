
#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>

#define PI 3.14159265358979323846
using namespace std;
double f1(double x) {
    return x - sin(x) - 0.25;
}

double f2(double x) {
    double sign;
    if (x > 0)
        sign = 1;
    else if (x == 0)
        sign = 0;
    else
        sign = -1;

    return pow(x, 5) + 0.4 * sign * pow(x, 4) + 2;
}
double* genUniformNodes(double start, double end, int nodeCount) {
    double* nodes = new double[nodeCount];
    double step = (end - start) / (nodeCount - 1);
    for (int i = 0; i < nodeCount; ++i) {
        nodes[i] = start + i * step;
    }
    return nodes;
}

double* genChebNodes(double start, double end, int nodeCount) {
    double* nodes = new double[nodeCount];
    for (int i = 0; i < nodeCount; ++i) {
        nodes[i] = (start + end) / 2 + cos(PI * (2 * i + 1) / (2 * nodeCount)) * (end - start) / 2;
    }
    return nodes;
}
double* divDiff(double* xNodes, double* yNodes, int nodeCount) {
    double* divDif = new double[nodeCount];
    for (int i = 0; i < nodeCount; i++) {
        divDif[i] = yNodes[i];
    }

    for (int i = 1; i < nodeCount; i++) {
        for (int j = 0; j < nodeCount - i; j++) {
            divDif[j] = (divDif[j + 1] - divDif[j]) / (xNodes[i + j] - xNodes[j]);
        }
    }
    return divDif;
}
double newtonInterpBack(double x, double* xNodes, double* yNodes, int nodeCount, double* divDif) {

    double result = 0, t;

    for (int i = nodeCount - 1; i >= 0; i--) {
        t = 1;
        for (int j = nodeCount - 1; j > i; j--) {
            t *= (x - xNodes[j]);
        }
        result += divDif[i] * t;
    }


    return result;
}

void writeNodesToFile(double* x, double (*func)(double), int nodeCount, const string& filename) {
    ofstream file(filename);

    for (int i = 0; i < nodeCount; ++i) {
        file << x[i] << " " << func(x[i]) << endl;
    }

    file.close();
}

void writePolynomialsToFile(double* x, double* y, int nodeCount, const string& filename) {
    ofstream file(filename);
    double* divDif = divDiff(x, y, nodeCount);
    const int pointsCount = 1000;
    double step = (x[nodeCount - 1] - x[0]) / (pointsCount - 1);
    for (int i = 0; i < pointsCount; ++i) {
        double xi = x[0] + i * step;
        file << xi << " " << newtonInterpBack(xi, x, y, nodeCount,divDif) << endl;
    }

    file.close();
}

void WriteMaxErrors(double (*func)(double), const string& filename) {
    ofstream file(filename);

    const int maxNodes = 100;
    for (int nodeCount = 5; nodeCount <= maxNodes; ++nodeCount) {
        double* xNodes = genChebNodes(-5.0, 5.0, nodeCount);
        double* yNodes = new double[nodeCount];
        for (int i = 0; i < nodeCount; ++i) {
            yNodes[i] = func(xNodes[i]);
        }
        double* divDif = divDiff(xNodes, yNodes, nodeCount);
        double maxError = 0;
        for (int i = 0; i < nodeCount; ++i) {
            double interpolatedValue = newtonInterpBack(xNodes[i], xNodes, yNodes, nodeCount, divDif);
            double exactValue = func(xNodes[i]);
            double error = fabs(exactValue - interpolatedValue);
            if (error > maxError) {
                maxError = error;
            }
        }
        file << nodeCount << " " << maxError << endl;
        delete[] xNodes;
        delete[] yNodes;
    }

    file.close();
}
void writeSelectErrors(double (*func)(double), const string& filename) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file " << filename << endl;
        exit(EXIT_FAILURE);
    }

    const int maxNodes = 100;
    for (int nodeCount = 5; nodeCount <= maxNodes; nodeCount += 5) {

       /* double* divDif2_uni = divDiff(xk2_uniform, yk2_uniform, node);
        double* divDif1_cheb = divDiff(xk1_chebyshev, yk1_chebyshev, node);
        double* divDif2_cheb = divDiff(xk2_chebyshev, yk2_chebyshev, node);*/
        double* xNodes = genUniformNodes(-5.0, 5.0, nodeCount);
        double* yNodes = new double[nodeCount];
        double* divDif1_uni = divDiff(xNodes, yNodes, maxNodes);
        for (int i = 0; i < nodeCount; ++i) {
            yNodes[i] = func(xNodes[i]);
        }


       // double error1 = fabs(newtonInterpBack(1.0, xNodes, yNodes, nodeCount) - func(1.0));
        double error3 = fabs(newtonInterpBack(3.0, xNodes, yNodes, nodeCount,divDif1_uni) - func(3.0));


        file << nodeCount << " " << error3 /*<< " " << error3*/ << endl;

        delete[] xNodes;
        delete[] yNodes;
    }

    file.close();
}

void writeTheoryErrors(double (*func)(double), int nodeCount, const string& filename) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file " << filename << endl;
        exit(EXIT_FAILURE);
    }

    const int pointsCount = 1000;
    double step = 10.0 / (pointsCount - 1);

    double* xNodes = genUniformNodes(-5.0, 5.0, nodeCount);

    double factorial = 1.0;
    for (int i = 2; i <= nodeCount; ++i) {
        factorial *= i;
    }


    int n = 0;
    for (int i = 0; i < pointsCount; ++i) {
        double xi = -5.0 + i * step;
        double product = 1.0;
        for (int i = 0; i < nodeCount; ++i) {
         /*   if (xNodes[i] != -5.0 + i * step) {*/
                product *= (xi - xNodes[i]);

      /*      }
            else ++i;*/
        }

        double theoreticalError = fabs(1 * product / factorial);

        file << xi << " " << theoreticalError << endl;

    }

    delete[] xNodes;
    file.close();
}

int main() {
    double a = -5.0;
    double b = 5.0;
    int node = 10;

    double* xk1_uniform = genUniformNodes(a, b, node);
    double* xk1_chebyshev = genChebNodes(a, b, node);
    double* xk2_uniform = genUniformNodes(a, b, node);
    double* xk2_chebyshev = genChebNodes(a, b, node);

    writeNodesToFile(xk1_uniform, f1, node, "nodes_uniform_func1.txt");
    writeNodesToFile(xk1_chebyshev, f1, node, "nodes_chebyshev_func1.txt");
    writeNodesToFile(xk2_uniform, f2, node, "nodes_uniform_func2.txt");
    writeNodesToFile(xk2_chebyshev, f2, node, "nodes_chebyshev_func2.txt");
    double* x1 = genUniformNodes(a, b, node);
    double* x2 = genChebNodes(a, b, node);
    // Рассчитываем значения полиномов Ньютона и записываем их в файлы
    double* yk1_uniform = new double[node];
    double* yk1_chebyshev = new double[node];
    double* yk2_uniform = new double[node];
    double* yk2_chebyshev = new double[node];
    double* y1 = new double[node];
    for (int i = 0; i < node; ++i) {
        yk1_uniform[i] = f1(xk1_uniform[i]);
        yk1_chebyshev[i] = f1(xk1_chebyshev[i]);
        yk2_uniform[i] = f2(xk2_uniform[i]);
        yk2_chebyshev[i] = f2(xk2_chebyshev[i]);



    }

    writePolynomialsToFile(xk1_uniform, yk1_uniform, node, "P1uni.txt");
    writePolynomialsToFile(xk1_chebyshev, yk1_chebyshev, node, "P1cheb.txt");
    writePolynomialsToFile(xk2_uniform, yk2_uniform, node, "P2uni.txt");
    writePolynomialsToFile(xk2_chebyshev, yk2_chebyshev, node, "P2cheb.txt");



    ofstream error_file1_uniform("errors_func1_uniform.txt");
    ofstream error_file1_chebyshev("errors_func1_chebyshev.txt");
   ofstream error_file2_uniform("errors_func2_uniform.txt");
   ofstream error_file2_chebyshev("errors_func2_chebyshev.txt");

    if (!error_file1_uniform.is_open() || !error_file1_chebyshev.is_open() || !error_file2_uniform.is_open() || !error_file2_chebyshev.is_open()) {
       cerr << "Error opening file(s) for writing errors" << endl;
       exit(EXIT_FAILURE);
    }

    double* divDif1_uni = divDiff(xk1_uniform, yk1_uniform, node);
    double* divDif2_uni = divDiff(xk2_uniform, yk2_uniform, node);
    double* divDif1_cheb = divDiff(xk1_chebyshev, yk1_chebyshev, node);
    double* divDif2_cheb = divDiff(xk2_chebyshev, yk2_chebyshev, node);



    const int points_count = 1000;
    double step1 = (xk1_uniform[node - 1] - xk1_uniform[0]) / (points_count - 1);
    double step2 = (xk2_uniform[node - 1] - xk2_uniform[0]) / (points_count - 1);
    for (int i = 0; i < points_count; ++i) {
        double xi1 = xk1_uniform[0] + i * step1;
        double xi2 = xk2_uniform[0] + i * step2;

        double y_exact1_uniform = f1(xi1);
        double y_exact1_chebyshev = f1(xi1);
        double y_exact2_uniform = f2(xi2);
        double y_exact2_chebyshev = f2(xi2);

        double y_interpolated1_uniform = newtonInterpBack(xi1, xk1_uniform, yk1_uniform, node, divDif1_uni);
        double y_interpolated1_chebyshev = newtonInterpBack(xi1, xk1_chebyshev, yk1_chebyshev, node, divDif1_cheb);
        double y_interpolated2_uniform = newtonInterpBack(xi2, xk2_uniform, yk2_uniform, node, divDif2_uni);
        double y_interpolated2_chebyshev = newtonInterpBack(xi2, xk2_chebyshev, yk2_chebyshev, node, divDif2_cheb);

        double error1_uniform = fabs(y_exact1_uniform - y_interpolated1_uniform);
        double error1_chebyshev = fabs(y_exact1_chebyshev - y_interpolated1_chebyshev);
        double error2_uniform = fabs(y_exact2_uniform - y_interpolated2_uniform);
        double error2_chebyshev = fabs(y_exact2_chebyshev - y_interpolated2_chebyshev);
        error_file1_uniform << xi1 << " " << error1_uniform << endl;
        error_file1_chebyshev << xi1 << " " << error1_chebyshev << endl;
        error_file2_uniform << xi2 << " " << error2_uniform << endl;
        error_file2_chebyshev << xi2 << " " << error2_chebyshev << endl;
    }

    error_file1_uniform.close();
    error_file1_chebyshev.close();
    error_file2_uniform.close();
    error_file2_chebyshev.close();
   //writeMaxErrors(f1, "max_errors_uniform_func1.txt");
   //writeMaxErrors(f1, "max_errors_chebyshev_func1.txt");
   //writeMaxErrors(f2, "max_errors_uniform_func2.txt");
   //writeMaxErrors(f2, "max_errors_chebyshev_func2.txt");
    writeSelectErrors(f1, "select_uniform_func1.txt");
    writeSelectErrors(f2, "select_uniform_func2.txt");
    //writeSelectErrors(f1, "select_chebyshev_func1.txt");
    //writeSelectErrors(f2, "select_chebyshev_func2.txt");
    writeTheoryErrors(f1, node, "theory_func1_uniform.txt");
   // writeTheoryErrors(f1, node, "theory_func1_chebyshev.txt");

    delete[] xk1_uniform;
    delete[] xk1_chebyshev;
    delete[] xk2_uniform;
    delete[] xk2_chebyshev;
    delete[] yk1_uniform;
    delete[] yk1_chebyshev;
    delete[] yk2_uniform;
    delete[] yk2_chebyshev;

    return 0;

}
