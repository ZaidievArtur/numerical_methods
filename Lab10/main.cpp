#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

// Функция для вычисления значения функции y = x - sin(x) - 0.25
double function1(double x) {
    return x - sin(x) - 0.25;
}

// Функция для вычисления значения функции y = x^5 + 0.4 * sign(x) * x^4 + 2
double function2(double x) {
    return pow(x, 5) + 0.4 * copysign(1.0, x) * pow(x, 4) + 2;
}

// Функция для генерации точек на отрезке
vector<double> generatePoints(double start, double end, int count) {
    vector<double> points(count);
    double step = (end - start) / (count - 1);
    for (int i = 0; i < count; ++i) {
        points[i] = start + i * step;
    }
    return points;
}

// Структура для хранения коэффициентов кубического сплайна на каждом сегменте
struct CubicSplineSegment {
    double a, b, c, d, x;
};

// Функция для вычисления кубического сплайна
vector<CubicSplineSegment> computeCubicSpline(const vector<double>& x, const vector<double>& y) {
    int n = x.size() - 1;
    vector<double> h(n), alpha(n), l(n + 1), mu(n), z(n + 1);

    for (int i = 0; i < n; ++i) {
        h[i] = x[i + 1] - x[i];
    }

    for (int i = 1; i < n; ++i) {
        alpha[i] = (3.0 / h[i]) * (y[i + 1] - y[i]) - (3.0 / h[i - 1]) * (y[i] - y[i - 1]);
    }

    l[0] = 1;
    mu[0] = 0;
    z[0] = 0;
    for (int i = 1; i < n; ++i) {
        l[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    vector<double> c(n + 1), b(n), d(n);
    l[n] = 1;
    z[n] = 0;
    c[n] = 0;
    for (int j = n - 1; j >= 0; --j) {
        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
        d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
    }

    vector<double> a(n), M(n);
    for (int i = 0; i < n; ++i) {
        a[i] = y[i];
        M[i] = c[i];
    }

    vector<CubicSplineSegment> segments;
    for (int i = 0; i < n; ++i) {
        CubicSplineSegment segment;
        segment.a = a[i];
        segment.b = b[i];
        segment.c = M[i];
        segment.d = d[i];
        segment.x = x[i];
        segments.push_back(segment);
    }

    return segments;
}

// Функция для вычисления значения кубического сплайна в точке
double evaluateCubicSpline(const vector<CubicSplineSegment>& segments, double x) {
    int n = segments.size();
    int i = 0;
    while (i < n && x > segments[i].x) {
        ++i;
    }
    if (i == n) {
        --i;
    }
    double dx = x - segments[i].x;
    return segments[i].a + segments[i].b * dx + segments[i].c * dx * dx + segments[i].d * dx * dx * dx;
}

int main() {
    // Задание параметров отрезка и количества точек
    double start = -5.0;
    double end = 5.0;
    int pointCount = 1000;

    // Генерация точек на отрезке
    vector<double> points = generatePoints(start, end, pointCount);

    // Вычисление значений функций на точках
    vector<double> y1(pointCount), y2(pointCount);
    for (int i = 0; i < pointCount; ++i) {
        y1[i] = function1(points[i]);
        y2[i] = function2(points[i]);
    }

    // Вычисление кубических сплайнов
    vector<CubicSplineSegment> spline1 = computeCubicSpline(points, y1);
    vector<CubicSplineSegment> spline2 = computeCubicSpline(points, y2);

    // Запись значений в файл
    ofstream outputFile("output.txt");
    for (int i = 0; i < pointCount; ++i) {
        double splineValue1 = evaluateCubicSpline(spline1, points[i]);
        double splineValue2 = evaluateCubicSpline(spline2, points[i]);
        outputFile << points[i] << " " << y1[i] << " " << y2[i] << " " << splineValue1 << " " << splineValue2 << endl;
    }
    outputFile.close();

    return 0;
}
