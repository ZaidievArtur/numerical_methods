% Загрузка данных из файлов
data_uniform_func1 = load('nodes_uniform_func1.txt');
data_chebyshev_func1 = load('nodes_chebyshev_func1.txt');
data_uniform_func2 = load('nodes_uniform_func2.txt');
data_chebyshev_func2 = load('nodes_chebyshev_func2.txt');
data_P1uni = load('P1uni.txt');
data_P1cheb = load('P1cheb.txt');
data_P2uni = load('P2uni.txt');
data_P2cheb = load('P2cheb.txt');

% Извлечение узлов и значений интерполяционных полиномов
x_uniform_func1 = data_uniform_func1(:, 1);
y_uniform_func1 = data_uniform_func1(:, 2);
x_chebyshev_func1 = data_chebyshev_func1(:, 1);
y_chebyshev_func1 = data_chebyshev_func1(:, 2);
x_uniform_func2 = data_uniform_func2(:, 1);
y_uniform_func2 = data_uniform_func2(:, 2);
x_chebyshev_func2 = data_chebyshev_func2(:, 1);
y_chebyshev_func2 = data_chebyshev_func2(:, 2);
x_P1uni = data_P1uni(:, 1);
y_P1uni = data_P1uni(:, 2);
x_P1cheb = data_P1cheb(:, 1);
y_P1cheb = data_P1cheb(:, 2);
x_P2uni = data_P2uni(:, 1);
y_P2uni = data_P2uni(:, 2);
x_P2cheb = data_P2cheb(:, 1);
y_P2cheb = data_P2cheb(:, 2);

% Определение функций Func1 и Func2 через @(x)
Func1 = @(x) x - sin(x) - 0.25;
Func2 = @(x) x .^ 5 + 0.4 * sign(x) .* x .^ 4 + 2;

% Построение проверочной сетки
x_check = linspace(min(x_uniform_func1), max(x_uniform_func1), 10000);

% Вычисление значений функций на проверочной сетке
y_check_func1 = Func1(x_check);
y_check_func2 = Func2(x_check);

% Построение графиков интерполяционных полиномов
figure;
% График для функции 1 на равномерной сетке
subplot(2, 2, 1);
plot(x_uniform_func1, y_uniform_func1, 'ko', x_P1uni, y_P1uni, 'LineWidth', 1.3, 'r-', x_check, y_check_func1, 'g-');
legend('Nodes', 'Interpolation', 'x - sin(x) - 0.25','Location', 'southeast');

xlabel('x');
ylabel('y');
% График для функции 1 на сетке Чебышёва
subplot(2, 2, 2);
plot(x_chebyshev_func1, y_chebyshev_func1, 'ko', x_P1cheb, y_P1cheb, 'LineWidth', 1.3, 'r-', x_check, y_check_func1, 'g-');
legend('Nodes', 'Interpolation', 'x - sin(x) - 0.25', 'Location', 'southeast');
xlabel('x');
ylabel('y');
% График для функции 2 на равномерной сетке
subplot(2, 2, 3);
plot(x_uniform_func2, y_uniform_func2, 'ko', x_P2uni, y_P2uni, 'LineWidth', 1.3, 'r-', x_check, y_check_func2, 'g-');
legend('Nodes', 'Interpolation', 'x^5 + 0.4 * sign(x) * x^4 + 2','Location', 'southeast');
xlabel('x');
ylabel('y');
% График для функции 2 на сетке Чебышёва
subplot(2, 2, 4);
plot(x_chebyshev_func2, y_chebyshev_func2, 'ko', x_P2cheb, y_P2cheb, 'LineWidth', 1.3, 'r-', x_check, y_check_func2, 'g-');
legend('Nodes', 'Interpolation', 'x^5 + 0.4 * sign(x) * x^4 + 2','Location', 'southeast');
xlabel('x');
ylabel('y');

