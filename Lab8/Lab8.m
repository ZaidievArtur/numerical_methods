% Загрузка данных из файлов
data_uniform_func1 = load('errors_func1_uniform.txt');
data_chebyshev_func1 = load('errors_func1_chebyshev.txt');
data_uniform_func2 = load('errors_func2_uniform.txt');
data_chebyshev_func2 = load('errors_func2_chebyshev.txt');

% Загрузка значений узлов из файлов
nodes_uniform_func1 = load('nodes_uniform_func1.txt');
nodes_chebyshev_func1 = load('nodes_chebyshev_func1.txt');
nodes_uniform_func2 = load('nodes_uniform_func2.txt');
nodes_chebyshev_func2 = load('nodes_chebyshev_func2.txt');

% Создание нового окна для всех графиков
figure;

% График поточечной ошибки для функции 1 на равномерной сетке
subplot(2, 2, 1);
plot(data_uniform_func1(:, 1), data_uniform_func1(:, 2), 'b-', 'LineWidth', 0.5); % График ошибки для равномерной сетки
hold on;
plot(data_uniform_func1(:, 1), zeros(size(data_uniform_func1(:, 1))), 'g-'); % Горизонтальная линия на уровне нуля
plot(nodes_uniform_func1(:, 1), zeros(size(nodes_uniform_func1(:, 1))), 'ro', 'MarkerSize', 5); % Узлы равномерной сетки
%title('Рис. 7','FontSize', 20);
xlabel('x');
ylabel('|Rn(x)|');
legend('Error x - sin(x) - 0.25','Teoretic error line', 'Nodes', 'Location', 'north');
grid on;
hold off;

% График поточечной ошибки для функции 1 на сетке Чебышёва
subplot(2, 2, 2);
plot(data_chebyshev_func1(:, 1), data_chebyshev_func1(:, 2), 'b-', 'LineWidth', 0.5); % График ошибки для сетки Чебышёва
hold on;
plot(data_chebyshev_func1(:, 1), zeros(size(data_chebyshev_func1(:, 1))), 'g-'); % Горизонтальная линия на уровне нуля
plot(nodes_chebyshev_func1(:, 1), zeros(size(nodes_chebyshev_func1(:, 1))), 'ro', 'MarkerSize', 5); % Узлы сетки Чебышёва
%title('Рис. 8','FontSize', 20);
xlabel('x');
ylabel('|Rn(x)|');
legend('Error x - sin(x) - 0.25', 'Teoretic error line', 'Nodes', 'Location', 'northwest');
grid on;
hold off;

% График поточечной ошибки для функции 2 на равномерной сетке
subplot(2, 2, 3);
plot(data_uniform_func2(:, 1), data_uniform_func2(:, 2), 'b-', 'LineWidth', 0.5); % График ошибки для равномерной сетки
hold on;
plot(data_uniform_func2(:, 1), zeros(size(data_uniform_func2(:, 1))), 'g-'); % Горизонтальная линия на уровне нуля
plot(nodes_uniform_func2(:, 1), zeros(size(nodes_uniform_func2(:, 1))), 'ro', 'MarkerSize', 5); % Узлы равномерной сетки
xlabel('x');
ylabel('|Rn(x)|');
legend('Error x^5 + 0.4*sign(x) * x^4 + 2', 'Teoretic error line','Nodes', 'Location', 'north');
grid on;
hold off;

% График поточечной ошибки для функции 2 на сетке Чебышёва
subplot(2, 2, 4);
plot(data_chebyshev_func2(:, 1), data_chebyshev_func2(:, 2), 'b-', 'LineWidth', 0.5); % График ошибки для сетки Чебышёва
hold on;
plot(data_chebyshev_func2(:, 1), zeros(size(data_chebyshev_func2(:, 1))), 'g-'); % Горизонтальная линия на уровне нуля
plot(nodes_chebyshev_func2(:, 1), zeros(size(nodes_chebyshev_func2(:, 1))), 'ro', 'MarkerSize', 5); % Узлы сетки Чебышёва
xlabel('x');
ylabel('|Rn(x)|');
legend('Error x^5 + 0.4*sign(x) * x^4 + 2', 'Teoretic error line','Nodes','Location', 'northwest');
grid on;
hold off;

