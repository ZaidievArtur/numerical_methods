% Загрузка данных
data_uniform_func1 = load("max_errors_uniform_func1.txt");
data_chebyshev_func1 = load("max_errors_chebyshev_func1.txt");
data_uniform_func2 = load("max_errors_uniform_func2.txt");
data_chebyshev_func2 = load("max_errors_chebyshev_func2.txt");

% Извлечение данных
nodes_uniform_func1 = data_uniform_func1(:, 1);
errors_uniform_func1 = data_uniform_func1(:, 2);

nodes_chebyshev_func1 = data_chebyshev_func1(:, 1);
errors_chebyshev_func1 = data_chebyshev_func1(:, 2);

nodes_uniform_func2 = data_uniform_func2(:, 1);
errors_uniform_func2 = data_uniform_func2(:, 2);

nodes_chebyshev_func2 = data_chebyshev_func2(:, 1);
errors_chebyshev_func2 = data_chebyshev_func2(:, 2);

% Построение графика
plot(nodes_uniform_func1, errors_uniform_func1, "-o", 'DisplayName', "Uniform Nodes - Function 1");
hold on;
plot(nodes_chebyshev_func1, errors_chebyshev_func1, "-o", 'DisplayName', "Chebyshev Nodes - Function 1");
plot(nodes_uniform_func2, errors_uniform_func2, "-o", 'DisplayName', "Uniform Nodes - Function 2");
plot(nodes_chebyshev_func2, errors_chebyshev_func2, "-o", 'DisplayName', "Chebyshev Nodes - Function 2");

xlabel("Number of Nodes");
ylabel("Max Interpolation Error");
title("Dependence of Interpolation Error on Number of Nodes");
legend;
grid on;

