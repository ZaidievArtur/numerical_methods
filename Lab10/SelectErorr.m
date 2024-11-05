% Загрузка данных из файлов
data_uniform_func1 = load('select_uniform_func1.txt');
data_chebyshev_func1 = load('select_chebyshev_func1.txt');
data_uniform_func2 = load('select_uniform_func2.txt');
data_chebyshev_func2 = load('select_chebyshev_func2.txt');

% Извлечение данных
nodes_uniform_func1 = data_uniform_func1(:, 1);
errors_uniform_func1 = data_uniform_func1(:, 2);

nodes_chebyshev_func1 = data_chebyshev_func1(:, 1);
errors_chebyshev_func1 = data_chebyshev_func1(:, 2);

nodes_uniform_func2 = data_uniform_func2(:, 1);
errors_uniform_func2 = data_uniform_func2(:, 2);

nodes_chebyshev_func2 = data_chebyshev_func2(:, 1);
errors_chebyshev_func2 = data_chebyshev_func2(:, 2);

% Построение графиков
figure;
hold on;

semilogy(nodes_uniform_func1, errors_uniform_func1, 'b-');
semilogy(nodes_chebyshev_func1, errors_chebyshev_func1, 'r-');
semilogy(nodes_uniform_func2, errors_uniform_func2, 'g-');
semilogy(nodes_chebyshev_func2, errors_chebyshev_func2, 'm-');

hold off;

xlabel('Nodes','FontSize', 15);
ylabel('(|Rn(3)|)','FontSize', 15 );
legend('x - sin(x) - 0.25 Uniform','x - sin(x) - 0.25 Chebyshev','x^5 + 0.4 * sign(x) * x^4 + 2 Uniform','x^5 + 0.4 * sign(x) * x^4 + 2 Chebyshev', 'Location', 'north');
grid on;

