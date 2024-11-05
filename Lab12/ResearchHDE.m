% Load data from the file
data = dlmread('data.txt');

% Extract different sections of the data
error_data = data(1:9, :);
h1_data = data(10:14, :);
h2_data = data(15:23, :);

% Extract variables for exact solution
a = 1;
b = 3;
x_exact = linspace(a, b, 100);
y_exact = exp(x_exact) .* (log(abs(x_exact)) + 1);

% Plot the exact solution
figure;
plot(x_exact, y_exact, 'r-', 'LineWidth', 2);
hold on;

% Plot the numerical solutions for n = 4 (h = 0.375)
plot(h1_data(:, 1), h1_data(:, 2), 'g--', 'LineWidth', 2);

% Plot the numerical solutions for n = 8 (h = 0.1875)
plot(h2_data(:, 1), h2_data(:, 2), 'b-.', 'LineWidth', 2);

% Set plot properties
xlabel('x');
ylabel('y');
legend('Точное решение', 'Кутты-Мерсона для n = 4', 'Кутты-Мерсона для n = 8');
title('Сравнение точного и численных решений');
grid on;

% Save the figure
saveas(gcf, 'comparison_plot.png');

% Plot the error data
figure;
semilogy(error_data(:, 2), error_data(:, 1), 'o-');
xlabel('epsilon');
ylabel('error');
title('Ошибка метода Кутты-Мерсона');
grid on;

% Save the figure
saveas(gcf, 'error_plot.png');

