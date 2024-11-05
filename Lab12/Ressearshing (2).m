% Чтение данных из файла data.txt
data = load('data.txt');

% Извлечение данных для построения графиков
error = data(:, 1);
iterations = data(:, 2);
log_n = data(:, 3);

% Логарифм по основанию 2 фактической ошибки
log_error = log2(error);

% Построение графика зависимости log_2(фактической ошибки) от log_2(длины отрезка разбиения)
figure;
plot(log_n, log_error, '-o', 'LineWidth', 1.5);
xlabel('log_2(Длина отрезка разбиения)');
ylabel('log_2(Фактическая ошибка)');
title('Зависимость log_2(фактической ошибки) от log_2(длины отрезка разбиения)');
grid on;

% Линейная регрессия
coeffs = polyfit(log_n, log_error, 1);
k = coeffs(1); % Порядок точности
log_C = coeffs(2); % Логарифм константы

% Добавление линии регрессии на график
hold on;
regression_line = polyval(coeffs, log_n);
plot(log_n, regression_line, '--', 'LineWidth', 1.5);
legend('Данные', 'Линейная регрессия');

% Вывод результатов
fprintf('Порядок точности (k): %f\n', k);
fprintf('Константа (C): %f\n', 2^log_C);

