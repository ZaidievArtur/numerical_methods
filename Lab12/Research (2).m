% Загрузка данных из файла
filename = 'output.txt';

data = dlmread(filename, ' ');
% Разделение данных на отдельные столбцы
x = data(:,1);
y1 = data(:, 2);
y2 = data(:, 3);
spline1 = data(:, 4);
spline2 = data(:, 5);

% Построение графиков
figure;

% График функции y1 и ее сплайна
subplot(2, 1, 1);
semilogx(x, y1, 'b-', x, spline1, 'r--', 'LineWidth', 1.5);
legend('Функция y1', 'Сплайн y1');
xlabel('x');
ylabel('y');
title('График функции y1 и ее сплайна');

% График функции y2 и ее сплайна
subplot(2, 1, 2);
semilogx(x, y2, 'b-', x, spline2, 'r--', 'LineWidth', 1.5);
legend('Функция y2', 'Сплайн y2');
xlabel('x');
ylabel('y');
title('График функции y2 и ее сплайна');

