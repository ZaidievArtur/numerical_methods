data = load('data.txt');

% Extracting data for plotting
x_exact = linspace(1, 3, 100);
y_exact = exp(x_exact) .* (log(abs(x_exact)) + 1);

x_h1 = data(1:5, 1);
y_h1 = data(1:5, 2);

x_h2 = data(6:14, 1);
y_h2 = data(6:14, 2);

% Load error data
errors = data(15:end, 1);
epsilons = data(15:end, 2);
step_sizes = data(15:end, 3);

% Extracting data for the first set of errors (n=4)
x_errors_h1 = errors(1:4);
errors_h1 = errors(1:4);

% Extracting data for the second set of errors (n=8)
x_errors_h2 = errors(5:9);
errors_h2 = errors(5:9);

% Plotting
subplot(2, 1, 1);
plot(x_exact, y_exact, 'r-', 'LineWidth', 2);
hold on;
plot(x_h1, y_h1, 'g--', 'LineWidth', 2);
plot(x_h2, y_h2, 'b-.', 'LineWidth', 2);
xlabel('x');
ylabel('y');
legend('Exact solution', 'Kutta-Merson for n = 4', 'Kutta-Merson for n = 8');
grid on;

subplot(2, 1, 2);
loglog(epsilons, errors, 'o-', 'LineWidth', 2);
xlabel('/epsilon');
ylabel('error');
grid on;

figure;
loglog(x_errors_h1, errors_h1, '-g', 'LineWidth', 2);
hold on;
loglog(x_errors_h2, errors_h2, '-b', 'LineWidth', 2);
xlabel('Step Size');
ylabel('Error');
legend('Error for n = 4', 'Error for n = 8');
grid on;

