% Load data
solution5 = load("nodes5.txt");
solution10 = load("nodes10.txt");
solution_refined = load("nodes_refined.txt");
error_data = load("error.txt");
delta_data = load("delta.txt");

% Graph for exact and numerical solutions for fixed step sizes
figure;
plot(solution5(:, 1), solution5(:, 2), 'b-', solution5(:, 1), solution5(:, 3), 'r-', solution10(:, 1), solution10(:, 2), 'g-');
grid on;
legend('n=5', 'accurate', 'n=10','Location', 'northwest');
xlabel('x');
ylabel('y');


% Separate graph for refined grid solution
figure;
plot(solution_refined(:, 1), solution_refined(:, 2), 'g-', solution_refined(:, 1), solution_refined(:, 3), 'r-');
grid on;
legend('condensed', 'accurate');
xlabel('x');
ylabel('y');


% Graph for error with fixed step sizes
figure;
plot(solution5(:, 1), solution5(:, 4), solution10(:, 1), solution10(:, 4));
grid on;
legend('n=5', 'n=10','Location', 'northwest');
xlabel('x');
ylabel('Error');


% Separate graph for error with refined grid
figure;
plot(solution_refined(:, 1), solution_refined(:, 4));
grid on;
xlabel('x');
ylabel('Error');


% Graph for dependence of actual error on given accuracy
figure;
loglog(error_data(:, 1), error_data(:, 2));
grid on;
xlabel('Epsilon');
ylabel('Error');



