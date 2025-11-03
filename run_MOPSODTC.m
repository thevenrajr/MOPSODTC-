% =========================================================================
% Main script to run the MOPSODTC algorithm
% =========================================================================
clc;
clear;
close all;

%% Problem Definition
problem_name = 'ZDT3'; % Select the problem: 'ZDT1', 'ZDT2', 'ZDT3', 'ZDT4', 'ZDT6'
problem = TestProblems(problem_name);
CostFunction = problem.CostFunction;
nVar = problem.nVar;
VarMin = problem.VarMin;
VarMax = problem.VarMax;

%% MOPSODTC Parameters
MaxIt = 200;      % Maximum Number of Iterations (adjust based on function evaluations, e.g., 10000 evals / 200 pop_size = 50 iterations)
N = 200;          % Total Population Size
archive_size = 200; % Archive Size

gamma = 50;       % Initial size control parameter (from paper, Section 4.3)
c = 1.0;          % Perturbation step size for ES (from paper, Section 3.2)
w = 0.4;          % Inertia Weight for DS
c1 = 2.0;         % Personal Learning Coefficient for DS
c2 = 2.0;         % Global Learning Coefficient for DS

% Perturbation Parameters
pCrossover = 0.9;   % Crossover Probability
pMutation = 1/nVar; % Mutation Probability
eta_c = 20;         % Distribution Index for Crossover
eta_m = 20;         % Distribution Index for Mutation

%% Run MOPSODTC Algorithm
fprintf('Running MOPSODTC on %s...\n', problem_name);
[archive, execution_time] = MOPSODTC(problem, MaxIt, N, archive_size, gamma, c, w, c1, c2, pCrossover, pMutation, eta_c, eta_m);
fprintf('Execution finished in %.2f seconds.\n', execution_time);

%% Results
% Extract final solution set
final_costs = vertcat(archive.Cost);

% Plotting the results
figure;
hold on;
plot(final_costs(:,1), final_costs(:,2), 'r*', 'MarkerSize', 8);

% Plot true Pareto Front if available
if isfield(problem, 'TruePF')
    plot(problem.TruePF(:,1), problem.TruePF(:,2), 'k-', 'LineWidth', 2);
    legend('MOPSODTC Solutions', 'True Pareto Front');
else
    legend('MOPSODTC Solutions');
end

title(['MOPSODTC on ' problem_name]);
xlabel('f_1');
ylabel('f_2');
grid on;
box on;
