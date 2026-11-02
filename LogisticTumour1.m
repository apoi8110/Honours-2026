%% PURE LOGISTIC TUMOUR MODEL

% Model: dC/dt = r_C*C*(1 - C/K_C)
% with r_C and K_C fitted to untreated tumour-volume data.

% Figure 1 - untreated data + logistic fit


clear;
clc;
close all;

col = honoursPlotStyle();

% Loading untreated tumour data

data = readtable('Data.csv');

t_data = [4; 7; 11; 14; 18; 21; 25];
C_data = data{:,2};

t_data = t_data(:);
C_data = C_data(:);

%Initial condition

C0 = 10 + 1e-6;

% Initial guess to sweep around 

r_guess = 0.25;       % per day
K_guess = 600;        % mm^3

% Log-transform to make sure it is positive
theta0 = log([r_guess K_guess]);

% Fit r and K

options = optimset('Display', 'off', ...
    'MaxFunEvals', 5000, ...
    'MaxIter', 5000);

theta_fit = fminsearch( ...
    @(theta) logisticObjective(theta, ...
    t_data, C_data, C0), ...
    theta0, options);

params = exp(theta_fit);

r_C = params(1);
K_C = params(2);

% Smooth fitted trajectory

t_plot = linspace(t_data(1), t_data(end), 500);

[t_model, C_model] = ode45( ...
    @(t,C) logisticRHS(t, C, r_C, K_C), ...
    t_plot, C0);

% Predictions at same times as experiment

[~, C_at_data] = ode45( ...
    @(t,C) logisticRHS(t, C, r_C, K_C), ...
    t_data, C0);

% Residuals

residuals = C_at_data - C_data;

SSE  = sum(residuals.^2);
RMSE = sqrt(mean(residuals.^2));

% Print results

fprintf('\nPURE LOGISTIC MODEL\n');
fprintf('r_C  = %.6f /day\n', r_C);
fprintf('K_C  = %.6f mm^3\n', K_C);
fprintf('SSE  = %.6f\n', SSE);
fprintf('RMSE = %.6f mm^3\n\n', RMSE);


%% FIGURE 1: Logistic fit

figure('Position', [100 100 760 480]);

plot(t_model, C_model, ...
    'Color', col.total, ...
    'LineWidth', 3);

hold on;

plot(t_data, C_data, 'o', ...
    'Color', col.data, ...
    'MarkerFaceColor', 'w', ...
    'MarkerSize', 8, ...
    'LineWidth', 1.8);

xlabel('Time (days)');
ylabel('Tumour volume $C(t)$ (mm$^3$)');

title('Logistic tumour growth model');

legend( ...
    'Logistic model', ...
    'Experimental data', ...
    'Location', 'northwest');

xlim([4 25]);

%% FUNCTION

function error = logisticObjective(theta, ...
    t_data, C_data, C0)

    parameters = exp(theta);

    r = parameters(1);
    K = parameters(2);

    [~, C] = ode45( ...
        @(t,C) logisticRHS(t, C, r, K), ...
        t_data, C0);

    residuals = C - C_data;

    error = sum(residuals.^2);

end


function dCdt = logisticRHS(~, C, r, K)

    dCdt = r*C*(1 - C/K);

end