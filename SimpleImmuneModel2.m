%% SIMPLE SINGLE-POPULATION TUMOUR-IMMUNE MODEL
%  dC/dt = r_C*C*(1 - C/K_C) - k_C*T*C
%  dT/dt = s - d*T + a*C*T

% r_C and K_C are fitted.
% Immune parameters are fixed. 

% Figures:
%   Figure 1 - tumour fit + T-cell trajectory
%   Figure 2 - residuals


clear;
clc;
close all;

col = honoursPlotStyle();

% Load data

data = readtable('Data.csv');

t_data = [4; 7; 11; 14; 18; 21; 25];
C_data = data{:,2};

t_data = t_data(:);
C_data = C_data(:);

% Immune parameters

k_C = 1e-5;

s = 10;
d = 0.2;
a = 1e-4;

% Initial conditions

C0 = 10 + 1e-6;

T0 = s/d;

y0 = [C0; T0];

% Fit tumour parameters

r_guess = 0.25;
K_guess = 600;

theta0 = log([r_guess K_guess]);

options = optimset( ...
    'Display', 'off', ...
    'MaxFunEvals', 5000, ...
    'MaxIter', 5000);

theta_fit = fminsearch( ...
    @(theta) tumourImmuneObjective(theta, ...
    t_data, C_data, y0, ...
    k_C, s, d, a), ...
    theta0, options);

parameters = exp(theta_fit);

r_C = parameters(1);
K_C = parameters(2);

% Smooth solution

t_plot = linspace(4, 25, 500);

[t_model, y_model] = ode45( ...
    @(t,y) tumourImmuneRHS(t, y, ...
    r_C, K_C, k_C, s, d, a), ...
    t_plot, y0);

C_model = y_model(:,1);
T_model = y_model(:,2);

% Predictions at data points

[~, y_data_model] = ode45( ...
    @(t,y) tumourImmuneRHS(t, y, ...
    r_C, K_C, k_C, s, d, a), ...
    t_data, y0);

C_at_data = y_data_model(:,1);

% Fit quality

residuals = C_at_data - C_data;

SSE  = sum(residuals.^2);
RMSE = sqrt(mean(residuals.^2));

fprintf('\nSIMPLE TUMOUR-IMMUNE MODEL\n');
fprintf('--------------------------------\n');
fprintf('r_C  = %.6f /day\n', r_C);
fprintf('K_C  = %.6f mm^3\n', K_C);
fprintf('SSE  = %.6f\n', SSE);
fprintf('RMSE = %.6f mm^3\n\n', RMSE);


%% FIGURE 1: tumour and immune trajectories

figure('Position', [100 100 780 650]);

tiledlayout(2,1, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

% Tumour
nexttile;

plot(t_model, C_model, ...
    'Color', col.total, ...
    'LineWidth', 3);

hold on;

plot(t_data, C_data, 'o', ...
    'Color', col.data, ...
    'MarkerFaceColor', 'w', ...
    'MarkerSize', 8, ...
    'LineWidth', 1.8);

ylabel('Tumour volume $C(t)$ (mm$^3$)');

title('Single population tumour-immune model');

legend( ...
    'Model', ...
    'Experimental data', ...
    'Location', 'northwest');

xlim([4 25]);

% T cells
nexttile;

plot(t_model, T_model, ...
    'Color', col.T, ...
    'LineWidth', 3);

xlabel('Time (days)');
ylabel('Effective T cells $T(t)$');

xlim([4 25]);
ylim([0, 1.1*max(T_model)]);


% FIGURE 2: residuals

figure('Position', [100 100 760 400]);

stem(t_data, residuals, ...
    'filled', ...
    'Color', col.total, ...
    'LineWidth', 1.8);

hold on;

yline(0, '--', ...
    'Color', [0.3 0.3 0.3]);

xlabel('Time (days)');
ylabel('Residual (mm$^3$)');

title('Residuals: single population tumour-immune model');

xlim([4 25]);


%% FUNCTIONS

function error = tumourImmuneObjective(theta, ...
    t_data, C_data, y0, ...
    k, s, d, a)

    p = exp(theta);

    r = p(1);
    K = p(2);

    [~, y] = ode45( ...
        @(t,y) tumourImmuneRHS(t, y, ...
        r, K, k, s, d, a), ...
        t_data, y0);

    residuals = y(:,1) - C_data;

    error = sum(residuals.^2);

end


function dydt = tumourImmuneRHS(~, y, ...
    r, K, k, s, d, a)

    C = y(1);
    T = y(2);

    dCdt = r*C*(1 - C/K) ...
           - k*T*C;

    dTdt = s ...
           - d*T ...
           + a*C*T;

    dydt = [dCdt; dTdt];

end