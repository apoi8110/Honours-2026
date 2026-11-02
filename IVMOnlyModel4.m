%% IVM-ONLY TUMOUR-IMMUNE MODEL
% Fitting just delta_IVM with the remaining parameters fixed.

clear;
clc;
close all;

col = honoursPlotStyle();

% IVM-only data

data = readtable('Data_IVM.csv');

t_data = [4; 7; 11; 14; 18; 21; 25];
C_data = data{:,2};
C_data = C_data(:);

% Fixed tumour and immune parameters

r_A = 0.365371;
r_B = 0.365371;
K_A = 58.772570;
K_B = 413.826598;
eps_mut = 0.05;

k_A = 1e-5;
k_B = 1e-6;
s = 10;
d = 0.2;
a = 1e-4;

% ICD parameters and treatment window

eta_IVM = 0.08;
d_D = 0.25;
p_D = 0.05;
h_D = 10;

IVM_start = 4;
IVM_end = 10;

% Initial conditions

A0 = 8.5;
B0 = 1.5;
T0 = s/d;
D0 = 0;

y0 = [A0; B0; T0; D0];

% Numerical settings

ode_options = odeset( ...
    'RelTol', 1e-8, ...
    'AbsTol', 1e-10, ...
    'MaxStep', 0.05);

options = optimset( ...
    'Display', 'off', ...
    'MaxFunEvals', 5000, ...
    'MaxIter', 5000, ...
    'TolX', 1e-9, ...
    'TolFun', 1e-9);

% The same equations are used with and without IVM.

model = @(t,y,delta,treated) ivmRHS( ...
    t, y, delta, treated, ...
    r_A, K_A, r_B, K_B, eps_mut, ...
    k_A, k_B, s, d, a, eta_IVM, d_D, p_D, h_D, ...
    IVM_start, IVM_end);

% Fit delta_IVM in log parameters to retain positivity.

delta_guess = 0.08;
theta0 = log(delta_guess);

theta_fit = fminsearch( ...
    @(theta) ivmObjective(theta, ...
    t_data, C_data, y0, model, ode_options), ...
    theta0, options);

delta_IVM = exp(theta_fit);

% Fitted and matched untreated trajectories

t_plot = unique([linspace(4,25,600)'; IVM_end]);

[t_model, y_ivm] = ode45( ...
    @(t,y) model(t,y,delta_IVM,true), ...
    t_plot, y0, ode_options);

[~, y_no] = ode45( ...
    @(t,y) model(t,y,delta_IVM,false), ...
    t_plot, y0, ode_options);

A_ivm = y_ivm(:,1);
B_ivm = y_ivm(:,2);
T_ivm = y_ivm(:,3);
D_ivm = y_ivm(:,4);
C_ivm = A_ivm + B_ivm;

C_no = y_no(:,1) + y_no(:,2);
T_no = y_no(:,3);

% Predictions at observation times and reported quantities

[~, y_at_data] = ode45( ...
    @(t,y) model(t,y,delta_IVM,true), ...
    t_data, y0, ode_options);

C_at_data = y_at_data(:,1) + y_at_data(:,2);
residuals = C_at_data - C_data;

SSE = sum(residuals.^2);
RMSE = sqrt(mean(residuals.^2));
RMSE_after_IC = sqrt(mean(residuals(2:end).^2));

AUC_ivm = trapz(t_model, C_ivm);
AUC_no = trapz(t_model, C_no);
ICD_half_life = log(2)/d_D;

% Print results

fprintf('\nIVM-ONLY MODEL\n');
fprintf('delta_IVM = %.6f /day\n', delta_IVM);
fprintf('SSE = %.6f\n', SSE);
fprintf('RMSE = %.6f mm^3\n', RMSE);
fprintf('RMSE, days 7-25 = %.6f mm^3\n', RMSE_after_IC);
fprintf('ICD signal half-life = %.6f days\n', ICD_half_life);

fprintf('\nNo treatment C(25) = %.6f mm^3\n', C_no(end));
fprintf('IVM C(25) = %.6f mm^3\n', C_ivm(end));
fprintf('No treatment AUC = %.6f mm^3 day\n', AUC_no);
fprintf('IVM AUC = %.6f mm^3 day\n\n', AUC_ivm);

%% FIGURE 1: Immune and ICD dynamics

figure('Position', [100 100 800 650]);

tiledlayout(2,1, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

nexttile;

plot(t_model, T_no, '--', ...
    'Color', col.no, ...
    'LineWidth', 2.4);

hold on;

plot(t_model, T_ivm, ...
    'Color', col.T, ...
    'LineWidth', 3);

ylabel('Effective T cells $T(t)$');
title('Immune dynamics under IVM');

legend( ...
    'No treatment', ...
    'IVM', ...
    'Location', 'northwest');

xlim([4 25]);
ylim([0, 10*ceil(1.10*max([T_no; T_ivm])/10)]);

nexttile;

plot(t_model, D_ivm, ...
    'Color', col.ivm, ...
    'LineWidth', 3);

xlabel('Time (days)');
ylabel('DAMP signal $D(t)$');

xlim([4 25]);
ylim([0, 1.10*max(D_ivm)]);

%% FUNCTIONS

function error = ivmObjective(theta, ...
    t_data, C_data, y0, model, ode_options)

    delta_IVM = exp(theta);

    [~, y] = ode45( ...
        @(t,y) model(t,y,delta_IVM,true), ...
        t_data, y0, ode_options);

    C_model = y(:,1) + y(:,2);
    residuals = C_model - C_data;

    error = sum(residuals.^2);

end


function dydt = ivmRHS(t, y, deltaIVM, treated, ...
    rA, KA, rB, KB, eps_mut, ...
    kA, kB, s, d, a, etaIVM, dD, pD, hD, ...
    IVM_start, IVM_end)

    A = y(1);
    B = y(2);
    T = y(3);
    D = y(4);

    C = A + B;
    uIVM = treated && t >= IVM_start && t < IVM_end;

    G_A = rA*A*(1 - C/KA);
    G_B = rB*B*(1 - C/KB);

    dAdt = (1-eps_mut)*G_A + eps_mut*G_B ...
        - kA*T*A - deltaIVM*uIVM*A;

    dBdt = (1-eps_mut)*G_B + eps_mut*G_A ...
        - kB*T*B - deltaIVM*uIVM*B;

    dTdt = s - d*T + a*C*T + pD*D/(hD+D)*T;
    dDdt = etaIVM*uIVM*C - dD*D;

    dydt = [dAdt; dBdt; dTdt; dDdt];

end