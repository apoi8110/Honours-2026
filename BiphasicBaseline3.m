%% FIT BIPHASIC MUTATION-UPON-DIVISION BASELINE MODEL
% Fits: Fitting r_A, K_A, r_B, K_B
% Fixing: eps_mut, initial A/B composition,immune parameters

%       r_A > r_B
%       K_A < K_B


clear;
clc;
close all;

col = honoursPlotStyle();

% Load untreated tumour data

data = readtable('Data.csv');

t_data = [4; 7; 11; 14; 18; 21; 25];

C_data = data{:,2};

t_data = t_data(:);
C_data = C_data(:);


% Initial tumour composition
% Keeping total initial tumour volume equal to 10 mm^3:
C0 = 10;

% Baseline assumption: 15% population B and 85% population A - initially
fB0 = 0.15;

B0 = fB0*C0;
A0 = (1-fB0)*C0;

fprintf('\nInitial tumour composition:\n');
fprintf('A0 = %.4f mm^3\n',A0);
fprintf('B0 = %.4f mm^3\n',B0);
fprintf('B fraction = %.1f %%\n',100*fB0);


% Fixing henotypic mixing fraction
eps_mut = 0.05;


% Fixing immune parameters 

k_A = 1e-5;
k_B = 1e-6;

% Baseline effective T-cell input.
s = 10;

% Effective T-cell loss rate.
d = 0.2;

% Tumour-driven immune stimulation.
a = 1e-4;

% Baseline immune population.
T0 = s/d;

y0 = [A0; B0; T0];

% Biologically:
% A = faster-growing, lower-capacity phenotype
% B = slower-growing, higher-capacity phenotype

rA_guess = 0.30;
rB_guess = 0.18;

KA_guess = 100;
KB_guess = 500;

% transforming parameters 

theta0 = log([rB_guess, ...
    rA_guess-rB_guess, ...
    KA_guess, ...
    KB_guess-KA_guess]);


% FIT USING FMINSEARCH

options = optimset('Display','iter', ...
    'MaxFunEvals',20000, ...
    'MaxIter',10000, ...
    'TolX',1e-9, ...
    'TolFun',1e-9);

theta_fit = fminsearch(@(theta) biphasicMutationObjective( ...
        theta, ...
        t_data,C_data,y0, ...
        eps_mut, ...
        k_A,k_B,s,d,a), ...
    theta0, ...
    options);


% EXTRACT FITTED PARAMETERS

r_B = exp(theta_fit(1));

delta_r = exp(theta_fit(2));

r_A = r_B + delta_r;

K_A = exp(theta_fit(3));

delta_K = exp(theta_fit(4));

K_B = K_A + delta_K;


% SOLVE FITTED MODEL ON SMOOTH GRID


t_plot = linspace(4,25,500);

[t_model,y_model] = ode45( ...
    @(t,y) biphasicMutationRHS( ...
        t,y, ...
        r_A,K_A,r_B,K_B, ...
        eps_mut, ...
        k_A,k_B,s,d,a), ...
    t_plot,y0);

A_model = y_model(:,1);
B_model = y_model(:,2);
T_model = y_model(:,3);

C_model = A_model+B_model;


% Model outputs

[~,y_at_data] = ode45( ...
    @(t,y) biphasicMutationRHS( ...
        t,y, ...
        r_A,K_A,r_B,K_B, ...
        eps_mut, ...
        k_A,k_B,s,d,a), ...
    t_data,y0);

A_at_data = y_at_data(:,1);
B_at_data = y_at_data(:,2);

C_at_data = A_at_data+B_at_data;

% Fit quality and error measurement

residuals = C_at_data-C_data;

SSE = sum(residuals.^2);

RMSE = sqrt(mean(residuals.^2));

% Results - parameter fits 


fprintf('\nFitted tumour parameters:\n');

fprintf('r_A = %.6f /day\n',r_A);
fprintf('r_B = %.6f /day\n',r_B);

fprintf('K_A = %.6f mm^3\n',K_A);
fprintf('K_B = %.6f mm^3\n',K_B);

fprintf('\nFixed phenotypic switching fraction:\n');
fprintf('epsilon = %.4f\n',eps_mut);

fprintf('\nInitial conditions:\n');
fprintf('A0 = %.6f mm^3\n',A0);
fprintf('B0 = %.6f mm^3\n',B0);
fprintf('T0 = %.6f\n',T0);

fprintf('\nFit quality:\n');
fprintf('SSE  = %.6f\n',SSE);
fprintf('RMSE = %.6f mm^3\n',RMSE);

fprintf('\nFinal populations:\n');
fprintf('A(25) = %.6f mm^3\n',A_model(end));
fprintf('B(25) = %.6f mm^3\n',B_model(end));
fprintf('C(25) = %.6f mm^3\n',C_model(end));

fprintf('\nFinal tumour composition:\n');

fprintf('A fraction = %.2f %%\n', ...
    100*A_model(end)/C_model(end));

fprintf('B fraction = %.2f %%\n', ...
    100*B_model(end)/C_model(end));

% FIGURE 1: FITTED BIPHASIC MODEL

figure('Position',[100 100 780 500]);

plot(t_model,C_model, ...
    'Color',col.total, ...
    'LineWidth',3.2);

hold on;

plot(t_model,A_model,'--', ...
    'Color',col.A, ...
    'LineWidth',2.4);

plot(t_model,B_model,'--', ...
    'Color',col.B, ...
    'LineWidth',2.4);

plot(t_data,C_data,'o', ...
    'Color',col.data, ...
    'MarkerFaceColor','w', ...
    'MarkerSize',8, ...
    'LineWidth',1.8);

xlabel('Time (days)');
ylabel('Tumour volume (mm$^3$)');

title('Fitted biphasic tumour immune model');

legend('Total tumour $C=A+B$', ...
    'Population $A$', ...
    'Population $B$', ...
    'Experimental data', ...
    'Location','northwest');

xlim([4 25]);

ylim([0 1.10*max([C_model;C_data])]);

% FIGURE 2: SEMILOG PLOT - SHOWS BIPHASIC DYNAMICS

figure('Position',[100 100 780 500]);

semilogy(t_model,max(C_model,1e-8), ...
    'Color',col.total, ...
    'LineWidth',3.2);

hold on;

semilogy(t_model,max(A_model,1e-8),'--', ...
    'Color',col.A, ...
    'LineWidth',2.4);

semilogy(t_model,max(B_model,1e-8),'--', ...
    'Color',col.B, ...
    'LineWidth',2.4);

semilogy(t_data,max(C_data,1e-8),'o', ...
    'Color',col.data, ...
    'MarkerFaceColor','w', ...
    'MarkerSize',8, ...
    'LineWidth',1.8);

xlabel('Time (days)');
ylabel('Tumour volume - log scale (mm$^3$)');

title('Biphasic tumour dynamics on a semilog plot');

legend('$C=A+B$', ...
    '$A$', ...
    '$B$', ...
    'Experimental data', ...
    'Location','northwest');

xlim([4 25]);


% OBJECTIVE FUNCTION

function error = biphasicMutationObjective( ...
    theta, ...
    t_data,C_data,y0, ...
    eps_mut, ...
    kA,kB,s,d,a)

    % Recovering biologically constrained parameters

    % Slower population B growth rate
    rB = exp(theta(1));

    % Positive difference ensures rA > rB
    delta_r = exp(theta(2));

    rA = rB + delta_r;

    % Lower carrying capacity
    KA = exp(theta(3));

    % Positive difference ensures KB > KA
    delta_K = exp(theta(4));

    KB = KA + delta_K;

    % Solve model at experimental times

    [~,y] = ode45( ...
        @(t,y) biphasicMutationRHS( ...
            t,y, ...
            rA,KA,rB,KB, ...
            eps_mut, ...
            kA,kB,s,d,a), ...
        t_data,y0);

    % Total tumour prediction

    A = y(:,1);
    B = y(:,2);

    C_model = A+B;

    % Residual sum of squares

    residuals = C_model-C_data;

    error = sum(residuals.^2);

end


% BIPHASIC MUTATION-UPON-DIVISION ODE MODEL

function dydt = biphasicMutationRHS( ...
    ~,y, ...
    rA,KA,rB,KB, ...
    eps_mut, ...
    kA,kB,s,d,a)

    % State variables

    A = y(1);
    B = y(2);
    T = y(3);

    C = A+B;


    % Logistic growth of each tumour phenotype

    GA = rA*A*(1-C/KA);

    GB = rB*B*(1-C/KB);


    % Population A
    dAdt = (1-eps_mut)*GA ...
        + eps_mut*GB ...
        - kA*T*A;


    % Population B
    dBdt = (1-eps_mut)*GB ...
        + eps_mut*GA ...
        - kB*T*B;


    % Effective T-cell compartment

    dTdt = s ...
        - d*T ...
        + a*C*T;


    % Return derivatives

    dydt = [dAdt;
        dBdt;
        dTdt];

end

