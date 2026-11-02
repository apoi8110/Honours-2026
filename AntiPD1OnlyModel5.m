%% ANTI-PD-1-ONLY TUMOUR-IMMUNE MODEL

% Fitting mu_PD1 and as usual, keeping remaining parameters fixed.

clear;
clc;
close all;

col = honoursPlotStyle();

% Anti-PD-1 data

data = readtable('Data_PD1_Actual.csv');
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

% Treatment window and initial conditions

PD1_start = 4;
PD1_end = 25;

A0 = 8.5;
B0 = 1.5;
T0 = s/d;
y0 = [A0; B0; T0];

% Numerical settings

ode_options = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.05);
options = optimset('Display','off','MaxFunEvals',5000,'MaxIter',5000,'TolX',1e-9,'TolFun',1e-9);

model = @(t,y,mu,treated) pd1RHS(t,y,mu,treated,r_A,K_A,r_B,K_B,eps_mut,k_A,k_B,s,d,a,PD1_start,PD1_end);

% Fit mu_PD1 in log parameters to retain positivity.

mu_guess = 20;
theta0 = log(mu_guess);
theta_fit = fminsearch(@(theta) pd1Objective(theta,t_data,C_data,y0,model,ode_options),theta0,options);
mu_PD1 = exp(theta_fit);

% Fitted and matched untreated trajectories

t_plot = linspace(4,25,600);

[t_model, y_pd1] = ode45(@(t,y) model(t,y,mu_PD1,true),t_plot,y0,ode_options);
[~, y_no] = ode45(@(t,y) model(t,y,mu_PD1,false),t_plot,y0,ode_options);

A_pd1 = y_pd1(:,1);
B_pd1 = y_pd1(:,2);
C_pd1 = A_pd1 + B_pd1;
C_no = y_no(:,1) + y_no(:,2);

% Predictions at observation times and reported quantities

[~, y_at_data] = ode45(@(t,y) model(t,y,mu_PD1,true),t_data,y0,ode_options);

C_at_data = y_at_data(:,1) + y_at_data(:,2);
residuals = C_at_data - C_data;

SSE = sum(residuals.^2);
RMSE = sqrt(mean(residuals.^2));
RMSE_after_IC = sqrt(mean(residuals(2:end).^2));

AUC_pd1 = trapz(t_model,C_pd1);
AUC_no = trapz(t_model,C_no);

killing_multiplier = 1 + mu_PD1;
k_A_eff = k_A*killing_multiplier;
k_B_eff = k_B*killing_multiplier;
killing_ratio = k_B_eff/k_A_eff;

% Print results

fprintf('\nANTI-PD-1-ONLY MODEL\n');
fprintf('mu_PD1 = %.6f\n',mu_PD1);
fprintf('Killing multiplier = %.6f\n',killing_multiplier);
fprintf('k_A_eff = %.6e\n',k_A_eff);
fprintf('k_B_eff = %.6e\n',k_B_eff);
fprintf('k_B_eff / k_A_eff = %.6f\n',killing_ratio);
fprintf('SSE = %.6f\n',SSE);
fprintf('RMSE = %.6f mm^3\n',RMSE);
fprintf('RMSE, days 7-25 = %.6f mm^3\n',RMSE_after_IC);

fprintf('\nNo treatment C(25) = %.6f mm^3\n',C_no(end));
fprintf('Anti-PD-1 C(25) = %.6f mm^3\n',C_pd1(end));
fprintf('Observed anti-PD-1 C(25) = %.6f mm^3\n',C_data(end));
fprintf('No treatment AUC = %.6f mm^3 day\n',AUC_no);
fprintf('Anti-PD-1 AUC = %.6f mm^3 day\n\n',AUC_pd1);

%% FIGURE 1: Anti-PD-1 tumour fit

figure('Position',[100 100 800 520]);

plot(t_model,C_pd1,'Color',col.pd1,'LineWidth',3);
hold on;
plot(t_model,A_pd1,'--','Color',col.A,'LineWidth',2.3);
plot(t_model,B_pd1,'--','Color',col.B,'LineWidth',2.3);
plot(t_data,C_data,'o','Color',col.data,'MarkerFaceColor','w','MarkerSize',8,'LineWidth',1.8);

xlabel('Time (days)');
ylabel('Tumour volume (mm$^3$)');
title('Biphasic tumour response to anti-PD-1');
legend('Fitted total tumour $C=A+B$','Population $A$','Population $B$','Experimental data','Location','northwest');

xlim([4 25]);
ylim([0, 1.10*max([C_pd1; C_data])]);

%% FUNCTIONS

function error = pd1Objective(theta,t_data,C_data,y0,model,ode_options)

    mu_PD1 = exp(theta);

    [~, y] = ode45(@(t,y) model(t,y,mu_PD1,true),t_data,y0,ode_options);

    C_model = y(:,1) + y(:,2);
    residuals = C_model - C_data;
    error = sum(residuals.^2);

end

function dydt = pd1RHS(t,y,muPD1,treated,rA,KA,rB,KB,eps_mut,kA,kB,s,d,a,PD1_start,PD1_end)

    A = y(1);
    B = y(2);
    T = y(3);
    C = A+B;

    uPD1 = treated && t>=PD1_start && t<=PD1_end;
    m = 1+muPD1*uPD1;

    G_A = rA*A*(1-C/KA);
    G_B = rB*B*(1-C/KB);

    dAdt = (1-eps_mut)*G_A + eps_mut*G_B - m*kA*T*A;
    dBdt = (1-eps_mut)*G_B + eps_mut*G_A - m*kB*T*B;
    dTdt = s-d*T+a*C*T;

    dydt = [dAdt; dBdt; dTdt];

end