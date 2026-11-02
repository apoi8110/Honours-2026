%% ADDITIVE IVM + ANTI-PD1 MODEL

% Use the single-treatment estimates - no fitting extra data, 
% Want to test the fit

clear;
clc;
close all;

col = honoursPlotStyle();

% Experimental data

t_data = [4; 7; 11; 14; 18; 21; 25];

data_no = readtable('Data.csv');
data_ivm = readtable('Data_IVM.csv');
data_pd1 = readtable('Data_PD1_Actual.csv');
data_combo = readtable('Data_IVM_PD1_Actual.csv');

C_data_no = data_no{:,2};
C_data_ivm = data_ivm{:,2};
C_data_pd1 = data_pd1{:,2};
C_data_combo = data_combo{:,2};
C_data = [C_data_no, C_data_ivm, C_data_pd1, C_data_combo];

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

% Previously calibrated treatment parameters

delta_IVM = 0.132649;
eta_IVM = 0.08;
d_D = 0.25;
p_D = 0.05;
h_D = 10;
mu_PD1 = 1048.95;

IVM_start = 4;
IVM_end = 10;
PD1_start = 4;
PD1_end = 25;

% Initial conditions and numerical settings

A0 = 8.5;
B0 = 1.5;
T0 = s/d;
D0 = 0;
y0 = [A0; B0; T0; D0];

ode_options = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.05);
t_plot = unique([linspace(4,25,600)'; t_data; IVM_end]);

% The same model is used with the two treatments switched on or off.

model = @(t,y,useIVM,usePD1) combinationRHS(t,y,useIVM,usePD1,r_A,K_A,r_B,K_B,eps_mut,k_A,k_B,s,d,a,...
    delta_IVM,eta_IVM,d_D,p_D,h_D,mu_PD1,IVM_start,IVM_end,PD1_start,PD1_end);

[t_model, y_no] = ode45(@(t,y) model(t,y,false,false),t_plot,y0,ode_options);
[~, y_ivm] = ode45(@(t,y) model(t,y,true,false),t_plot,y0,ode_options);
[~, y_pd1] = ode45(@(t,y) model(t,y,false,true),t_plot,y0,ode_options);
[~, y_combo] = ode45(@(t,y) model(t,y,true,true),t_plot,y0,ode_options);

C_no = y_no(:,1) + y_no(:,2);
C_ivm = y_ivm(:,1) + y_ivm(:,2);
C_pd1 = y_pd1(:,1) + y_pd1(:,2);
C_combo = y_combo(:,1) + y_combo(:,2);
C_model = [C_no, C_ivm, C_pd1, C_combo];

A_combo = y_combo(:,1);
B_combo = y_combo(:,2);
T_combo = y_combo(:,3);

% Errors, endpoint volumes and cumulative burden

[~, data_rows] = ismember(t_data,t_model);
C_at_data = C_model(data_rows,:);
residuals = C_at_data - C_data;

SSE = sum(residuals.^2,1)';
RMSE = sqrt(SSE/numel(t_data));
RMSE_after_IC = sqrt(mean(residuals(2:end,:).^2,1))';
C25 = C_model(end,:)';
AUC = trapz(t_model,C_model)';

Treatment = {'No treatment'; 'IVM'; 'Anti-PD1'; 'Combination'};
results = table(Treatment,SSE,RMSE,RMSE_after_IC,C25,AUC);

residuals_combo = residuals(:,4);
final_residual = residuals_combo(end);
final_residual_share = 100*final_residual^2/SSE(4);
endpoint_reduction = 100*(1-C25(4)/C25(1));
AUC_reduction = 100*(1-AUC(4)/AUC(1));

% Composition and growth/removal terms

q_no = y_no(:,2)./C_no;
q_combo = B_combo./C_combo;
day10 = find(t_model==10,1);

growth = r_A*A_combo.*(1-C_combo/K_A) + r_B*B_combo.*(1-C_combo/K_B);
multiplier = 1 + mu_PD1*(t_model>=PD1_start & t_model<=PD1_end);
immune_loss = multiplier.*T_combo.*(k_A*A_combo + k_B*B_combo);
IVM_loss = delta_IVM*(t_model>=IVM_start & t_model<IVM_end).*C_combo;
net_rate = growth - immune_loss - IVM_loss;

% Print results

fprintf('\nADDITIVE COMBINATION MODEL\n');
disp(results);

fprintf('Observed combination C(25) = %.6f mm^3\n',C_data_combo(end));
fprintf('Final combination residual = %.6f mm^3\n',final_residual);
fprintf('Final residual share of SSE = %.2f %%\n',final_residual_share);
fprintf('Combination endpoint reduction = %.2f %%\n',endpoint_reduction);
fprintf('Combination AUC reduction = %.2f %%\n',AUC_reduction);

fprintf('\nUntreated B/C at day 10 = %.6f\n',q_no(day10));
fprintf('Combination B/C at day 10 = %.6f\n',q_combo(day10));
fprintf('Untreated B fraction at day 25 = %.2f %%\n',100*q_no(end));
fprintf('Combination B fraction at day 25 = %.2f %%\n',100*q_combo(end));

fprintf('\nCombination rates just before day 25:\n');
fprintf('Growth contribution = %.6f mm^3/day\n',growth(end));
fprintf('Immune removal = %.6f mm^3/day\n',immune_loss(end));
fprintf('Direct IVM removal = %.6f mm^3/day\n',IVM_loss(end));
fprintf('Net tumour rate = %.6f mm^3/day\n\n',net_rate(end));

%% FIGURE 1: Four-treatment comparison

figure('Position',[100 100 900 560]);

plot(t_model,C_no,'Color',col.no,'LineWidth',2.8);
hold on;
plot(t_model,C_ivm,'Color',col.ivm,'LineWidth',3);
plot(t_model,C_pd1,'Color',col.pd1,'LineWidth',3);
plot(t_model,C_combo,'Color',col.combo,'LineWidth',3.2);

plot(t_data,C_data_no,'o','Color',col.no,'MarkerFaceColor','w','MarkerSize',7,'LineWidth',1.5,'HandleVisibility','off');
plot(t_data,C_data_ivm,'o','Color',col.ivm,'MarkerFaceColor','w','MarkerSize',7,'LineWidth',1.5,'HandleVisibility','off');
plot(t_data,C_data_pd1,'o','Color',col.pd1,'MarkerFaceColor','w','MarkerSize',7,'LineWidth',1.5,'HandleVisibility','off');
plot(t_data,C_data_combo,'o','Color',col.combo,'MarkerFaceColor','w','MarkerSize',7,'LineWidth',1.5,'HandleVisibility','off');

xlabel('Time (days)');
ylabel('Total tumour volume $C(t)$ (mm$^3$)');
title('Comparison of treatment responses');
subtitle('Solid lines: model predictions; circles: experimental data');
legend('No treatment','IVM only','Anti-PD-1 only','IVM + anti-PD-1','Location','northwest');

xlim([4 25]);
ylim([0, 1.10*max([C_model(:); C_data(:)])]);

%% FIGURE 2: Combination prediction

figure('Position',[100 100 820 520]);

plot(t_model,C_combo,'Color',col.combo,'LineWidth',3.2);
hold on;
plot(t_model,A_combo,'--','Color',col.A,'LineWidth',2.2);
plot(t_model,B_combo,'--','Color',col.B,'LineWidth',2.2);
plot(t_data,C_data_combo,'o','Color',col.data,'MarkerFaceColor','w','MarkerSize',8,'LineWidth',1.8);

xlabel('Time (days)');
ylabel('Tumour volume (mm$^3$)');
title('IVM + anti-PD-1 combination prediction');
legend('Predicted total tumour $C=A+B$','Population $A$','Population $B$','Combination data','Location','northwest');

xlim([4 25]);
ylim([0, 1.10*max([C_combo; C_data_combo])]);

%% FUNCTION

function dydt = combinationRHS(t,y,useIVM,usePD1,rA,KA,rB,KB,eps_mut,kA,kB, ...
    s,d,a,deltaIVM,etaIVM,dD,pD,hD,muPD1,IVM_start,IVM_end,PD1_start,PD1_end)

    A = y(1);
    B = y(2);
    T = y(3);
    D = y(4);
    C = A+B;

    uIVM = useIVM && t>=IVM_start && t<IVM_end;
    uPD1 = usePD1 && t>=PD1_start && t<=PD1_end;

    G_A = rA*A*(1-C/KA);
    G_B = rB*B*(1-C/KB);
    multiplier = 1+muPD1*uPD1;

    dAdt = (1-eps_mut)*G_A + eps_mut*G_B - multiplier*kA*T*A - deltaIVM*uIVM*A;
    dBdt = (1-eps_mut)*G_B + eps_mut*G_A - multiplier*kB*T*B - deltaIVM*uIVM*B;
    dTdt = s-d*T+a*C*T+pD*D/(hD+D)*T;
    dDdt = etaIVM*uIVM*C-dD*D;

    dydt = [dAdt; dBdt; dTdt; dDdt];

end