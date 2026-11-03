%% PARAMETER SENSITIVITY

clear;
clc;
close all;
format longG;

col = honoursPlotStyle();
rng(1,'twister');

% Combination observations

data = readtable('Data_IVM_PD1_Actual.csv');

t_data = [4; 7; 11; 14; 18; 21; 25];
C_data = data{:,2};
C_data = C_data(:);

% Sequential reference parameters

base.rA = 0.365371;
base.rB = 0.365371;
base.KA = 58.772570;
base.KB = 413.826598;
base.eps_mut = 0.05;
base.kA = 1e-5;
base.kB = 1e-6;
base.s = 10;
base.d = 0.2;
base.a = 1e-4;
base.deltaIVM = 0.132649;
base.etaIVM = 0.08;
base.dD = 0.25;
base.pD = 0.05;
base.hD = 10;
base.muPD1 = 1048.95;
base.C0 = 10;
base.fB0 = 0.15;
base.IVM_start = 4;
base.IVM_end = 10;
base.PD1_start = 4;
base.PD1_end = 25;

paramNames = {'rA','rB','KA','KB','eps_mut','kA','kB','s','d','a','deltaIVM','etaIVM','dD','pD','hD','muPD1','fB0'};
prettyNames = {'$r_A$','$r_B$','$K_A$','$K_B$','$\varepsilon$','$k_A$','$k_B$','$s$','$d$','$a$','$\delta_{\mathrm{IVM}}$','$\eta_{\mathrm{IVM}}$','$d_D$','$p_D$','$h_D$','$\mu_{\mathrm{PD1}}$','$f_{B,0}$'};

p0 = cellfun(@(name) base.(name),paramNames);
nPar = numel(p0);

lowerBounds = 0.8*p0;
upperBounds = 1.2*p0;

rangeTable = table(string(paramNames(:)),p0(:),lowerBounds(:),upperBounds(:),'VariableNames',{'Parameter','Baseline','Lower','Upper'});

% Baseline simulation and trapezoidal AUC

t_grid = (4:0.1:25)';
t_PRCC = [7 11 14 18 21 25];

ode_options = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',0.1,'NonNegative',1:4);

[baseMetrics, C_base] = simulateCombination(base,t_grid,t_data,C_data,ode_options);

fprintf('\nBASELINE ADDITIVE COMBINATION\n');
fprintf('C(25) = %.6f mm^3\n',baseMetrics.C25);
fprintf('AUC = %.6f mm^3 day\n',baseMetrics.AUC);
fprintf('B fraction at day 25 = %.4f %%\n',baseMetrics.Bfrac25);
fprintf('RMSE after day 4 = %.6f mm^3\n',baseMetrics.RMSEafterIC);

% Central one-percent perturbations, varying each input independently

h_local = 0.01;

S_C25 = zeros(nPar,1);
S_AUC = zeros(nPar,1);
S_Bfrac = zeros(nPar,1);

for j = 1:nPar
    pPlus = base;
    pMinus = base;

    pPlus.(paramNames{j}) = p0(j)*(1+h_local);
    pMinus.(paramNames{j}) = p0(j)*(1-h_local);

    mPlus = simulateCombination(pPlus,t_grid,t_data,C_data,ode_options);
    mMinus = simulateCombination(pMinus,t_grid,t_data,C_data,ode_options);

    S_C25(j) = (mPlus.C25-mMinus.C25)/(2*h_local*baseMetrics.C25);
    S_AUC(j) = (mPlus.AUC-mMinus.AUC)/(2*h_local*baseMetrics.AUC);
    S_Bfrac(j) = (mPlus.Bfrac25-mMinus.Bfrac25)/(2*h_local*baseMetrics.Bfrac25);
end

localTable = table(string(paramNames(:)),S_C25,S_AUC,S_Bfrac,abs(S_C25),'VariableNames',{'Parameter','S_C25','S_AUC','S_Bfrac25','Abs_S_C25'});
localTable = sortrows(localTable,'Abs_S_C25','descend');

fprintf('\nLOCAL NORMALISED SENSITIVITIES\n');
disp(localTable);

% LHS draws over +/-20%; this permutation reproduces the original sample.

N = 1000;
X = zeros(N,nPar);

for j = 1:nPar
    u = ((0:N-1)'+rand(N,1))/N;
    [~, order] = sort(rand(N,1));
    X(:,j) = lowerBounds(j)+(upperBounds(j)-lowerBounds(j))*u(order);
end

Y_C25 = zeros(N,1);
Y_AUC = zeros(N,1);
Y_Bfrac25 = zeros(N,1);
Y_RMSE = zeros(N,1);
Y_RMSEafter = zeros(N,1);

C_store = zeros(N,numel(t_grid));
C_PRCC_time = zeros(N,numel(t_PRCC));

fprintf('\nRunning %d LHS simulations\n',N);

for n = 1:N
    p = base;

    for j = 1:nPar
        p.(paramNames{j}) = X(n,j);
    end

    [metrics, C] = simulateCombination(p,t_grid,t_data,C_data,ode_options);

    Y_C25(n) = metrics.C25;
    Y_AUC(n) = metrics.AUC;
    Y_Bfrac25(n) = metrics.Bfrac25;
    Y_RMSE(n) = metrics.RMSEall;
    Y_RMSEafter(n) = metrics.RMSEafterIC;

    C_store(n,:) = C';
    C_PRCC_time(n,:) = interp1(t_grid,C,t_PRCC,'linear');

    if mod(n,100)==0
        fprintf('Completed %d / %d simulations\n',n,N);
    end
end

% PRCCs for the three outputs and six observation times

[PRCC, pvalues] = computePRCC(X,[Y_C25 Y_AUC Y_Bfrac25 C_PRCC_time]);

PRCC_C25 = PRCC(:,1);
PRCC_AUC = PRCC(:,2);
PRCC_B = PRCC(:,3);

pval_C25 = pvalues(:,1);
pval_AUC = pvalues(:,2);
pval_B = pvalues(:,3);

PRCC_time = PRCC(:,4:end);
P_time = pvalues(:,4:end);

PRCC_table = table(string(paramNames(:)),PRCC_C25,pval_C25,PRCC_AUC,pval_AUC,PRCC_B,pval_B,max(abs(PRCC(:,1:2)),[],2),'VariableNames',{'Parameter','PRCC_C25','p_C25','PRCC_AUC','p_AUC','PRCC_Bfrac25','p_Bfrac25','MaxAbsTumourPRCC'});
PRCC_table = sortrows(PRCC_table,'MaxAbsTumourPRCC','descend');

timeTable = [table(string(paramNames(:)),'VariableNames',{'Parameter'}), array2table(PRCC_time,'VariableNames',{'Day7','Day11','Day14','Day18','Day21','Day25'})];

fprintf('\nGLOBAL SENSITIVITIES\n');
disp(PRCC_table);

fprintf('\nTIME-DEPENDENT PRCC\n');
disp(timeTable);

%% FIGURE 1: Time-dependent PRCC

figure('Position',[100 100 920 680]);

% Columns represent the six observation times, which are not equally spaced.

imagesc(1:numel(t_PRCC),1:nPar,PRCC_time);
set(gca,'YDir','normal','XTick',1:numel(t_PRCC),'XTickLabel',string(t_PRCC),'YTick',1:nPar,'YTickLabel',prettyNames,'TickLabelInterpreter','latex');

map = interp1([-1 0 1],[col.pd1; col.no; col.ivm],linspace(-1,1,256));
colormap(gca,map);
caxis([-1 1]);

cb = colorbar;
cb.Label.String = 'PRCC';

xlabel('Time (days)');
ylabel('Parameter');
title('Time-dependent PRCC for total tumour volume $C(t)$');

%% FUNCTIONS

function [metrics, C] = simulateCombination(p,t_grid,t_data,observations,options)

    % Recompute T0=s/d whenever either parameter changes.

    y0 = [p.C0*(1-p.fB0); p.C0*p.fB0; p.s/p.d; 0];
    [t, y] = ode45(@(t,y) combinationRHS(t,y,p),t_grid,y0,options);

    C = y(:,1)+y(:,2);
    residuals = interp1(t,C,t_data,'linear')-observations;

    metrics.C25 = C(end);
    metrics.AUC = trapz(t,C);
    metrics.Bfrac25 = 100*y(end,2)/C(end);
    metrics.RMSEall = sqrt(mean(residuals.^2));
    metrics.RMSEafterIC = sqrt(mean(residuals(2:end).^2));

end

function dydt = combinationRHS(t,y,p)

    A = y(1);
    B = y(2);
    T = y(3);
    D = y(4);
    C = A+B;

    uIVM = double(t>=p.IVM_start && t<=p.IVM_end);
    uPD1 = double(t>=p.PD1_start && t<=p.PD1_end);

    GA = p.rA*A*(1-C/p.KA);
    GB = p.rB*B*(1-C/p.KB);
    multiplier = 1+p.muPD1*uPD1;

    dA = (1-p.eps_mut)*GA+p.eps_mut*GB-multiplier*p.kA*T*A-p.deltaIVM*uIVM*A;
    dB = (1-p.eps_mut)*GB+p.eps_mut*GA-multiplier*p.kB*T*B-p.deltaIVM*uIVM*B;
    dT = p.s-p.d*T+p.a*C*T+p.pD*(D/(p.hD+D))*T;
    dD = p.etaIVM*uIVM*C-p.dD*D;

    dydt = [dA; dB; dT; dD];

end

function [prcc, pval] = computePRCC(X,Y)

    % Rank each variable, then remove the effects of the other parameters.

    [N, k] = size(X);
    XR = rankColumns(X);
    YR = rankColumns(Y);

    prcc = zeros(k,size(Y,2));
    pval = zeros(size(prcc));
    df = N-k-1;

    for j = 1:k
        Z = [ones(N,1), XR(:,setdiff(1:k,j))];

        rx = XR(:,j)-Z*(Z\XR(:,j));
        ry = YR-Z*(Z\YR);

        denominator = sqrt((rx'*rx)*sum(ry.^2,1));

        r = zeros(1,size(Y,2));
        valid = denominator>0;
        r(valid) = (rx'*ry(:,valid))./denominator(valid);
        r = max(-0.999999999,min(0.999999999,r));

        prcc(j,:) = r;

        tStatistic = abs(r).*sqrt(df./(1-r.^2));
        pval(j,:) = betainc(df./(df+tStatistic.^2),df/2,0.5);
    end

end

function ranks = rankColumns(X)

    % Tied values receive their average rank.

    [N, k] = size(X);
    ranks = zeros(N,k);

    for j = 1:k
        [values, order] = sort(X(:,j));
        [~, ~, group] = unique(values);
        averageRank = accumarray(group,(1:N)',[],@mean);
        ranks(order,j) = averageRank(group);
    end

end