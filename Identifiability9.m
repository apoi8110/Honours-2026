%% PARAMETER IDENTIFIABILITY

clear;
clc;
close all;
format longG;

col = honoursPlotStyle();

% Use the fitted profiles from the thesis and recalculate predictions.

S = load('Identifiability_Ordered_Checked.mat', 'fixed', 'thetaNames', 'thetaMLE', 'lowerBounds', 'upperBounds', 'profiles', 'costMin');

fixed = S.fixed;
paramNames = cellstr(string(S.thetaNames(:)))';
thetaJoint = S.thetaMLE(:)';
profiles = S.profiles;
costMin = S.costMin;

thetaSequential = [0.365371, 0.365371, 413.826598, 0.05, 0.15, 1e-6, 10, 0.2, 0.132649, 1048.95];

prettyNames = {'$r_A$', '$r_B$', '$K_B$', '$\varepsilon$', '$f_{B,0}$', '$k_B$', '$s$', '$d$', '$\delta_{\mathrm{IVM}}$', '$\mu_{\mathrm{PD1}}$'};

% Observations and numerical settings

t_data = [4; 7; 11; 14; 18; 21; 25];
t_plot = unique([(4:0.1:25)'; t_data; 10]);

[~, dataRows] = ismember(t_data,t_plot);

fitRows = 2:numel(t_data);
nFit = numel(fitRows);

files = {'Data.csv', 'Data_IVM.csv', 'Data_PD1_Actual.csv', 'Data_IVM_PD1_Actual.csv'};
observations = zeros(7,4);

for j = 1:4
    data = readmatrix(files{j});
    observations(:,j) = data(isfinite(data(:,2)),2);
end

ode_options = odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',0.1);

% Sequential and joint predictions; combination data remain withheld.

parameterVectors = [thetaSequential; thetaJoint];
vectorNames = {'Sequential', 'Joint'};
armNames = {'Untreated', 'IVM', 'AntiPD1', 'Combination'};

uIVM = [0 1 0 1];
uPD1 = [0 0 1 1];

states = zeros(numel(t_plot),5,4,2);
metricRows = cell(8,7);
calibrationSSE = zeros(2,3);
row = 0;

for k = 1:2
    p = fixed;

    for j = 1:numel(paramNames)
        p.(paramNames{j}) = parameterVectors(k,j);
    end

    for g = 1:4
        Y = solveArm(p,t_plot,uIVM(g),uPD1(g),ode_options);
        states(:,:,g,k) = Y;

        C = Y(:,1)+Y(:,2);
        residuals = C(dataRows)-observations(:,g);

        % Exclude day 4, which sets the initial condition.

        SSE6 = sum(residuals(fitRows).^2);
        RMSE6 = sqrt(SSE6/nFit);

        row = row+1;
        metricRows(row,:) = {vectorNames{k}, armNames{g}, SSE6, RMSE6, C(end), Y(end,5), Y(end,2)/C(end)};

        if g<=3
            calibrationSSE(k,g) = SSE6;
        end
    end
end

predictionTable = cell2table(metricRows,'VariableNames',{'ParameterVector','Arm','SSE6','RMSE6','C25','AUC','q25'});

% Each calibration group has its own profiled residual variance.

variances = calibrationSSE/nFit;
profiledCost = sum(nFit*log(variances),2);

likelihoodTable = table(string(vectorNames(:)),variances(:,1),variances(:,2),variances(:,3),profiledCost,'VariableNames',{'ParameterVector','VarianceUntreated','VarianceIVM','VarianceAntiPD1','ProfiledCost'});

fprintf('\nSEQUENTIAL AND JOINT PREDICTIONS\n');
disp(predictionTable);

fprintf('\nCALIBRATION LIKELIHOOD\n');
fprintf('Saved joint objective = %.10f\n',costMin);
disp(likelihoodTable);

% Interpolated profile crossings; rA was not profiled.

profileIdx = 2:10;
threshold = 3.841459;
referenceTolerance = 5e-4;

lowerCrossing = nan(10,1);
upperCrossing = nan(10,1);
lowerFlagged = false(10,1);
upperFlagged = false(10,1);

for j = profileIdx
    P = profiles.(paramNames{j});

    [lowerCrossing(j), upperCrossing(j), lowerFlagged(j), upperFlagged(j)] = profileCrossings(P.grid(:),P.cost(:)-costMin,P.exitFlags(:),thetaJoint(j),threshold,j~=5);
end

parameterTable = table(string(paramNames(:)),thetaJoint(:),S.lowerBounds(:),S.upperBounds(:),lowerCrossing,upperCrossing,lowerFlagged,upperFlagged,'VariableNames',{'Parameter','Estimate','LowerBound','UpperBound','LowerCrossing','UpperCrossing','LowerBracketNonconverged','UpperBracketNonconverged'});

fprintf('\nPARAMETER ESTIMATES AND PROFILE CROSSINGS\n');
fprintf('NaN means no crossing was found, or the parameter was not profiled.\n');
disp(parameterTable);

%% FIGURE 1: Profile likelihoods

figure('Position',[50 50 1200 850]);
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');

for j = profileIdx
    P = profiles.(paramNames{j});
    delta = P.cost(:)-costMin;
    nonconverged = P.exitFlags(:)<=0;

    nexttile;

    plot(P.grid,delta,'o-','Color',col.combo,'LineWidth',1.8);
    hold on;

    plot(P.grid(nonconverged),delta(nonconverged),'x','Color',col.data,'LineWidth',1.8,'MarkerSize',9);
    yline(threshold,'--','Color',col.ivm,'LineWidth',1.2);
    xline(thetaJoint(j),':','Color',col.data,'LineWidth',1.2);

    if j~=5
        set(gca,'XScale','log');
    end

    xlabel(prettyNames{j});
    ylabel('$\Delta J$');
    title(paramNames{j},'Interpreter','none');

    grid on;
end

sgtitle('Profile likelihoods: crosses mark nonconverged searches');

%% FIGURE 2: Parameter compensation

xParameter = [10 7 8 2];
yParameter = [6 10 10 3];

figure('Position',[100 60 1050 720]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

for k = 1:4
    P = profiles.(paramNames{xParameter(k)});
    delta = P.cost(:)-costMin;

    lowLoss = delta<=threshold & delta>=-referenceTolerance;
    nonconverged = P.exitFlags(:)<=0;

    nexttile;

    loglog(P.grid,P.optimisedTheta(:,yParameter(k)),'-','Color',col.no,'LineWidth',1.5,'DisplayName','All stored points');
    hold on;

    loglog(P.grid(lowLoss),P.optimisedTheta(lowLoss,yParameter(k)),'o','Color',col.combo,'MarkerFaceColor',col.combo,'DisplayName','Below LR reference');
    loglog(P.grid(nonconverged),P.optimisedTheta(nonconverged,yParameter(k)),'x','Color',col.data,'LineWidth',1.5,'DisplayName','Nonconverged search');

    xlabel(prettyNames{xParameter(k)});
    ylabel(['Optimised ',prettyNames{yParameter(k)}]);

    grid on;

    if k==1
        legend('Location','best');
    end
end

%% FUNCTIONS

function Y = solveArm(p,t_eval,uIVM,uPD1,options)

    % Restart when IVM stops; the fifth state integrates tumour burden.

    y0 = [p.C0*(1-p.fB0); p.C0*p.fB0; p.s/p.d; 0; 0];

    options1 = odeset(options,'Jacobian',@(t,y) modelJacobian(y,p,uIVM,uPD1));
    options2 = odeset(options,'Jacobian',@(t,y) modelJacobian(y,p,0,uPD1));

    sol1 = ode15s(@(t,y) modelRHS(y,p,uIVM,uPD1),[4 10],y0,options1);
    sol2 = ode15s(@(t,y) modelRHS(y,p,0,uPD1),[10 25],sol1.y(:,end),options2);

    first = t_eval<=10;

    Y = zeros(numel(t_eval),5);
    Y(first,:) = deval(sol1,t_eval(first))';
    Y(~first,:) = deval(sol2,t_eval(~first))';

end

function dydt = modelRHS(y,p,uIVM,uPD1)

    A = y(1);
    B = y(2);
    T = y(3);
    D = y(4);
    C = A+B;

    GA = p.rA*A*(1-C/p.KA);
    GB = p.rB*B*(1-C/p.KB);
    m = 1+p.muPD1*uPD1;

    dA = (1-p.eps_mut)*GA+p.eps_mut*GB-m*p.kA*T*A-p.deltaIVM*uIVM*A;
    dB = p.eps_mut*GA+(1-p.eps_mut)*GB-m*p.kB*T*B-p.deltaIVM*uIVM*B;
    dT = p.s-p.d*T+p.a*C*T+p.pD*D/(p.hD+D)*T;
    dD = p.etaIVM*uIVM*C-p.dD*D;

    dydt = [dA; dB; dT; dD; C];

end

function J = modelJacobian(y,p,uIVM,uPD1)

    A = y(1);
    B = y(2);
    T = y(3);
    D = y(4);
    C = A+B;
    m = 1+p.muPD1*uPD1;

    growth = [p.rA*(1-(2*A+B)/p.KA), -p.rA*A/p.KA; -p.rB*B/p.KB, p.rB*(1-(A+2*B)/p.KB)];
    mixing = [1-p.eps_mut, p.eps_mut; p.eps_mut, 1-p.eps_mut];

    J = zeros(5);

    J(1:2,1:2) = mixing*growth-diag([m*p.kA*T+p.delta*uIVM, m*p.kB*T+p.delta*uIVM]);
    J(1:2,3) = [-m*p.kA*A; -m*p.kB*B];

    J(3,:) = [p.a*T, p.a*T, -p.d+p.a*C+p.pD*D/(p.hD+D), p.pD*p.hD*T/(p.hD+D)^2, 0];
    J(4,:) = [p.etaIVM*uIVM, p.etaIVM*uIVM, 0, -p.dD, 0];
    J(5,1:2) = 1;

end

function [low, high, lowFlag, highFlag] = profileCrossings(grid,delta,flags,estimate,threshold,isLog)

    low = NaN;
    high = NaN;
    lowFlag = false;
    highFlag = false;

    [~, centre] = min(abs(grid-estimate));

    if isLog
        grid = log(grid);
    end

    for j = centre-1:-1:1
        if delta(j)>=threshold && delta(j+1)<threshold
            low = interp1(delta(j:j+1),grid(j:j+1),threshold);
            lowFlag = any(flags(j:j+1)<=0);
            break;
        end
    end

    for j = centre+1:numel(grid)
        if delta(j)>=threshold && delta(j-1)<threshold
            high = interp1(delta(j-1:j),grid(j-1:j),threshold);
            highFlag = any(flags(j-1:j)<=0);
            break;
        end
    end

    if isLog
        low = exp(low);
        high = exp(high);
    end

end