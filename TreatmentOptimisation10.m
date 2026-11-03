%% TREATMENT OPTIMISATION

clear;
clc;
close all;
format longG;

col = honoursPlotStyle();

% Anti-PD-1 data

data = readtable('Data_PD1_Actual.csv');
t_data = [4; 7; 11; 14; 18; 21; 25];
C_data = data{:,2};
C_data = C_data(:);

% Tumour, immune and IVM parameters

p.rA = 0.365371;
p.rB = 0.365371;
p.KA = 58.772570;
p.KB = 413.826598;
p.eps = 0.05;
p.kA = 1e-5;
p.kB = 1e-6;
p.s = 10;
p.d = 0.2;
p.a = 1e-4;
p.deltaIVM = 0.132649;
p.etaIVM = 0.08;
p.dD = 0.25;
p.pD = 0.05;
p.hD = 10;
p.C0 = 10;
p.fB0 = 0.15;

ode_options = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.1);
fit_options = optimset('Display','off','TolX',1e-8,'TolFun',1e-8,'MaxIter',1500,'MaxFunEvals',5000);

% Fit checkpoint potency and half-life using three weekly doses.

lower = [1e-3; 0.25];
upper = [1e6; 60];
mu_grid = logspace(-1,5,26);
tau_grid = logspace(log10(0.5),log10(30),26);
SSE_grid = zeros(numel(tau_grid),numel(mu_grid));
objective = @(z) blockadeSSE(z,p,t_data,C_data,ode_options,lower,upper);

fprintf('\nCALIBRATING CHECKPOINT BLOCKADE\n');
for i = 1:numel(tau_grid)
    for j = 1:numel(mu_grid)
        SSE_grid(i,j) = objective(log([mu_grid(j); tau_grid(i)]));
    end
end

[~, order] = sort(SSE_grid(:));
fit_rows = zeros(7,4);

for k = 1:5
    [i,j] = ind2sub(size(SSE_grid),order(k));
    z0 = log([mu_grid(j); tau_grid(i)]);
    [z, cost, flag] = fminsearch(objective,z0,fit_options);
    fit_rows(k,:) = [exp(z(:))', cost, flag];
end

% Include both half-life bounds in the constrained fit.

scalar_options = optimset('Display','off','TolX',1e-9);
tau_bounds = [lower(2), upper(2)];

for k = 1:2
    tau = tau_bounds(k);
    [z, cost, flag] = fminbnd(@(z) objective([z; log(tau)]),log(lower(1)),log(upper(1)),scalar_options);
    fit_rows(k+5,:) = [exp(z), tau, cost, flag];
end

[fit_cost, best] = min(fit_rows(:,3));
assert(isfinite(fit_cost) && fit_cost<1e90,'Checkpoint calibration failed.');
p.muP = fit_rows(best,1);
p.tauP = fit_rows(best,2);
fitTable = array2table(fit_rows,'VariableNames',{'muP','tauP','SSE','ExitFlag'});

fit = simulateTreatment(p,NaN,4,t_data,ode_options);
residuals = fit.C-C_data;
SSE = sum(residuals.^2);
RMSE = sqrt(mean(residuals.^2));
RMSE_after_IC = sqrt(mean(residuals(2:end).^2));
calibrationTable = table(t_data,C_data,fit.C,residuals,'VariableNames',{'Day','Observed','Predicted','Residual'});

fprintf('muP = %.9f; tauP = %.9f days\n',p.muP,p.tauP);
fprintf('SSE = %.9f; RMSE = %.9f; RMSE days 7-25 = %.9f\n',SSE,RMSE,RMSE_after_IC);

% Six IVM days and three weekly checkpoint doses; only start times vary.

fine_IVM = (4:0.125:19)';
fine_PD1 = (4:0.125:11)';
fine_AUC = zeros(numel(fine_IVM),numel(fine_PD1));
fine_C25 = fine_AUC;
fine_Bpercent = fine_AUC;

fprintf('\nSEARCHING TREATMENT START TIMES\n');
for i = 1:numel(fine_IVM)
    for j = 1:numel(fine_PD1)
        sim = simulateTreatment(p,fine_IVM(i),fine_PD1(j),[4;25],ode_options);
        fine_AUC(i,j) = sim.AUC;
        fine_C25(i,j) = sim.C(end);
        fine_Bpercent(i,j) = 100*sim.B(end)/sim.C(end);
    end
    if mod(i,10)==0 || i==numel(fine_IVM)
        fprintf('Completed IVM start %d / %d\n',i,numel(fine_IVM));
    end
end

% Every second point of the fine grid gives the original 0.25-day grid.

steps = [0.25 0.125];
grid_rows = zeros(2,13);

for k = 1:2
    stride = round(steps(k)/0.125);
    tI = fine_IVM(1:stride:end);
    tP = fine_PD1(1:stride:end);
    AUC = fine_AUC(1:stride:end,1:stride:end);
    C25 = fine_C25(1:stride:end,1:stride:end);

    [minimum_AUC, index] = min(AUC(:));
    [ia,ja] = ind2sub(size(AUC),index);
    [minimum_C25, index] = min(C25(:));
    [ie,je] = ind2sub(size(C25),index);
    [ni,np] = find(AUC<=1.01*minimum_AUC);

    grid_rows(k,:) = [steps(k),tI(ia),tP(ja),minimum_AUC,C25(ia,ja),tI(ie),tP(je),minimum_C25,numel(ni),min(tI(ni)),max(tI(ni)),min(tP(np)),max(tP(np))];

    if k==1
        t_IVM = tI;
        t_PD1 = tP;
        AUC_grid = AUC;
        best_AUC_times = [tI(ia),tP(ja)];
        best_endpoint_times = [tI(ie),tP(je)];
    end
end

grid_names = {'Step','AUC_IVM','AUC_PD1','MinimumAUC','C25_at_AUC','Endpoint_IVM','Endpoint_PD1','MinimumC25','NearCount','NearIVM_min','NearIVM_max','NearPD1_min','NearPD1_max'};
gridSummary = array2table(grid_rows,'VariableNames',grid_names);
nearTable = table(fine_IVM(ni),fine_PD1(np),'VariableNames',{'IVM_Start','PD1_FirstDose'});
disp(gridSummary);

% Toxicity surrogate, normalised by the reference combination.

tox.rho = log(2)/3;
tox.lambda = log(2)/p.tauP;
tox.kI = tox.rho/2;
tox.kP = tox.rho/2;
tox.Pref = 1;
[reference_peak, reference_peak_day] = toxicityPeak(4,4,tox);

% Independently optimise each fixed-course monotherapy.

accurate_options = odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',0.1);
start_grids = {fine_IVM,fine_PD1};
arm_names = ["IVM","Anti-PD-1"];
arm = strings(numel(fine_IVM)+numel(fine_PD1),1);
mono_rows = zeros(numel(arm),7);
row = 0;

for g = 1:2
    for start = start_grids{g}'
        tI = NaN;
        tP = NaN;
        if g==1
            tI = start;
        else
            tP = start;
        end

        sim = simulateTreatment(p,tI,tP,[4;25],accurate_options);
        [peak, peak_day] = toxicityPeak(tI,tP,tox);
        row = row+1;
        arm(row) = arm_names(g);
        mono_rows(row,:) = [start,sim.AUC,sim.C(end),sim.B(end),100*sim.B(end)/sim.C(end),peak/reference_peak,peak_day];
    end
end

mono_names = {'Start','AUC','C25','B25','Bpercent','RelativePeak','PeakDay'};
monoGrid = [table(arm,'VariableNames',{'Arm'}),array2table(mono_rows,'VariableNames',mono_names)];
best_rows = zeros(2,1);

for g = 1:2
    rows = find(monoGrid.Arm==arm_names(g));
    [~, k] = min(monoGrid.AUC(rows));
    best_rows(g) = rows(k);
end

monoBest = monoGrid(best_rows,:);
fprintf('\nOPTIMISED MONOTHERAPIES\n');
disp(monoBest);

% Compare untreated, reference, combination optima and monotherapies.

names = {'Untreated';'Reference combination';'AUC optimum';'Endpoint optimum';'Optimised IVM';'Optimised anti-PD-1'};
schedules = [NaN NaN;4 4;best_AUC_times;best_endpoint_times;monoBest.Start(1) NaN;NaN monoBest.Start(2)];
t_plot = unique([(4:0.05:25)';t_data]);
simulations = cell(6,1);
metrics = zeros(6,10);

for k = 1:6
    sim = simulateTreatment(p,schedules(k,1),schedules(k,2),t_plot,accurate_options);
    [peak, peak_day] = toxicityPeak(schedules(k,1),schedules(k,2),tox);
    simulations{k} = sim;
    metrics(k,:) = [sim.AUC,sim.C(end),sim.A(end),sim.B(end),100*sim.B(end)/sim.C(end),6*isfinite(schedules(k,1)),3*isfinite(schedules(k,2)),peak,peak_day,peak/reference_peak];
end

metric_names = {'IVM_Start','PD1_Start','AUC','C25','A25','B25','Bpercent','IVMdays','PD1units','Peak','PeakDay','RelativePeak'};
summary = [table(names,'VariableNames',{'Schedule'}),array2table([schedules,metrics],'VariableNames',metric_names)];

fprintf('\nTREATMENT COMPARISON\n');
disp(summary);
fprintf('AUC reduction from reference = %.6f %%\n',100*(1-metrics(3,1)/metrics(2,1)));
fprintf('Endpoint reduction from AUC optimum = %.6f %%\n',100*(1-metrics(4,2)/metrics(3,2)));
fprintf('AUC increase at endpoint optimum = %.6f %%\n',100*(metrics(4,1)/metrics(3,1)-1));
fprintf('Fine-grid schedules within 1%% of minimum = %d\n',height(nearTable));
fprintf('Near-optimal starts: IVM %.3f-%.3f; anti-PD-1 %.3f-%.3f days\n',min(nearTable.IVM_Start),max(nearTable.IVM_Start),min(nearTable.PD1_FirstDose),max(nearTable.PD1_FirstDose));

%% FIGURE 1: Treatment timing

figure('Position',[100 100 850 600]);
imagesc(t_PD1,t_IVM,100*(AUC_grid/metrics(2,1)-1));
axis xy;
colorbar;
hold on;
plot(4,4,'o','Color',col.data,'MarkerFaceColor','w','MarkerSize',9,'LineWidth',1.5);
plot(best_AUC_times(2),best_AUC_times(1),'p','Color',col.data,'MarkerFaceColor',col.data,'MarkerSize',13);
contour(t_PD1,t_IVM,AUC_grid,[1.01*min(AUC_grid(:)),1.01*min(AUC_grid(:))],'w','LineWidth',2);
xlabel('First anti-PD-1 dose (day)');
ylabel('IVM course start (day)');
title('Tumour AUC change relative to reference (%)');
grid on;

%% FIGURE 2: AUC and endpoint objectives

figure('Position',[100 100 850 520]);
colours = [col.no;col.combo;col.T;col.pd1];
styles = {'--','-','-',':'};

for k = 1:4
    sim = simulations{k};
    plot(sim.t,sim.C,styles{k},'Color',colours(k,:),'LineWidth',2.5);
    hold on;
end

xlabel('Time (days)');
ylabel('Tumour volume (mm$^3$)');
legend(names(1:4),'Location','northwest');
xlim([4 25]);
grid on;

%% FUNCTIONS

function SSE = blockadeSSE(z,p,t_data,C_data,options,lower,upper)

    v = exp(z(:));
    SSE = 1e100;
    if any(~isfinite(v)) || any(v<lower*(1-1e-12)) || any(v>upper*(1+1e-12))
        return;
    end

    p.muP = v(1);
    p.tauP = v(2);

    try
        sim = simulateTreatment(p,NaN,4,t_data,options);
        SSE = sum((sim.C-C_data).^2);
        if ~isfinite(SSE)
            SSE = 1e100;
        end
    catch
        SSE = 1e100;
    end

end

function sim = simulateTreatment(p,tI,tP,t_eval,options)

    % Integrate between treatment events; the sixth state accumulates AUC.

    doses = tP+(0:2)'*7;
    doses = doses(isfinite(doses));
    events = unique([4;25;tI;tI+6;doses]);
    events = events(isfinite(events));
    t_eval = t_eval(:);
    y = [p.C0*(1-p.fB0);p.C0*p.fB0;p.s/p.d;0;0;0];
    Y = nan(numel(t_eval),6);

    for k = 1:numel(events)-1
        left = events(k);
        right = events(k+1);
        y(5) = y(5)+nnz(doses==left);
        uIVM = double(isfinite(tI) && (left+right)/2>=tI && (left+right)/2<tI+6);
        sol = ode45(@(t,y) treatmentRHS(y,p,uIVM),[left right],y,options);
        assert(abs(sol.x(end)-right)<1e-8 && all(isfinite(sol.y(:))) && min(min(sol.y(1:5,:)))>=-1e-6,'Invalid ODE solution.');

        take = t_eval>=left & t_eval<right;
        if any(take)
            Y(take,:) = deval(sol,t_eval(take))';
        end
        y = sol.y(:,end);
    end

    y(5) = y(5)+nnz(doses==25);
    Y(t_eval==25,:) = repmat(y',nnz(t_eval==25),1);
    sim.t = t_eval;
    sim.A = Y(:,1);
    sim.B = Y(:,2);
    sim.T = Y(:,3);
    sim.D = Y(:,4);
    sim.P = Y(:,5);
    sim.C = sim.A+sim.B;
    sim.AUC = y(6);

end

function dydt = treatmentRHS(y,p,uIVM)

    A = y(1);
    B = y(2);
    T = y(3);
    D = y(4);
    P = y(5);
    C = A+B;
    GA = p.rA*A*(1-C/p.KA);
    GB = p.rB*B*(1-C/p.KB);
    m = 1+p.muP*P;

    dA = (1-p.eps)*GA+p.eps*GB-m*p.kA*T*A-p.deltaIVM*uIVM*A;
    dB = p.eps*GA+(1-p.eps)*GB-m*p.kB*T*B-p.deltaIVM*uIVM*B;
    dT = p.s-p.d*T+p.a*C*T+p.pD*D/(p.hD+D)*T;
    dD = p.etaIVM*uIVM*C-p.dD*D;
    dP = -log(2)/p.tauP*P;
    dydt = [dA;dB;dT;dD;dP;C];

end

function [peak, peak_day] = toxicityPeak(tI,tP,tox)

    % Check interval endpoints and stationary points, including after day 25.

    doses = tP+(0:2)'*7;
    doses = doses(isfinite(doses));
    events = unique([4;25;tI;tI+6;doses]);
    events = events(isfinite(events));
    Z = 0;
    P = 0;
    peak = 0;
    peak_day = 4;

    for k = 1:numel(events)
        left = events(k);
        P = P+nnz(doses==left);
        L = Inf;
        if k<numel(events)
            L = events(k+1)-left;
        end

        uIVM = double(isfinite(tI) && left>=tI && left<tI+6);
        b = tox.kP*P/tox.Pref;
        c = tox.kI*uIVM;
        times = 0;
        if isfinite(L)
            times(end+1) = L;
        end

        derivative = c+b-tox.rho*Z;
        if b>0 && derivative>0
            scale = derivative/(tox.lambda*b);
            ratio = (tox.rho-tox.lambda)*scale;
            if ratio>-1
                h = scale;
                if ratio~=0
                    h = scale*log1p(ratio)/ratio;
                end
                if h>0 && h<L
                    times(end+1) = h;
                end
            end
        end

        [value, index] = max(toxicityValue(times,Z,b,c,tox));
        if value>peak
            peak = value;
            peak_day = left+times(index);
        end

        if isfinite(L)
            Z = toxicityValue(L,Z,b,c,tox);
            P = P*exp(-tox.lambda*L);
        end
    end

end

function Z = toxicityValue(h,Z0,b,c,tox)

    gap = abs(tox.rho-tox.lambda);
    if gap==0
        kernel = h.*exp(-tox.rho*h);
    else
        kernel = exp(-min(tox.rho,tox.lambda)*h).*(-expm1(-gap*h))/gap;
    end

    Z = Z0*exp(-tox.rho*h)+c/tox.rho*(-expm1(-tox.rho*h))+b*kernel;

end