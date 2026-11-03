%% OPTIMISATION PART 2 
% WITH TOXICITY HERE AND LETTING TIMING AND ANTI-PD-1 ALLOCATION VARY

clear;
clc;
close all;
format longG;

col = honoursPlotStyle();

% Tumour, immune and treatment parameters

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

% Persistent-blockade parameters calibrated before

p.muP = 587.0996052791231;
p.tauP = 60;

tox.lambda = log(2)/p.tauP;
tox.rho = log(2)/3;
tox.kI = tox.rho/2;
tox.kP = tox.rho/2;

settings.dt = 0.05;
settings.maxIterations = 100;
settings.maxPairedPasses = 12;
settings.improvementTolerance = 1e-7;
settings.feasibilityTolerance = 1e-10;
settings.batchSize = 512;

% IVM: 21 daily indicators. PD1: 22 daily amounts, half-units

reference = zeros(1,43);
reference(1:6) = 1;
reference(21+[1 8 15]) = 2;

[referencePeak, referencePeakTime] = toxicityPeak(reference,tox);
[referenceMetric, referenceTrajectory] = simulateSchedule(reference,p);
weeklyPD1 = reference;
weeklyPD1(1:21) = 0;
[weeklyPD1Metric, ~] = simulateSchedule(weeklyPD1,p);

fprintf('\nREFERENCE COMBINATION\n');
fprintf('AUC = %.6f; C(25) = %.6f\n',referenceMetric.AUC,referenceMetric.C25);
fprintf('Toxicity-surrogate peak = %.9f at day %.6f\n',referencePeak,referencePeakTime);

% IVM-day allowance, checkpoint-unit allowance and relative peak cap

caseNames = {'I6','P3','IP6','IP6_cap090','IP12','IP21'};
caseSettings = [6 0 1; 0 3 1; 6 3 1; 6 3 0.9; 12 3 1; 21 3 1];
results = struct();
summaryRows = zeros(8,11);
summaryRows(1,:) = [6 3 1 6 3 referenceMetric.AUC referenceMetric.C25 referenceMetric.B25 1 referencePeakTime 0];

for c = 1:numel(caseNames)
    name = caseNames{c};
    limits.BI = caseSettings(c,1);
    limits.BP = caseSettings(c,2);
    limits.cap = caseSettings(c,3);
    limits.referencePeak = referencePeak;

    fprintf('\n%s: IVM <= %g days, PD1 <= %g units, peak cap %.2f\n',name,limits.BI,limits.BP,limits.cap);
    starts = startingSchedules(name);
    best = optimiseSchedules(starts,limits,settings,p,tox);

    % Recalculate the selected trajectory using adaptive integration.

    [metric, trajectory] = simulateSchedule(best,p);
    [peak, peakTime] = toxicityPeak(best,tox);
    relativePeak = peak/referencePeak;
    assert(relativePeak<=limits.cap+settings.feasibilityTolerance,'The selected schedule exceeds its peak cap.');
    metric.Peak = peak;
    metric.PeakTime = peakTime;
    metric.RelativePeak = relativePeak;

    summaryRows(c+1,:) = [limits.BI limits.BP limits.cap sum(best(1:21)) sum(best(22:43))/2 metric.AUC metric.C25 metric.B25 relativePeak peakTime 100*(1-metric.AUC/referenceMetric.AUC)];
    results.(name).control = best;
    results.(name).metrics = metric;
    results.(name).trajectory = trajectory;

    doseIndex = find(best(22:43)>0);
    fprintf('IVM-active days: %s\n',mat2str(3+find(best(1:21)>0)));
    fprintf('PD1 days: %s; amounts: %s\n',mat2str(3+doseIndex),mat2str(best(21+doseIndex)/2));
    fprintf('AUC = %.6f; C(25) = %.6f; relative peak = %.6f\n',metric.AUC,metric.C25,relativePeak);
end

% Recomputing the 17-day schedule reported in the thesis, separately from the search.

reported = zeros(1,43);
reported([4:11 13:14 16:17 19:20 22:24]-3) = 1;
reported(22:23) = 2;
[metric, trajectory] = simulateSchedule(reported,p);
[peak, peakTime] = toxicityPeak(reported,tox);
metric.Peak = peak;
metric.PeakTime = peakTime;
metric.RelativePeak = peak/referencePeak;
assert(metric.RelativePeak<=1+settings.feasibilityTolerance,'The reported schedule exceeds its peak cap.');
results.IP21_Thesis.control = reported;
results.IP21_Thesis.metrics = metric;
results.IP21_Thesis.trajectory = trajectory;
summaryRows(8,:) = [21 3 1 sum(reported(1:21)) sum(reported(22:43))/2 metric.AUC metric.C25 metric.B25 metric.RelativePeak peakTime 100*(1-metric.AUC/referenceMetric.AUC)];

summary = array2table(summaryRows,'VariableNames',{'IVMAllowance','PD1Allowance','PeakCap','IVMDaysUsed','PD1UnitsUsed','AUC','C25','B25','RelativePeak','PeakTime','AUCImprovementPercent'});
summaryNames = ["Reference"; string(caseNames(:)); "IP21_Thesis"];
summaryNames(summaryNames=="IP21") = "IP21_Search";
summary = addvars(summary,summaryNames,'Before',1,'NewVariableNames','Case');

fprintf('\nFLEXIBLE SCHEDULES: BEST FOUND AND REPORTED IP21\n');
disp(summary);
fprintf('\nWeekly anti-PD-1: AUC = %.6f; C(25) = %.6f\n',weeklyPD1Metric.AUC,weeklyPD1Metric.C25);
fprintf('Flexible anti-PD-1 AUC improvement over weekly dosing = %.2f%%\n',100*(1-results.P3.metrics.AUC/weeklyPD1Metric.AUC));

%% FIGURE 1: Flexible treatment comparison

figure('Position',[100 100 950 500]);
plot(referenceTrajectory.Time,referenceTrajectory.A+referenceTrajectory.B,'--','Color',col.no,'LineWidth',2);
hold on;

chosen = {'I6','P3','IP6','IP6_cap090'};
colours = [col.ivm; col.pd1; col.combo; col.T];

for j = 1:numel(chosen)
    tr = results.(chosen{j}).trajectory;
    plot(tr.Time,tr.A+tr.B,'Color',colours(j,:),'LineWidth',2.2);
end

xlabel('Time (days)');
ylabel('Tumour volume (mm$^3$)');
legend('Reference combination','Flexible IVM only','Flexible anti-PD-1','Flexible combination, cap 1','Flexible combination, cap 0.9','Location','northwest');
xlim([4 25]);
grid on;

%% FIGURE 2: Selected schedules

figure('Position',[100 100 1000 570]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

% Use the reported IP21 candidate for Figure 4.

selected = {'IP6','IP12','IP21_Thesis'};
titles = {'6-day IVM allowance','12-day IVM allowance','21-day IVM allowance (reported)'};

for j = 1:3
    x = results.(selected{j}).control;
    nexttile;
    yyaxis left;
    stairs(4:25,[x(1:21) 0],'Color',col.ivm,'LineWidth',2);
    ylabel('IVM indicator');
    ylim([-0.05 1.15]);
    yyaxis right;
    stem(4:25,x(22:43)/2,'Color',col.pd1,'LineWidth',1.6);
    ylabel('PD1 increment');
    ylim([0 1.15]);
    title(titles{j});
    xlim([4 25]);
    grid on;
end

xlabel('Time (days)');

%% FUNCTIONS

function best = optimiseSchedules(starts,limits,settings,p,tox)

    best = [];
    bestAUC = Inf;

    for j = 1:size(starts,1)
        [x, value] = localSearch(starts(j,:),limits,settings,p,tox);
        if value<bestAUC
            best = x;
            bestAUC = value;
        end
        fprintf('Start %d/%d: AUC %.6f\n',j,size(starts,1),value);
    end

    assert(~isempty(best),'No feasible starting schedule found.');
    if limits.BI==0 || limits.BP==0
        return;
    end

    % Test simultaneous changes to IVM and checkpoint allocation.

    for pass = 1:settings.maxPairedPasses
        NI = neighbours(best,true,false);
        NP = neighbours(best,false,true);
        candidates = repelem(NI,size(NP,1),1);
        candidates(:,22:43) = repmat(NP(:,22:43),size(NI,1),1);
        values = scheduleCost(candidates,limits,settings,p,tox);
        [value, index] = min(values);

        if bestAUC-value<=settings.improvementTolerance
            return;
        end
        [best, bestAUC] = localSearch(candidates(index,:),limits,settings,p,tox);
    end

    warning('Paired-move iteration limit reached.');

end

function [x, value] = localSearch(x,limits,settings,p,tox)

    value = scheduleCost(x,limits,settings,p,tox);
    if ~isfinite(value)
        return;
    end

    for iteration = 1:settings.maxIterations
        candidates = neighbours(x,limits.BI>0,limits.BP>0);
        values = scheduleCost(candidates,limits,settings,p,tox);
        [newValue, index] = min(values);

        if newValue>=value-settings.improvementTolerance
            return;
        end
        x = candidates(index,:);
        value = newValue;
    end

    warning('Local-search iteration limit reached.');

end

function N = neighbours(x,includeI,includeP)

    % Add, remove or relocate one treatment increment.

    N = zeros(0,43);

    if includeI
        for i = 1:21
            y = x;
            y(i) = 1-y(i);
            N(end+1,:) = y;
        end
        for i = find(x(1:21)>0)
            for j = find(x(1:21)==0)
                y = x;
                y(i) = 0;
                y(j) = 1;
                N(end+1,:) = y;
            end
        end
    end

    if includeP
        for i = 22:43
            if x(i)<2
                y = x;
                y(i) = y(i)+1;
                N(end+1,:) = y;
            end
            if x(i)>0
                y = x;
                y(i) = y(i)-1;
                N(end+1,:) = y;
            end
        end
        for i = 21+find(x(22:43)>0)
            for j = 21+find(x(22:43)<2)
                if i==j
                    continue;
                end
                y = x;
                y(i) = y(i)-1;
                y(j) = y(j)+1;
                N(end+1,:) = y;
            end
        end
    end

end

function cost = scheduleCost(X,limits,settings,p,tox)

    valid = all(X(:,1:21)>=0 & X(:,1:21)<=1,2) & all(X(:,22:43)>=0 & X(:,22:43)<=2,2) & all(X==round(X),2);
    valid = valid & sum(X(:,1:21),2)<=limits.BI & sum(X(:,22:43),2)/2<=limits.BP;
    cost = inf(size(X,1),1);
    rows = find(valid);

    if isempty(rows)
        return;
    end

    peaks = toxicityPeak(X(rows,:),tox);
    rows = rows(peaks/limits.referencePeak<=limits.cap+settings.feasibilityTolerance);

    for j = 1:settings.batchSize:numel(rows)
        take = rows(j:min(j+settings.batchSize-1,numel(rows)));
        cost(take) = tumourAUC(X(take,:),settings.dt,p);
    end

end

function auc = tumourAUC(X,dt,p)

    % RK4 with exact checkpoint decay at each stage.

    N = size(X,1);
    nSteps = round(1/dt);
    h = 1/nSteps;
    decay = exp(-log(2)/p.tauP*h/2);
    Y = repmat([p.C0*(1-p.fB0),p.C0*p.fB0,p.s/p.d,0,0],N,1);
    P = zeros(N,1);
    invalid = false(N,1);

    for day = 1:21
        P = P+X(:,21+day)/2;
        u = X(:,day);

        for step = 1:nSteps
            k1 = tumourRHS(Y,P,u,p);
            k2 = tumourRHS(Y+h*k1/2,P*decay,u,p);
            k3 = tumourRHS(Y+h*k2/2,P*decay,u,p);
            k4 = tumourRHS(Y+h*k3,P*decay^2,u,p);
            Y = Y+h/6*(k1+2*k2+2*k3+k4);
            P = P*decay^2;
            invalid = invalid | any(~isfinite(Y),2) | any(Y(:,1:4)<-1e-6,2);
        end
    end

    auc = Y(:,5);
    auc(invalid) = Inf;

end

function F = tumourRHS(Y,P,u,p)

    A = Y(:,1);
    B = Y(:,2);
    T = Y(:,3);
    D = Y(:,4);
    C = A+B;
    GA = p.rA*A.*(1-C/p.KA);
    GB = p.rB*B.*(1-C/p.KB);
    m = 1+p.muP*P;

    dA = (1-p.eps)*GA+p.eps*GB-m*p.kA.*T.*A-p.deltaIVM*u.*A;
    dB = p.eps*GA+(1-p.eps)*GB-m*p.kB.*T.*B-p.deltaIVM*u.*B;
    dT = p.s-p.d*T+p.a*C.*T+p.pD*D./(p.hD+D).*T;
    dD = p.etaIVM*u.*C-p.dD*D;
    F = [dA dB dT dD C];

end

function [peak, when] = toxicityPeak(X,tox)

    % Include interval endpoints, interior maxima and the tail after day 25.

    Z = zeros(size(X,1),1);
    P = Z;
    peak = Z;
    when = 4+Z;

    for day = 1:22
        P = P+X(:,21+day)/2;
        b = tox.kP*P;
        c = zeros(size(Z));
        H = Inf;
        if day<=21
            c = tox.kI*X(:,day);
            H = 1;
        end

        h = inf(size(Z));
        rise = c+b-tox.rho*Z;
        good = b>0 & rise>0;
        scale = zeros(size(Z));
        scale(good) = rise(good)./(tox.lambda*b(good));
        ratio = (tox.rho-tox.lambda)*scale;
        good = good & ratio>-1;
        equal = good & ratio==0;
        h(equal) = scale(equal);
        other = good & ratio~=0;
        h(other) = scale(other).*log1p(ratio(other))./ratio(other);

        rows = find(isfinite(h) & h>0 & h<H);
        value = toxicityValue(h(rows),Z(rows),b(rows),c(rows),tox);
        better = value>peak(rows);
        rows = rows(better);
        peak(rows) = value(better);
        when(rows) = 4+(day-1)+h(rows);

        if day<=21
            Z = toxicityValue(1,Z,b,c,tox);
            P = P*exp(-tox.lambda);
            better = Z>peak;
            peak(better) = Z(better);
            when(better) = 4+day;
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
    Z = Z0.*exp(-tox.rho*h)+c/tox.rho.*(-expm1(-tox.rho*h))+b.*kernel;

end

function [metric, trajectory] = simulateSchedule(x,p)

    options = odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',0.1);
    y = [p.C0*(1-p.fB0);p.C0*p.fB0;p.s/p.d;0;0;0];
    time = (4:0.05:25)';
    Y = zeros(numel(time),6);

    for day = 1:21
        left = day+3;
        y(5) = y(5)+x(21+day)/2;
        sol = ode45(@(t,y) scheduleRHS(t,y,p,x(day)),[left left+1],y,options);
        assert(abs(sol.x(end)-left-1)<1e-8 && all(isfinite(sol.y(:))) && min(min(sol.y(1:5,:)))>=-1e-6,'Invalid numerical trajectory.');
        rows = (day-1)*20+(1:20);
        Y(rows,:) = deval(sol,time(rows))';
        y = sol.y(:,end);
    end

    y(5) = y(5)+x(43)/2;
    Y(end,:) = y';
    metric.AUC = y(6);
    metric.C25 = y(1)+y(2);
    metric.B25 = y(2);
    trajectory = array2table([time Y],'VariableNames',{'Time','A','B','T','D','P','Q'});

end

function F = scheduleRHS(~,y,p,u)

    f = tumourRHS([y(1:4)' y(6)],y(5),u,p);
    F = [f(1:4)';-log(2)/p.tauP*y(5);f(5)];

end

function X = startingSchedules(name)

    % Original starting schedules: IVM days and allowing half-unit checkpoint days.
    % Repeated checkpoint days represent a full unit.

    common = {[],[]; 4:9,[]; [],[4 4 11 11 18 18]; 4:9,[4 4 11 11 18 18]; 4:9,[5 5 12 12 19 19]};
    mixed6 = {[7 8 14 19 22],[10 16 20]; 24,[7 10 22 23 25]; [8 14 17 19],[7 10 11 19]; [14 18],[12 18 21]; [4 9 11 14 18],[6 8 21]};

    switch name
        case 'I6'
            seeds = {[],[]; 4:9,[]; 7:12,[]; 19:24,[]; [7 8 14 19 22],[]; [10 16 18],[]; [7 10 18 20 23],[]; 10,[]; [19 22],[]; 21,[]};
        case 'P3'
            seeds = {[],[]; [],[4 4 11 11 18 18]; [],[5 5 12 12 19 19]; [],[4 8 12 16 20 24]; [],[7 9 16 23 23]; [],[11 11 20 20 24]; [],[7 16]; [],[7 10 19 22 23 25]; [],[4 8 9 22]; [],[10 21]};
        case 'IP6'
            seeds = [common; mixed6];
        case 'IP6_cap090'
            seeds = [common(1:2,:); mixed6; {24,[20 22]; 17,21; [9 10 17 21],[8 11 18 19]}];
        case 'IP12'
            seeds = [common; {4:15,[]; [9 18 22 24],[17 23]; [7 8 11 13 14 15 19 21 23],[4 13 16 18 19]; [4 5 7 8 10 12 13 14 19],[6 6 9 17 22 23]; [20 24],[12 17 18 25 25]}];
        case 'IP21'
            seeds = [common; {4:24,[]; [5 6 7 8 10 11 12 13 14 15 16 17 18 19 21 22 23],[8 17 20 23]; [4 5 6 7 8 9 10 11 13 14 15 16 17 18 19 21 22 23],[6 10 15 22]; [4 5 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 24],[7 15 18 19]; [5 6 7 8 9 10 11 13 15 17 18 19 20 22 23],[11 18]}];
    end

    X = zeros(size(seeds,1),43);
    for j = 1:size(seeds,1)
        X(j,seeds{j,1}-3) = 1;
        for day = seeds{j,2}
            X(j,day+18) = X(j,day+18)+1;
        end
    end

end