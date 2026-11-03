%% STABILITY ANALYSIS

clear;
clc;
close all;
format longG;

col = honoursPlotStyle();

% Homogeneous parameters

r_C = 0.21372030053324256;
K_C = 556.4803820509597;
k_C = 1e-5;
s = 10;
d = 0.2;
a = 1e-4;
Tb = s/d;

kcrit = r_C*d/s;
b = r_C/K_C + a*r_C/d;

% Tumour-free and admissible coexistence equilibria

caseNames = {'Reference', '0.5 kcrit', '1.5 kcrit', 'Threshold'};
kValues = [k_C, 0.5*kcrit, 1.5*kcrit, kcrit];
homogeneousRows = cell(0,9);

for j = 1:numel(kValues)
    k = kValues(j);
    states = [0; Tb];
    equilibriumNames = {'Tumour-free'};

    if r_C-k*Tb > 1e-12*max([r_C, k*Tb, d])
        discriminant = (K_C-d/a)^2 + 4*K_C*k*s/(r_C*a);
        Cstar = (2*K_C/a)*(d-k*s/r_C)/(K_C+d/a+sqrt(discriminant));
        Tstar = s/(d-a*Cstar);
        states(:,end+1) = [Cstar; Tstar];
        equilibriumNames{end+1} = 'Coexistence';
    end

    for n = 1:size(states,2)
        Cstar = states(1,n);
        Tstar = states(2,n);
        J = [r_C*(1-2*Cstar/K_C)-k*Tstar, -k*Cstar; a*Tstar, -d+a*Cstar];
        eigenvalues = eig(J);
        [~, order] = sortrows([real(eigenvalues), imag(eigenvalues)], [-1 -2]);
        eigenvalues = eigenvalues(order);
        tolerance = 1e-10*max([r_C, d, k*Tb, abs(trace(J)), sqrt(abs(det(J))), 1e-12]);
        realParts = real(eigenvalues);

        if any(realParts > tolerance) && any(realParts < -tolerance)
            classification = 'Saddle';
        elseif all(realParts < -tolerance)
            if any(abs(imag(eigenvalues)) > tolerance)
                classification = 'Stable focus';
            else
                classification = 'Stable node';
            end
        elseif any(realParts > tolerance)
            classification = 'Unstable';
        else
            classification = 'Nonhyperbolic';
        end

        homogeneousRows(end+1,:) = {caseNames{j}, k, k/kcrit, equilibriumNames{n}, Cstar, Tstar, eigenvalues(1), eigenvalues(2), classification};
    end
end

homogeneousTable = cell2table(homogeneousRows, 'VariableNames', {'Case','k_C','k_over_kcrit','Equilibrium','Cstar','Tstar','Lambda1','Lambda2','Classification'});

fprintf('\nHOMOGENEOUS STABILITY\n');
fprintf('kcrit = %.12g\n', kcrit);
fprintf('b = %.12g\n', b);
disp(homogeneousTable);

%% FIGURE 1: Homogeneous bifurcation

ratio = unique([linspace(0.001,1,500), linspace(1,1.6,300), k_C/kcrit]);
k_branch = ratio*kcrit;
discriminant = (K_C-d/a)^2 + 4*K_C*k_branch*s/(r_C*a);
C_branch = (2*K_C/a)*(d-k_branch*s/r_C)./(K_C+d/a+sqrt(discriminant));
positive = ratio <= 1;

figure('Position', [100 100 800 510]);
hold on;

plot(ratio(positive), C_branch(positive)/K_C, 'Color', col.total, 'LineWidth', 3, 'DisplayName', 'Stable coexistence');
plot([0 1], [0 0], '--', 'Color', col.B, 'LineWidth', 2, 'DisplayName', 'Unstable tumour-free');
plot([1 1.6], [0 0], 'Color', col.T, 'LineWidth', 3, 'DisplayName', 'Stable tumour-free');
plot(ratio(~positive), C_branch(~positive)/K_C, ':', 'Color', col.no, 'LineWidth', 2, 'DisplayName', 'Nonphysical continuation');
plot(1, 0, 'o', 'Color', col.data, 'MarkerFaceColor', 'w', 'MarkerSize', 7, 'HandleVisibility', 'off');

xline(k_C/kcrit, ':', 'Reference', 'Color', col.data, 'HandleVisibility', 'off');

xlabel('Killing ratio $k_C/k_{\mathrm{crit}}$');
ylabel('Equilibrium burden $C_*/K_C$');

xlim([0 1.6]);
ylim([min(C_branch/K_C)-0.04, 1.05]);

legend('Location', 'northeast');

%% HETEROGENEOUS STABILITY AND TREATMENT THRESHOLDS

% Constant inputs represent sustained treatment.

p.rA = 0.365371;
p.rB = 0.365371;
p.KA = 58.772570;
p.KB = 413.826598;
p.eps = 0.05;
p.kA = 1e-5;
p.kB = 1e-6;
p.s = s;
p.d = d;
p.a = a;
p.delta = 0.132649;
p.eta = 0.08;
p.dD = 0.25;
p.pD = 0.05;
p.hD = 10;
p.mu = 1048.95;

m_cal = 1+p.mu;
scale = [p.KB; p.KB; p.s/p.d; p.hD];

% Critical multipliers and local transcritical coefficients

mcrit = zeros(2,1);
q_limit = zeros(2,1);
a_H = zeros(2,1);
b_H = zeros(2,1);

mixing = [1-p.eps, p.eps; p.eps, 1-p.eps];
killing = diag([p.kA, p.kB]);

for j = 1:2
    uI = j-1;
    mcrit(j) = fzero(@(m) tumourFree(m,uI,p),[1 8000],optimset('TolX',1e-10,'Display','off'));
    [~, B] = tumourFree(mcrit(j),uI,p);

    [V, L] = eig(B);
    [~, index] = min(abs(diag(L)));
    v = real(V(:,index));
    v = v/sum(v);

    [W, L] = eig(B');
    [~, index] = min(abs(diag(L)));
    w = real(W(:,index));
    w = w/(w'*v);

    d1 = p.eta*uI/p.dD;
    t1 = Tb*(p.a+p.pD*d1/p.hD)/p.d;

    a_H(j) = Tb*(w'*killing*v);
    b_H(j) = w'*(mixing*diag([p.rA/p.KA, p.rB/p.KB])+mcrit(j)*t1*killing)*v;
    q_limit(j) = v(2);
end

slope_at_threshold = -a_H./b_H;
thresholds = table([0;1],mcrit,100*q_limit,a_H,b_H,slope_at_threshold,'VariableNames',{'uI','mcrit','LimitingPercentB','aH','bH','dC_dm'});

% Continue the branches and enumerate all admissible roots at each m.

base_grid = unique([1:2:8000, 4000, 8000, m_cal]);
branch_rows = zeros(0,9);
root_counts = zeros(0,4);

for j = 1:2
    uI = j-1;
    m_grid = unique([base_grid, mcrit(j), mcrit(j)-[100 10 1 0.1 0.01 0.001]]);
    previous = [];

    for m = m_grid
        states = positiveEquilibria(m,uI,p);

        if ~isempty(previous)
            [~, J] = heterogeneousModel(previous,previous_m,uI,p);
            Jscaled = (J.*scale')./(p.rB*scale);
            fm = [-p.kA*previous(3)*previous(1); -p.kB*previous(3)*previous(2); 0; 0];

            if rcond(Jscaled)>1e-14
                tangent = -Jscaled\(fm./(p.rB*scale));
                predicted = previous+(m-previous_m)*scale.*tangent;
                [tracked, ok] = correctEquilibrium(predicted,m,uI,p);

                if ok && sum(tracked(1:2))>1e-10 && size(states,2)==1
                    assert(norm((tracked-states)./scale,inf)<1e-7,'Continuation and root enumeration disagree.');
                end
            end
        end

        root_counts(end+1,:) = [uI, m, size(states,2), tumourFree(m,uI,p)];

        for k = 1:size(states,2)
            y = states(:,k);
            [~, J] = heterogeneousModel(y,m,uI,p);
            eigenvalues = eig((J.*scale')./scale);
            alpha = max(real(eigenvalues));
            Cstar = sum(y(1:2));

            branch_rows(end+1,:) = [uI, m, Cstar, y(2)/Cstar, y', alpha];
        end

        previous = [];

        if size(states,2)==1
            previous = states;
            previous_m = m;
        end
    end
end

branches = array2table(branch_rows,'VariableNames',{'uI','m','Cstar','qstar','Astar','Bstar','Tstar','Dstar','AlphaStar'});

% Four reference cases from the same equilibrium calculation

names = {'Untreated'; 'Sustained IVM'; 'Sustained anti-PD-1'; 'Sustained combination'};
u_cases = [0; 1; 0; 1];
m_cases = [1; 1; m_cal; m_cal];
reference_rows = zeros(4,4);

for j = 1:4
    index = find(branches.uI==u_cases(j) & branches.m==m_cases(j));
    assert(numel(index)==1,'Expected one admissible reference equilibrium.');

    reference_rows(j,:) = [branches.Cstar(index), 100*branches.qstar(index), tumourFree(m_cases(j),u_cases(j),p), branches.AlphaStar(index)];
end

reference_table = table(names,reference_rows(:,1),reference_rows(:,2),reference_rows(:,3),reference_rows(:,4),'VariableNames',{'Condition','Cstar','PercentB','LambdaPlus','AlphaStar'});

fprintf('\nSUSTAINED-TREATMENT EQUILIBRIA\n');
disp(reference_table);

fprintf('\nTREATMENT THRESHOLDS\n');
fprintf('Calibrated multiplier = %.6f\n',m_cal);
disp(thresholds);

fprintf('\nEQUILIBRIUM BURDENS AT m = 4000\n');
disp(branches(branches.m==4000,{'uI','m','Cstar'}));

% These summaries concern the sampled branches, not a global uniqueness proof.

for j = 1:2
    counts = root_counts(root_counts(:,1)==j-1,:);
    path = branch_rows(branch_rows(:,1)==j-1,:);
    below = counts(:,2)<mcrit(j);

    one_root = all(counts(below,3)==1);
    no_root_above = all(counts(~below,3)==0);
    stable = all(path(:,9)<0);
    decreasing = one_root && all(diff(path(:,3))<0);

    fprintf('\nuI = %d: one positive root below threshold = %d; none at/above = %d\n',j-1,one_root,no_root_above);
    fprintf('All sampled positive equilibria stable = %d; burden decreasing = %d\n',stable,decreasing);
end

%% FIGURE 2: Tumour-free stability boundary

m_boundary = linspace(1,8000,1200);
boundary = arrayfun(@(m) tumourFree(m,0,p),m_boundary);

figure('Position',[100 100 800 510]);
hold on;

fill([m_boundary fliplr(m_boundary)],[max(boundary,0) 0.4*ones(size(boundary))],col.T,'FaceAlpha',0.13,'EdgeColor','none','HandleVisibility','off');

keep = boundary>=0;

plot(m_boundary(keep),boundary(keep),'Color',col.total,'LineWidth',3,'DisplayName','$\lambda_+=0$');
plot([1 1 m_cal m_cal],[0 p.delta 0 p.delta],'o','Color',col.data,'MarkerFaceColor','w','MarkerSize',7,'DisplayName','Retained input cases');

text(4500,0.31,'Tumour-free stable');
text(2200,0.06,'Tumour-free unstable');

xlabel('Checkpoint killing multiplier $m$');
ylabel('Direct IVM loss $\ell$ (day$^{-1}$)');

xlim([1 8000]);
ylim([0 0.4]);

legend('Location','northeast');

%% FIGURE 3: Equilibrium tumour burden

figure('Position',[100 100 1000 460]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

titles = {'No IVM','Sustained IVM'};
colours = [col.pd1; col.combo];

for j = 1:2
    nexttile;
    hold on;

    path = branch_rows(branch_rows(:,1)==j-1,:);
    stable_burden = path(:,3);
    unstable_burden = path(:,3);

    stable_burden(path(:,9)>=0) = NaN;
    unstable_burden(path(:,9)<0) = NaN;

    hs = plot(path(:,2),stable_burden,'Color',colours(j,:),'LineWidth',3);
    plot(path(:,2),unstable_burden,'--','Color',colours(j,:),'LineWidth',3,'HandleVisibility','off');

    hu = plot([1 mcrit(j)],[0 0],'--','Color',colours(j,:),'LineWidth',1.8);
    plot([mcrit(j) 8000],[0 0],'Color',colours(j,:),'LineWidth',3,'HandleVisibility','off');
    plot(mcrit(j),0,'o','Color',colours(j,:),'MarkerFaceColor','w','HandleVisibility','off');

    xline(m_cal,':','Calibrated $m$','Color',col.data,'HandleVisibility','off');

    xlabel('Checkpoint multiplier $m$');
    ylabel('Equilibrium burden $C_*$ (mm$^3$)');
    title(titles{j});

    xlim([1 8000]);
    ylim([-10 440]);

    legend([hs hu],{'Stable','Unstable'},'Location','northeast');
end

%% FUNCTIONS

function [lambda, B] = tumourFree(m,uI,p)

    Tb = p.s/p.d;
    B = [(1-p.eps)*p.rA-m*p.kA*Tb-p.delta*uI, p.eps*p.rB; p.eps*p.rA, (1-p.eps)*p.rB-m*p.kB*Tb-p.delta*uI];
    lambda = 0.5*(trace(B)+sqrt((B(1,1)-B(2,2))^2+4*B(1,2)*B(2,1)));

end

function states = positiveEquilibria(m,uI,p)

    % det H(C)=0 gives a polynomial in x=C/KB, using ascending powers.

    z = p.eta*uI*p.KB/(p.dD*p.hD);
    N = [1 z];
    Q = polyAdd(conv([p.d,-p.a*p.KB],N),[0,-p.pD*z]);

    gA = [p.rA,-p.rA*p.KB/p.KA];
    gB = [p.rB,-p.rB];

    PA = polyAdd(conv(polyAdd((1-p.eps)*gA,-p.delta*uI),Q),-m*p.kA*p.s*N);
    PB = polyAdd(conv(polyAdd((1-p.eps)*gB,-p.delta*uI),Q),-m*p.kB*p.s*N);
    P = polyAdd(conv(PA,PB),-p.eps^2*conv(conv(gA,gB),conv(Q,Q)));

    while numel(P)>1 && P(end)==0
        P(end) = [];
    end

    candidates = roots(fliplr(P/max(abs(P))));
    x = sort(real(candidates(abs(imag(candidates))<1e-7.*max(1,abs(candidates)))));
    x = x(x>1e-12 & x<1+1e-9);

    states = zeros(4,0);
    scale = [p.KB; p.KB; p.s/p.d; p.hD];

    for j = 1:numel(x)
        C = p.KB*x(j);
        D = p.eta*uI*C/p.dD;
        omega = p.d-p.a*C-p.pD*D/(p.hD+D);

        if omega<=0
            continue;
        end

        T = p.s/omega;
        H11 = (1-p.eps)*p.rA*(1-C/p.KA)-m*p.kA*T-p.delta*uI;
        H12 = p.eps*p.rB*(1-C/p.KB);
        v = [H12; -H11];

        if abs(sum(v))<1e-12
            continue;
        end

        y = [C*v/sum(v); T; D];

        if ~admissible(y,p)
            continue;
        end

        [y, ok] = correctEquilibrium(y,m,uI,p);
        assert(ok,'An equilibrium could not be refined.');

        if isempty(states) || all(max(abs((states-y)./scale),[],1)>=1e-8)
            states(:,end+1) = y;
        end
    end

end

function [y, ok] = correctEquilibrium(y,m,uI,p)

    % Damped Newton iteration in scaled state variables.

    scale = [p.KB; p.KB; p.s/p.d; p.hD];
    ok = false;

    if ~admissible(y,p)
        return;
    end

    z = y./scale;

    for iteration = 1:15
        [f, J] = heterogeneousModel(z.*scale,m,uI,p);
        f = f./(p.rB*scale);

        if norm(f,inf)<1e-13
            break;
        end

        J = (J.*scale')./(p.rB*scale);

        if rcond(J)<1e-14
            break;
        end

        step = -J\f;
        accepted = false;

        for k = 0:18
            trial = z+2^(-k)*step;
            state = trial.*scale;

            if admissible(state,p) && norm(heterogeneousModel(state,m,uI,p)./(p.rB*scale),inf)<norm(f,inf)
                z = trial;
                accepted = true;
                break;
            end
        end

        if ~accepted
            break;
        end
    end

    y = z.*scale;
    ok = admissible(y,p) && norm(heterogeneousModel(y,m,uI,p)./(p.rB*scale),inf)<1e-11;

end

function [f, J] = heterogeneousModel(y,m,uI,p)

    A = y(1);
    B = y(2);
    T = y(3);
    D = y(4);
    C = A+B;

    GA = p.rA*A*(1-C/p.KA);
    GB = p.rB*B*(1-C/p.KB);

    dA = (1-p.eps)*GA+p.eps*GB-m*p.kA*T*A-p.delta*uI*A;
    dB = (1-p.eps)*GB+p.eps*GA-m*p.kB*T*B-p.delta*uI*B;
    dT = p.s-p.d*T+p.a*C*T+p.pD*D/(p.hD+D)*T;
    dD = p.eta*uI*C-p.dD*D;

    f = [dA; dB; dT; dD];

    M = [1-p.eps, p.eps; p.eps, 1-p.eps];
    RG = [p.rA*(1-(2*A+B)/p.KA), -p.rA*A/p.KA; -p.rB*B/p.KB, p.rB*(1-(A+2*B)/p.KB)];

    J = zeros(4);
    J(1:2,1:2) = M*RG-diag([m*p.kA*T+p.delta*uI, m*p.kB*T+p.delta*uI]);
    J(1:2,3) = -m*[p.kA*A; p.kB*B];
    J(3,:) = [p.a*T, p.a*T, -p.d+p.a*C+p.pD*D/(p.hD+D), p.pD*p.hD*T/(p.hD+D)^2];
    J(4,:) = [p.eta*uI, p.eta*uI, 0, -p.dD];

end

function ok = admissible(y,p)

    ok = isreal(y) && all(isfinite(y)) && all(y>=0) && y(3)>0 && y(1)<=p.KA*(1+1e-10) && sum(y(1:2))<=p.KB*(1+1e-10);

end

function c = polyAdd(a,b)

    c = zeros(1,max(numel(a),numel(b)));
    c(1:numel(a)) = a;
    c(1:numel(b)) = c(1:numel(b))+b;

end