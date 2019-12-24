% This simulates the blood pressure regulation model bp_reg.m for Ang II infusion.
% 
% Steady state data is calculated by solve_ss_baseline.m or solve_ss_scenario.m.

% function [SSdata, f] = run_sim
function run_sim_AngII

close all

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Begin user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scenarios
% Normal - Normal conditions
% m_RSNA - male RSNA
% m_AT2R - male AT2R
% m_RAS  - male RAS pars
% m_Reab - male fractional sodium and water reabsorption
% m_RAS_&_m_Reab - male RAS pars & fractional sodium and water reabsorption
scenario = {'Normal', 'm_RSNA', 'm_AT2R', 'm_RAS', 'm_Reab', ...
            'm_RSNA_m_Reab'};
% scenario = {'Normal', 'm_RSNA', 'm_AT2R', 'm_RAS', 'm_Reab', ...
%             'm_RSNA_m_Reab', ...
%             'Pri_Hyp'};
num_scen = length(scenario);
fixed_ss = 1;

% Species
spe_ind = 2;

% Number of days to run simulation after change; Day at which to induce change;
days = 13; day_change = 1;
% Number of points for plotting resolution
N = ((days+1)*1440) / 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

species = {'human', 'rat'   };
sex     = {'male' , 'female'};

% Number of variables
num_vars = 93;

% Initialize variables.
% X = (variables, points, sex, scenario)
X = zeros(num_vars,N,2,num_scen);

for sce_ind = 1:num_scen % scenario
for sex_ind = 1:2        % sex

varargin_input = {scenario{sce_ind},true};

%% Parameters

% Parameter input
pars = get_pars(species{spe_ind}, sex{sex_ind}, varargin_input{:});

%% Drugs

% Ang II inf rate fmol/(ml min)
if     strcmp(sex{sex_ind}, 'male')
%     kappa_AngII = 2022; % Sampson 2008
%     kappa_AngII = 785; % Sullivan 2010
    kappa_AngII = 630; % Sullivan 2010
elseif strcmp(sex{sex_ind}, 'female')
%     kappa_AngII = 2060; % Sampson 2008
%     kappa_AngII = 475; % Sullivan 2010
    kappa_AngII = 630; % Sullivan 2010
end

varargin_input = {varargin_input{:}, 'AngII',kappa_AngII};

%% Solve DAE

% Initial value
% This initial condition is the steady state data value taken from
% solve_ss_scenario.m.

% Set name for data file to be loaded based upon sex and scenario.    
load_data_name = sprintf('%s_%s_ss_data_scenario_%s.mat', ...
                         species{spe_ind},sex{sex_ind},scenario{sce_ind});
% Load data for steady state initial value. 
load(load_data_name, 'SSdata');

% Variable names for plotting.
names  = {'$rsna$'; '$\alpha_{map}$'; '$\alpha_{rap}$'; '$R_{r}$'; ...
          '$\beta_{rsna}$'; '$\Phi_{rb}$'; '$\Phi_{gfilt}$'; '$P_{f}$'; ...
          '$P_{gh}$'; '$\Sigma_{tgf}$'; '$\Phi_{filsod}$'; ...
          '$\Phi_{pt-sodreab}$'; '$\eta_{pt-sodreab}$'; ...
          '$\gamma_{filsod}$'; '$\gamma_{at}$'; '$\gamma_{rsna}$'; ...
          '$\Phi_{md-sod}$'; '$\Phi_{dt-sodreab}$'; ...
          '$\eta_{dt-sodreab}$'; '$\psi_{al}$'; '$\Phi_{dt-sod}$'; ...
          '$\Phi_{cd-sodreab}$'; '$\eta_{cd-sodreab}$'; ...
          '$\lambda_{dt}$'; '$\lambda_{anp}$'; '$\lambda_{al}$'; ...
          '$\Phi_{u-sod}$'; '$\Phi_{sodin}$'; '$V_{ecf}$'; '$V_{b}$'; ...
          '$P_{mf}$'; '$\Phi_{vr}$'; '$\Phi_{co}$'; '$P_{ra}$'; ...
          '$vas$'; '$vas_{f}$'; '$vas_{d}$'; '$R_{a}$'; '$R_{ba}$'; ...
          '$R_{vr}$'; '$R_{tp}$'; '$P_{ma}$'; '$\epsilon_{aum}$'; ...
          '$a_{auto}$'; '$a_{chemo}$'; '$a_{baro}$'; '$C_{adh}$'; ...
          '$N_{adh}$'; '$N_{adhs}$'; '$\delta_{ra}$'; ...
          '$M_{sod}$'; '$C_{sod}$'; '$\nu_{md-sod}$'; '$\nu_{rsna}$'; ...
          '$C_{al}$'; '$N_{al}$'; '$N_{als}$'; '$\xi_{k/sod}$'; ...
          '$\xi_{map}$'; '$\xi_{at}$'; '$\hat{C}_{anp}$'; '$AGT$'; ...
          '$\nu_{AT1}$'; '$R_{sec}$'; '$PRC$'; '$PRA$'; '$Ang I$'; ...
          '$Ang II$'; '$Ang II_{AT1R-bound}$'; '$Ang II_{AT2R-bound}$'; ...
          '$Ang (1-7)$'; '$Ang IV$'; '$R_{aa}$'; '$R_{ea}$'; ...
          '$\Sigma_{myo}$'; '$\Psi_{AT1R-AA}$'; '$\Psi_{AT1R-EA}$'; ...
          '$\Psi_{AT2R-AA}$'; '$\Psi_{AT2R-EA}$'; ...
          '$\Phi_{pt-wreab}$'; '$\eta_{pt-wreab}$'; ...
          '$\mu_{pt-sodreab}$'; '$\Phi_{md-u}$'; '$\Phi_{dt-wreab}$'; ...
          '$\eta_{dt-wreab}$'; '$\mu_{dt-sodreab}$'; '$\Phi_{dt-u}$'; ...
          '$\Phi_{cd-wreab}$'; '$\eta_{cd-wreab}$'; ...
          '$\mu_{cd-sodreab}$'; '$\mu_{adh}$'; ...
          '$\Phi_{u}$'; '$\Phi_{win}$'};

% Initial condition for the variables and their derivatives. 
% System is initially at steady state, so the derivative is 0.
x0 = SSdata; x_p0 = zeros(num_vars,1);

% Time at which to keep steady state, change a parameter, etc.
tchange = day_change*1440;

% Initial time (min); Final time (min);
t0 = 0*1440; tend = tchange + days*1440;

% Time vector
tspan = linspace(t0,tend,N);

% ode options
options = odeset('MaxStep',1); % default is 0.1*abs(t0-tf)

% Solve dae
[t,x] = ode15i(@(t,x,x_p) ...
               bp_reg_mod(t,x,x_p,pars,tchange,varargin_input{:}), ...
               tspan, x0, x_p0, options);

% X = (variables, points, sex, scenario)
X(:,:,sex_ind,sce_ind) = x';

end % sex
end % scenario

%% Plot

% Retrieve male and female.
% X_m/f = (variables, points, scenario)
t = t';
X_m = reshape(X(:,:,1,:), [num_vars,N,num_scen]); 
X_f = reshape(X(:,:,2,:), [num_vars,N,num_scen]); 
% X_m = X_f;
% X_f = X_m;

% x-axis limits
xlower = t0; xupper = tend; 

% Convert minutes to days for longer simulations.
t = t/1440; tchange = tchange/1440; 
xlower = xlower/1440; xupper = xupper/1440; 

% y-axis limits
ylower = zeros(num_vars,1); yupper = ylower; 
for i = 1:length(ylower)
    ylower(i) = 0.95*min( min(X_m(i,:,fixed_ss)), min(X_f(i,:,fixed_ss)) );
    yupper(i) = 1.05*max( max(X_m(i,:,fixed_ss)), max(X_f(i,:,fixed_ss)) );
    if abs(yupper(i)) < eps*100
        ylower(i) = -10^(-5); yupper(i) = 10^(-5);
    end
end

% Interesting variables to plot.
var_ind = [33;41;42;9;73;74;6;7;27;92;93;29]; sub_var_num = length(var_ind);

% Plot all vars vs time. --------------------------------------------------

f = gobjects(7,1);
s = gobjects(7,15);
% Loop through each set of subplots.
for i = 1:7
    f(i) = figure('pos',[750 500 650 450]);
    % This is to avoid the empty plots in the last subplot set.
    if i == 7
        last_plot = mod(num_vars, 15);
    else
        last_plot = 15;
    end
    % Loop through each subplot within a set of subplots.
    for j = 1:last_plot
        s(i,j) = subplot(3,5,j);
        s(i,j).Position = s(i,j).Position + [0 0 0.01 0];
        
        plot(s(i,j), t,X_m((i-1)*15 + j,:,fixed_ss),'b', ...
                     t,X_f((i-1)*15 + j,:,fixed_ss),'r');
        
        xlim([xlower, xupper])
        ylim([ylower((i-1)*15 + j), yupper((i-1)*15 + j)])
        
%         Days
        ax.XTick = (tchange+0*(1) : 2 : tchange+days*(1));
        ax.XTickLabel = {'0','2','4','6','8','10','12','14'};

%         legend('Male', 'Female')
        title(names((i-1)*15 + j), 'Interpreter','latex', 'FontSize',15)
    end
end

% Plot interesting variables. ---------------------------------------------

f2 = figure('pos',[000 000 600 600], 'DefaultAxesFontSize',12);
s2 = gobjects(1,sub_var_num);
% Loop through each subplot within a set of subplots.
for j = 1:sub_var_num
    s2(j) = subplot(4,3,j);
    if     mod(j,3) == 1
        hshift = -0.05;
    elseif mod(j,3) == 0
        hshift = 0.05;
    else
        hshift = 0;
    end
    s2(j).Position = s2(j).Position + [hshift 0 0.01 0.01];

    plot(s2(j), t,X_m(var_ind(j),:,fixed_ss), 'Color',[0.203, 0.592, 0.835], 'LineWidth',2.5);
    hold(s2(j), 'on')
    plot(s2(j), t,X_f(var_ind(j),:,fixed_ss), 'Color',[0.835, 0.203, 0.576], 'LineWidth',2.5);
    hold(s2(j), 'off')

    xlim([xlower, xupper])
    ylim([ylower(var_ind(j)), yupper(var_ind(j))])
    
    set(s2(j), 'XTick', [tchange+0*(1) : 2 : tchange+days*(1)]);
    set(s2(j), 'XTickLabel', {'0','2','4','6','8','10','12','14'});
    
    if j == 10 || j == 11
        ax = gca;
        ax.YAxis.Exponent = -3;
    end

    ylabel(names(var_ind(j)), 'Interpreter','latex', 'FontSize',16)
end
legend(s2(1),'Male','Female', 'Location','east')
xlh = xlabel(s2(11),'Time (days)');
xlh.Position(2) = xlh.Position(2) - 0.0005;

% Plot Mean Arterial Pressure vs Time. ------------------------------------

% Data from Sampson 2008. MAP is in difference from baseline.
tdata     = [0+1  ,1+1  ,2+1  ,3+1  ,4+1  ,5+1  ,6+1  ,...
             7+1  ,8+1  ,9+1  ,10+1 ,11+1 ,12+1 ,13+1 ];
MAPdata_m = [0.035,7.218,18.33,19.48,17.76,14.59,19.58,...
             26.18,28.87,29.54,31.26,34.71,36.53,42.18];
MAPdata_f = [0.011,10.85,15.98,14.31,14.31,18.44,14.71,...
             13.91,17.31,17.04,18.37,19.63,23.23,24.42];
% % Data from Sullivan 2010. MAP is in difference from baseline.
% tdata     = [0+1 , 1+1 , 2+1 , 3+1 , 4+1 , 5+1 , 6+1 ,...
%              7+1 , 8+1 , 9+1 , 10+1, 11+1, 12+1, 13+1, 14+1];
% MAPdata_m = [0   , -1.3, 2.3 , 8.9 , 15.5, 18.3, 22.7, 22.6, ...
%              28.6, 31.2, 30.9, 32.8, 37.4, 41.4, 40.3];
% MAPdata_f = [0   , 5.2 ,  5.3,  3.9,  3.6,  5.9,    8,   13, ...
%              15.7, 17.4, 19.8, 23.7, 25.8,  23.5,  24];

% Substract MAP by baseline for each sex and all scenarios.
% X_m/f = (variable, points, scenario)
MAP_m = zeros(N,num_scen); MAP_f = zeros(N,num_scen);
for i = 1:num_scen
    MAP_m(:,i) = X_m(42,:,i) - X_m(42,1,i);
    MAP_f(:,i) = X_f(42,:,i) - X_f(42,1,i);
end
% MAP_m = reshape(X_m(42,:,i) - X_m(42,1,i), [N,num_scen]);
% MAP_f = reshape(X_f(42,:,i) - X_f(42,1,i), [N,num_scen]);

g(1) = figure('DefaultAxesFontSize',14);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 2.5]);
plot(t,MAP_m(:,fixed_ss),'-', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3);
xlim([xlower, xupper]); ylim([-1, 60]);
ax = gca;
ax.XTick = (tchange+0*(1) : 2 : tchange+days*(1));
ax.XTickLabel = {'0','2','4','6','8','10','12','14'};
xlabel('Time (days)'); ylabel('\DeltaMAP (mmHg)');
hold on
plot(t,MAP_f(:,fixed_ss),'-', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3)
plot(tdata,MAPdata_m,'o', 'DisplayName',  'Male data', 'Color',[0.203, 0.592, 0.835], 'MarkerSize',6, 'LineWidth',2)
plot(tdata,MAPdata_f,'o', 'DisplayName','Female data', 'Color',[0.835, 0.203, 0.576], 'MarkerSize',6, 'LineWidth',2)
[~, hobj, ~, ~] = legend({'Male sim','Female sim','Male data','Female data'}, 'FontSize',7,'Location','Northwest');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);

% Plot all other quantities of interest. ----------------------------------

% GFR; BV; RSNA; REA/RR for each sex and all scenarios.
% X_m/f = (variable, points, scenario)
GFR_m  = reshape(X_m( 7,:,:), [N,num_scen]);
GFR_f  = reshape(X_f( 7,:,:), [N,num_scen]);
BV_m   = reshape(X_m(30,:,:), [N,num_scen]);
BV_f   = reshape(X_f(30,:,:), [N,num_scen]);
RSNA_m = reshape(X_m( 1,:,:), [N,num_scen]);
RSNA_f = reshape(X_f( 1,:,:), [N,num_scen]);
R_m    = reshape(X_m(74,:,:) ./ X_m( 4,:,:), [N,num_scen]);
R_f    = reshape(X_f(74,:,:) ./ X_f( 4,:,:), [N,num_scen]);
% Plot as relative change in order to compare male and female.
GFR_m_bl = GFR_m(1,:);
GFR_f_bl = GFR_f(1,:);
BV_m_bl  = BV_m (1,:);
BV_f_bl  = BV_f (1,:);
R_m_bl   = R_m  (1,:);
R_f_bl   = R_f  (1,:);
for i = 1:N
    GFR_m(i,:) = GFR_m(i,:) ./ GFR_m_bl;
    GFR_f(i,:) = GFR_f(i,:) ./ GFR_f_bl;
    BV_m (i,:) = BV_m (i,:) ./ BV_m_bl ;
    BV_f (i,:) = BV_f (i,:) ./ BV_f_bl ;
    R_m  (i,:) = R_m  (i,:) ./ R_m_bl  ;
    R_f  (i,:) = R_f  (i,:) ./ R_f_bl  ;
end
% GFR_m = GFR_m ./ GFR_m(1,:);
% GFR_f = GFR_f ./ GFR_f(1,:);
% BV_m  = BV_m  ./ BV_m (1,:);
% BV_f  = BV_f  ./ BV_f (1,:);
% R_m   = R_m   ./ R_m  (1,:);
% R_f   = R_f   ./ R_f  (1,:);

% Filtration fraction for sodium and urine for each sex and all scenarios.
FRNA_m = reshape((X_m(11,:,:) - X_m(27,:,:)) ./ X_m(11,:,:), [N,num_scen]) * 100;
FRNA_f = reshape((X_f(11,:,:) - X_f(27,:,:)) ./ X_f(11,:,:), [N,num_scen]) * 100;
FRW_m  = reshape((X_m( 7,:,:) - X_m(92,:,:)) ./ X_m( 7,:,:), [N,num_scen]) * 100;
FRW_f  = reshape((X_f( 7,:,:) - X_f(92,:,:)) ./ X_f( 7,:,:), [N,num_scen]) * 100;
% Plot as relative change in order to compare male and female.
FRNA_m_bl = FRNA_m(1,:);
FRNA_f_bl = FRNA_f(1,:);
FRW_m_bl  = FRW_m (1,:);
FRW_f_bl  = FRW_f (1,:);
for i = 1:N
    FRNA_m(i,:) = FRNA_m(i,:) ./ FRNA_m_bl;
    FRNA_f(i,:) = FRNA_f(i,:) ./ FRNA_f_bl;
    FRW_m (i,:) = FRW_m (i,:) ./ FRW_m_bl ;
    FRW_f (i,:) = FRW_f (i,:) ./ FRW_f_bl ;
end
% FRNA_m = FRNA_m ./ FRNA_m(1,:);
% FRNA_f = FRNA_f ./ FRNA_f(1,:);
% FRW_m  = FRW_m  ./ FRW_m (1,:);
% FRW_f  = FRW_f  ./ FRW_f (1,:);

h(1) = figure('DefaultAxesFontSize',14);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 2.5]);
plot(t,RSNA_m(:,fixed_ss),'-', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3);
xlim([xlower, xupper]); ylim([0.50,1.05]);
ax = gca;
ax.XTick = (tchange+0*(1) : 2 : tchange+days*(1));
ax.XTickLabel = {'0','2','4','6','8','10','12','14'};
xlabel('Time (days)'); ylabel('RSNA');
hold on
plot(t,RSNA_f(:,fixed_ss),'-', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3)
[~, hobj, ~, ~] = legend({'Male','Female'}, 'FontSize',7,'Location','Northeast');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
% ---
h(2) = figure('DefaultAxesFontSize',14);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.15, 5]);
s_main(1) = subplot(2,2,1); 
s_main(2) = subplot(2,2,2); 
s_main(3) = subplot(2,2,3);
s_main(4) = subplot(2,2,4); 

plot(s_main(1), t,R_m   (:,fixed_ss), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim(s_main(1), [xlower, xupper]);
set(s_main(1), 'XTick', [tchange+0*(1) : 2 : tchange+days*(1)]);
set(s_main(1), 'XTickLabel', {'0','2','4','6','8','10','12','14'});
ylim(s_main(1), [0.75,1.05])
xlabel(s_main(1), 'Time (days)'); ylabel(s_main(1), 'R_{EA}/R_R (relative)');
hold(s_main(1), 'on')
plot(s_main(1), t,R_f   (:,fixed_ss), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold(s_main(1), 'off')
[~, hobj, ~, ~] = legend(s_main(1), {'Male','Female'}, 'FontSize',7,'Location','Northeast');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
title(s_main(1), 'A')

plot(s_main(2), t,FRNA_m(:,fixed_ss) ,'-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim(s_main(2), [xlower, xupper]);
set(s_main(2), 'XTick', [tchange+0*(1) : 2 : tchange+days*(1)]);
set(s_main(2), 'XTickLabel', {'0','2','4','6','8','10','12','14'});
% ylim(s_main(2), [97,100])
xlabel(s_main(2), 'Time (days)'); ylabel(s_main(2), 'FR (relative)');
hold(s_main(2), 'on')
plot(s_main(2), t,FRW_m (:,fixed_ss), '--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(s_main(2), t,FRNA_f(:,fixed_ss), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(s_main(2), t,FRW_f (:,fixed_ss), '--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
fakeplot = zeros(2, 1);
fakeplot(1) = plot(s_main(2), NaN,NaN, 'k-' );
fakeplot(2) = plot(s_main(2), NaN,NaN, 'k--');
[~, hobj, ~, ~] = legend(fakeplot, {'FR_{Na^+}','FR_{U}'}, 'FontSize',7,'Location','Northeast');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
hold(s_main(2), 'off')
title(s_main(2), 'B')

plot(s_main(3), t,GFR_m (:,fixed_ss), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim(s_main(3), [xlower, xupper]);
set(s_main(3), 'XTick', [tchange+0*(1) : 2 : tchange+days*(1)]);
set(s_main(3), 'XTickLabel', {'0','2','4','6','8','10','12','14'});
ylim(s_main(3), [0.75,1.35])
xlabel(s_main(3), 'Time (days)'); ylabel(s_main(3), 'GFR (relative)');
hold(s_main(3), 'on')
plot(s_main(3), t,GFR_f (:,fixed_ss), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold(s_main(3), 'off')
title(s_main(3), 'C')

plot(s_main(4), t,BV_m  (:,fixed_ss), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim(s_main(4), [xlower, xupper]);
set(s_main(4), 'XTick', [tchange+0*(1) : 2 : tchange+days*(1)]);
set(s_main(4), 'XTickLabel', {'0','2','4','6','8','10','12','14'});
ylim(s_main(4), [1,1.25])
xlabel(s_main(4), 'Time (days)'); ylabel(s_main(4), 'BV (relative)');
hold(s_main(4), 'on')
plot(s_main(4), t,BV_f  (:,fixed_ss), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold(s_main(4), 'off')
title(s_main(4), 'D')

% Plot male - female bar graph for each scenario --------------------------

% Retrieve last MAP value for each sex and all scenarios.
% X_m/f = (variable, points, scenario)
deltaMAP_m = reshape(X_m(42,end,1:end) - X_m(42,1,1:end), [1,num_scen]);
deltaMAP_f = reshape(X_f(42,end,1:end) - X_f(42,1,1:end), [1,num_scen]);
% Substract all MAP values from all scenarios from normal male MAP.
MAP_comp = deltaMAP_m(1) - [deltaMAP_m(1), deltaMAP_f];
% String for bar graph labels.
scen_comp = categorical({'M - M'              , 'M - F', ...
                         'M - F M RSNA', 'M - F M AT2R', ...
                         'M - F M RAS' , 'M - F M Reab', ...
                         'M - FM RSNA\newline& Reab'});
scen_comp = reordercats(scen_comp,{'M - M'              , 'M - F', ...
                                   'M - F M RSNA', 'M - F M AT2R', ...
                                   'M - F M RAS' , 'M - F M Reab', ...
                                   'M - FM RSNA\newline& Reab'});

k = figure('DefaultAxesFontSize',10);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.15, 3.5]);
s1(1) = subplot(1,2,1); 
s1(2) = subplot(1,2,2); 
% s1(2).Position = s1(2).Position + [0.0, 0.0, 0.0, -0.01];

plot(s1(1), t,MAP_m(:,fixed_ss),'-', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3);
xlim(s1(1), [xlower, xupper]); ylim(s1(1), [0, 60]);
xticks(s1(1), (tchange+0*(1) : 2 : tchange+days*(1))); xticklabels(s1(1), {'0','2','4','6','8','10','12','14'});
xlabel(s1(1), 'Time (days)', 'FontSize',14*1.1); ylabel(s1(1), '\DeltaMAP (mmHg)', 'FontSize',14*1.1);
hold(s1(1), 'on')
plot(s1(1), t,MAP_f(:,fixed_ss),'-', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3)
plot(s1(1), tdata,MAPdata_m,'o', 'DisplayName',  'Male data', 'Color',[0.203, 0.592, 0.835], 'MarkerSize',6, 'LineWidth',2)
plot(s1(1), tdata,MAPdata_f,'o', 'DisplayName','Female data', 'Color',[0.835, 0.203, 0.576], 'MarkerSize',6, 'LineWidth',2)
[~, hobj, ~, ~] = legend(s1(1), {'Male sim','Female sim','Male data','Female data'}, 'FontSize',7,'Location','Northwest');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
title(s1(1), 'A', 'FontSize',14)

bar(s1(2), scen_comp,MAP_comp,'k');
% set(gca,'xticklabel',scen_comp_text);
% xtickget = get(gca,'xticklabel');  
% set(gca,'xticklabel',xtickget,'fontsize',6)
% xtickangle(s1(2),90)
% xlim(s1(2), [1-1,6+1])
ylim(s1(2), [0,20])
xlabel(s1(2), 'Scenario', 'FontSize',14); ylabel(s1(2), '\DeltaMAP (mmHg)', 'FontSize',14);
% hAxes.XAxis.FontSize = 6;
title(s1(2), 'B', 'FontSize',14)

% Save figures. -----------------------------------------------------------

save_data_name = sprintf('all_vars_AngII_inf.fig');
save_data_name = strcat('Figures/', save_data_name);
savefig([f;f2;g;h';k], save_data_name)

end






























