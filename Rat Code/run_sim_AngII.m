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
scenario = {'Normal', 'm_RSNA', 'm_AT2R', 'm_RAS', 'm_Reab', 'm_RSNA_&_m_Reab'};
num_scen = length(scenario);
fixed_ss = 1;

% Number of days to run simulation after change; Day at which to induce change;
days = 13; day_change = 1;
% Number of points for plotting resolution
N = ((days+1)*1440) / 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of variables
num_vars = 92;

% Initialize variables.
% X = (variables, points, gender, scenario)
X = zeros(num_vars,N,2,num_scen);

gender = {'male', 'female'};

for ss = 1:num_scen % scenario
for gg = 1:2        % gender

%% Parameters

% Parameter input
pars = get_pars(gender{gg}, scenario{ss});

%% Drugs

% drugs = [Ang II inf rate fmol/(ml min), ACEi target level, ARB target level]

if     strcmp(gender{gg}, 'male')
    drugs = [2022, 0, 0]; % Sampson 2008 male + female; 13 days
%     drugs = [(3/3)*5492, 0, 0]; % Zimmerman 2015 male + female; 14 days
elseif strcmp(gender{gg}, 'female')
    drugs = [2060, 0, 0]; % Sampson 2008 male + female; 13 days
%     drugs = [(2/3)*5492, 0, 0]; % Zimmerman 2015 male + female; 14 days
end

%% Solve DAE

% Initial value
% This initial condition is the steady state data value taken from
% Karaaslan 2005 and Leete 2018. 

% Set name for data file to be loaded based upon gender and scenario.    
load_data_name = sprintf('%s_ss_data_scenario_%s.mat', gender{gg},scenario{ss});

% Load data for steady state initial value. 
load(load_data_name, 'SSdata');

% Retrieve and replace parameters in fixed variable equations.
% These are the shift parameters which ensure that effect variables are 1.
fixed_ind = [2, 10, 14, 24, 44, 49, 66, 71, 88];
fixed_var_pars = SSdata(fixed_ind);
SSdata(fixed_ind) = 1;

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
          '$\Phi_{u-sod}$'; '$\Phi_{win}$'; '$V_{ecf}$'; '$V_{b}$'; ...
          '$P_{mf}$'; '$\Phi_{vr}$'; '$\Phi_{co}$'; '$P_{ra}$'; ...
          '$vas$'; '$vas_{f}$'; '$vas_{d}$'; '$R_{a}$'; '$R_{ba}$'; ...
          '$R_{vr}$'; '$R_{tp}$'; '$P_{ma}$'; '$\epsilon_{aum}$'; ...
          '$a_{auto}$'; '$a_{chemo}$'; '$a_{baro}$'; '$C_{adh}$'; ...
          '$N_{adh}$'; '$N_{adhs}$'; '$\delta_{ra}$'; ...
          '$\Phi_{pt-wreab}$'; '$\eta_{pt-wreab}$'; ...
          '$\mu_{pt-sodreab}$'; '$\Phi_{md-u}$'; '$\Phi_{dt-wreab}$'; ...
          '$\eta_{dt-wreab}$'; '$\mu_{dt-sodreab}$'; '$\Phi_{dt-u}$'; ...
          '$\Phi_{cd-wreab}$'; '$\eta_{cd-wreab}$'; ...
          '$\mu_{cd-sodreab}$'; '$\mu_{adh}$'; '$\Phi_{u}$'; ...
          '$M_{sod}$'; '$C_{sod}$'; '$\nu_{md-sod}$'; '$\nu_{rsna}$'; ...
          '$C_{al}$'; '$N_{al}$'; '$N_{als}$'; '$\xi_{k/sod}$'; ...
          '$\xi_{map}$'; '$\xi_{at}$'; '$\hat{C}_{anp}$'; '$AGT$'; ...
          '$\nu_{AT1}$'; '$R_{sec}$'; '$PRC$'; '$PRA$'; '$Ang I$'; ...
          '$Ang II$'; '$Ang II_{AT1R-bound}$'; '$Ang II_{AT2R-bound}$'; ...
          '$Ang (1-7)$'; '$Ang IV$'; '$R_{aa}$'; '$R_{ea}$'; ...
          '$\Sigma_{myo}$'; '$\Psi_{AT1R-AA}$'; '$\Psi_{AT1R-EA}$'; ...
          '$\Psi_{AT2R-AA}$'; '$\Psi_{AT2R-EA}$'};

% Initial condition for the variables and their derivatives. 
% System is initially at steady state, so the derivative is 0.
x0 = SSdata; x_p0 = zeros(num_vars,1);

% Factor by which to change something.
fact = 1;
fact_var = '';

% Time at which to keep steady state, change a parameter, etc.
tchange = day_change*1440;
% tchange = 10;

% Initial time (min); Final time (min);
t0 = 0*1440; tend = tchange + days*1440;
% t0 = 0; tend = tchange + 1000;

% Time vector
tspan = linspace(t0,tend,N);

% ode options
% options = odeset('RelTol',1e-1, 'AbsTol',1e-4); % default is -3, -6
options = odeset('MaxStep',1); % default is 0.1*abs(t0-tf)
% options = odeset('RelTol',1e-2, 'AbsTol',1e-4, 'MaxStep',1e-0);

% Solve dae
[t,x] = ode15i(@(t,x,x_p) ...
               bp_reg_sim(t,x,x_p,pars,fixed_var_pars,SSdata,drugs,...
                          tchange,fact,fact_var,scenario{ss})     ,...
               tspan, x0, x_p0, options);

% X = (variables, points, gender, scenario)
X(:,:,gg,ss) = x';

end % gender
end % scenario

%% Plot

% Retrieve male and female.
% X_m/f = (variables, points, scenario)
t = t';
X_m = reshape(X(:,:,1,:), [num_vars,N,num_scen]);
X_f = reshape(X(:,:,2,:), [num_vars,N,num_scen]);

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

f = gobjects(7,1);
s = gobjects(7,15);
% Loop through each set of subplots.
for i = 1:7
%     f(i) = figure; 
    f(i) = figure('pos',[750 500 650 450]);
    % This is to avoid the empty plots in the last subplot set.
    if i == 7
        last_plot = 2;
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
        
        ax.XTick = (tchange+0*(1) : 2 : tchange+days*(1));
        ax.XTickLabel = {'0' ,'2' ,'4' ,'6' ,'8' ,'10','12',...
                         '14','16','18','20','22','24','26'};

%         legend('Male', 'Female')
        title(names((i-1)*15 + j), 'Interpreter','latex', 'FontSize',15)
    end
end

% Plot Mean Arterial Pressure vs Time. ------------------------------------

% Data from Sampson 2008. MAP is in difference from baseline.
tdata     = [0+1  ,1+1  ,2+1  ,3+1  ,4+1  ,5+1  ,6+1  ,...
             7+1  ,8+1  ,9+1  ,10+1 ,11+1 ,12+1 ,13+1 ];
MAPdata_m = [0.035,7.218,18.33,19.48,17.76,14.59,19.58,...
             26.18,28.87,29.54,31.26,34.71,36.53,42.18];
MAPdata_f = [0.011,10.85,15.98,14.31,14.31,18.44,14.71,...
             13.91,17.31,17.04,18.37,19.63,23.23,24.42];
% Substract MAP by baseline for each gender and all scenarios.
% X_m/f = (variable, points, scenario)
MAP_m = reshape(X_m(42,:,:) - X_m(42,1,:), [N,num_scen]);
MAP_f = reshape(X_f(42,:,:) - X_f(42,1,:), [N,num_scen]);

g(1) = figure('DefaultAxesFontSize',14);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 2.5]);
plot(t,MAP_m(:,fixed_ss),'-', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3);
xlim([xlower, xupper]); ylim([0, 60]);
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

% GFR; BV; RSNA; REA/RR for each gender and all scenarios.
% X_m/f = (variable, points, scenario)
GFR_m  = reshape(X_m( 7,:,:), [N,num_scen]);
GFR_f  = reshape(X_f( 7,:,:), [N,num_scen]);
BV_m   = reshape(X_m(30,:,:), [N,num_scen]);
BV_f   = reshape(X_f(30,:,:), [N,num_scen]);
RSNA_m = reshape(X_m( 1,:,:), [N,num_scen]);
RSNA_f = reshape(X_f( 1,:,:), [N,num_scen]);
R_m    = reshape(X_m(87,:,:) ./ X_m( 4,:,:), [N,num_scen]);
R_f    = reshape(X_f(87,:,:) ./ X_f( 4,:,:), [N,num_scen]);
% Plot as relative change in order to compare male and female.
GFR_m = GFR_m ./ GFR_m(1,:);
GFR_f = GFR_f ./ GFR_f(1,:);
BV_m  = BV_m  ./ BV_m (1,:);
BV_f  = BV_f  ./ BV_f (1,:);
R_m   = R_m   ./ R_m  (1,:);
R_f   = R_f   ./ R_f  (1,:);

% Filtration fraction for sodium and urine for each gender and all scenarios.
FRNA_m = reshape((X_m(11,:,:) - X_m(27,:,:)) ./ X_m(11,:,:), [N,num_scen]) * 100;
FRNA_f = reshape((X_f(11,:,:) - X_f(27,:,:)) ./ X_f(11,:,:), [N,num_scen]) * 100;
FRW_m  = reshape((X_m( 7,:,:) - X_m(63,:,:)) ./ X_m( 7,:,:), [N,num_scen]) * 100;
FRW_f  = reshape((X_f( 7,:,:) - X_f(63,:,:)) ./ X_f( 7,:,:), [N,num_scen]) * 100;

h(1) = figure('DefaultAxesFontSize',14);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 2.5]);
plot(t,RSNA_m(:,fixed_ss),'-', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3);
xlim([xlower, xupper]); ylim([0.55,1.05]);
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
ylim(s_main(2), [97,100])
xlabel(s_main(2), 'Time (days)'); ylabel(s_main(2), 'FR (%)');
hold(s_main(2), 'on')
plot(s_main(2), t,FRW_m (:,fixed_ss), '--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(s_main(2), t,FRNA_f(:,fixed_ss), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(s_main(2), t,FRW_f (:,fixed_ss), '--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
fakeplot = zeros(2, 1);
fakeplot(1) = plot(s_main(2), NaN,NaN, 'k-' );
fakeplot(2) = plot(s_main(2), NaN,NaN, 'k--');
[~, hobj, ~, ~] = legend(fakeplot, {'FR_{Na^+}','FR_{U}'}, 'FontSize',7,'Location','Southeast');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
hold(s_main(2), 'off')
title(s_main(2), 'B')

plot(s_main(3), t,GFR_m (:,fixed_ss), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim(s_main(3), [xlower, xupper]);
set(s_main(3), 'XTick', [tchange+0*(1) : 2 : tchange+days*(1)]);
set(s_main(3), 'XTickLabel', {'0','2','4','6','8','10','12','14'});
ylim(s_main(3), [0.7,1.25])
xlabel(s_main(3), 'Time (days)'); ylabel(s_main(3), 'GFR (relative)');
hold(s_main(3), 'on')
plot(s_main(3), t,GFR_f (:,fixed_ss), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold(s_main(3), 'off')
title(s_main(3), 'C')

plot(s_main(4), t,BV_m  (:,fixed_ss), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim(s_main(4), [xlower, xupper]);
set(s_main(4), 'XTick', [tchange+0*(1) : 2 : tchange+days*(1)]);
set(s_main(4), 'XTickLabel', {'0','2','4','6','8','10','12','14'});
xlabel(s_main(4), 'Time (days)'); ylabel(s_main(4), 'BV (relative)');
hold(s_main(4), 'on')
plot(s_main(4), t,BV_f  (:,fixed_ss), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold(s_main(4), 'off')
title(s_main(4), 'D')

% Plot male - female bar graph for each scenario --------------------------

% Retrieve last MAP value for each gender and all scenarios.
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
ylim(s1(2), [0,15])
xlabel(s1(2), 'Scenario', 'FontSize',14); ylabel(s1(2), '\DeltaMAP (mmHg)', 'FontSize',14);
% hAxes.XAxis.FontSize = 6;
title(s1(2), 'B', 'FontSize',14)

% Save figures. -----------------------------------------------------------

save_data_name = sprintf('all_vars_AngII_inf.fig');
save_data_name = strcat('Figures/', save_data_name);
savefig([f;g;h';k], save_data_name)

end






























