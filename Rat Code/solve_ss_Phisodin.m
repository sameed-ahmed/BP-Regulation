% This script vairies the sodium intake and then solves for the steady 
% state values of the blood pressure regulation model 
% bp_reg_solve_Phisodin.m for each inputted value. 
% 
% Steady state data for the intial guess is calculated by 
% solve_ss_baseline.m or solve_ss_scenario.m.
% 
% All variables are then plotted versus the relative change in input.

function solve_ss_Phisodin

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
% ACEi   - Angiotensin convernting enzyme inhibitor
% AngII  - Ang II infusion
scenario = {'Normal', 'm_RSNA', 'm_AT2R', 'm_RAS', 'm_Reab', 'ACEi', 'AngII'};
num_scen = length(scenario);
% Index of scenario to plot for all variables
fixed_ss = 1;

% Boolean to fix/vary water intake.
% win =  'fixed';
win = 'varied';

% Number of iterations below/above baseline.
iteration = 51; % must be odd number for symmetry
% Fold decrease/increase.
lower = 1/5; upper = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of variables changes based on whether water intake is
% fixed/varied.
if     strcmp(win,  'fixed')
    num_vars = 92-1;
elseif strcmp(win, 'varied')
    num_vars = 92;
end

% Baseline of water intake for each scenario if it is fixed.
Phi_win_bl = zeros(2,num_scen,2*iteration-1);

% Range for fold decrease/increase.
iter_range_l = linspace(lower, 1, iteration);
iter_range_u = linspace(1, upper, iteration);
% Delete overlapping baseline value.
iter_range_l(end) = '';
% Combine decrease/increase ranges into single vector.
iter_range   = [iter_range_l, iter_range_u];

% Initialize variables.
% X = (variable, iteration, gender, scenario)
X = zeros(num_vars,2*iteration-1,2,num_scen);
% Retrieve male/female. 
% X_m/f = (variable, iteration, scenario)
if     strcmp(win,  'fixed')
    X_m = zeros(num_vars+1,2*iteration-1,num_scen);
    X_f = zeros(num_vars+1,2*iteration-1,num_scen);
elseif strcmp(win, 'varied')
    X_m = zeros(num_vars,2*iteration-1,num_scen);
    X_f = zeros(num_vars,2*iteration-1,num_scen);
end

gender = {'male'    , 'female'  };
change = {'decrease', 'increase'};

for ss = 1:num_scen % scenario
for gg = 1:2        % gender
for cc = 1:2        % change

% Set name for data file to be loaded based upon gender and scenario.    
load_data_name = sprintf('%s_ss_data_scenario_%s.mat', gender{gg},scenario{ss});

% Retrieve and replace parameters in fixed variable equations.
% Load data for steady state initial guess. 
% Fixed variable parameters are only retrieved and replaced for scenarios
% which are solved for in the solve_ss_baseline because the parameter has
% changed for the fixed variable to remain 1.
% Otherwise scenarios which are solved for in solve_ss_scenario do not 
% require this because they load the fixed variable parameter from the
% baseline scenario, and they are perturbed scenarios in which the fixed
% variable is no longer 1.
fixed_ind = [2, 10, 14, 24, 44, 49, 66, 71, 88];
if   strcmp(scenario{ss}, 'ACEi') || strcmp(scenario{ss}, 'AngII')
    if     strcmp(gender{gg}, 'male')
        load(  'male_ss_data_scenario_Normal.mat', 'SSdata');
    elseif strcmp(gender{gg}, 'female')
        load('female_ss_data_scenario_Normal.mat', 'SSdata');
    end
    fixed_var_pars = SSdata(fixed_ind);
    phicophico = SSdata(33); cadhcadh = SSdata(47);
    fixed_var_pars = [fixed_var_pars; cadhcadh; phicophico];
    load(load_data_name, 'SSdata');
else
    load(load_data_name, 'SSdata');
    fixed_var_pars = SSdata(fixed_ind);
    phicophico = SSdata(33); cadhcadh = SSdata(47);
    fixed_var_pars = [fixed_var_pars; cadhcadh; phicophico];
    SSdata(fixed_ind) = 1;
end
SSdataIG = SSdata;
clear SSdata;

% Load data for baseline water intake if it is fixed.
Phi_win_bl(gg,ss,1) = SSdataIG(28);
% Input Phi_win if it is fixed.
Phi_win_input = Phi_win_bl(gg,ss,1);

% Delete Phi_win if it is fixed.
if     strcmp(win,  'fixed')
    SSdataIG(28) = '';
elseif strcmp(win, 'varied')
end

for iter = 1:iteration % range

%% Parameters

% Baseline/range of sodium intake.
Phi_sodin_bl_m = 1.2212;
Phi_sodin_bl_f = 1.2212; % CHECK IF DIFFERENT LATER.
Phi_sodin_range_m = Phi_sodin_bl_m * iter_range;
Phi_sodin_range_f = Phi_sodin_bl_f * iter_range;

% Vary sodium intake.
if     strcmp(gender{gg}, 'male')
    if     strcmp(change{cc}, 'decrease')
        Phi_sodin = Phi_sodin_range_m(iteration-iter+1);
    elseif strcmp(change{cc}, 'increase')
        Phi_sodin = Phi_sodin_range_m(iteration+iter-1);
    end
elseif strcmp(gender{gg}, 'female')
    if     strcmp(change{cc}, 'decrease')
        Phi_sodin = Phi_sodin_range_f(iteration-iter+1);
    elseif strcmp(change{cc}, 'increase')
        Phi_sodin = Phi_sodin_range_f(iteration+iter-1);
    end
end

% Parameter input
pars     = get_pars(gender{gg}, scenario{ss});
pars(18) = Phi_sodin;

%% Drugs

% drugs = [Ang II inf rate fmol/(ml min), ACEi target level, ARB target level]
if     strcmp(scenario{ss}, 'ACEi'  )
        drugs = [0   , 1, 0]; % Hall 2018
elseif strcmp(scenario{ss}, 'AngII' )
    if     strcmp(gender{gg}, 'male'  )
        drugs = [2022, 0, 0]; % Sampson 2008
    elseif strcmp(gender{gg}, 'female')
        drugs = [2060, 0, 0]; % Sampson 2008
    end
else
        drugs = [0   , 0, 0];
end

%% Variables initial guess

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

% Initial guess for the variables.
% Find the steady state solution, so the derivative is 0.
% Arbitrary value for time to input.
x0 = SSdataIG; x_p0 = zeros(num_vars,1); t = 0;

%% Find steady state solution

% options = optimset(); %options = optimset('MaxFunEvals',8100+10000);
options = optimset('Display','off');
[SSdata, ~, ...
exitflag, output] = fsolve(@(x) bp_reg_solve_Phisodin(t,x,x_p0,pars ,...
                                                      fixed_var_pars,...
                                                      drugs,win     ,...
                                                      Phi_win_input ,...
                                                      scenario{ss}) ,...
                           x0, options);

% Check for solver convergence.
if exitflag == 0
    disp('Solver did not converge.')
    disp(output)
end

% Check for imaginary solution.
if not (isreal(SSdata))
    disp('Imaginary number returned.')
end

% Set any values that are within machine precision of 0 equal to 0.
for i = 1:length(SSdata)
    if abs(SSdata(i)) < eps*100
        SSdata(i) = 0;
    end
end

% Store solution.
if     strcmp(change{cc}, 'decrease')
    X(:,iteration-iter+1,gg,ss) = SSdata;
elseif strcmp(change{cc}, 'increase')
    X(:,iteration+iter-1,gg,ss) = SSdata;
end

% To avoid the solver not converging, the initial guess for solution to the
% system is taken as the previous solution value. That is, IG_i = SOL_i-1.
% Update next initial guess as current solution.
SSdataIG = SSdata;

% Sanity check to see script's progress. Also a check for where to
% troubleshoot in case the solver does not converge.
fprintf('%s %s %s iteration = %s out of %s \n', ...
        scenario{ss},gender{gg},change{cc},num2str(iter),num2str(iteration))

end % change
end % range
end % gender

%% Retrieve data and visualize

% X = (variable, iteration, gender, scenario)
% Phi_win_range = [scenario, iteration]

% Add in Phi_win where it originally was if it is fixed.
if     strcmp(win,  'fixed')
    Phi_win_bl(:,ss,:) = Phi_win_bl(:,ss,1) .* ones(1,1,2*iteration-1);
    X_m(:,:,ss) = [X(1:27,:,1,ss); Phi_win_bl(1,ss,:); X(28:end,:,1,ss)];
    X_f(:,:,ss) = [X(1:27,:,2,ss); Phi_win_bl(2,ss,:); X(28:end,:,2,ss)];
elseif strcmp(win, 'varied')
    X_m(:,:,ss) = X(:,:,1,ss); X_f(:,:,ss) = X(:,:,2,ss);
end

end % scenario

% x-axis
xscale = iter_range;

% % Plot relative change/relative change vs relative change to see slope
% xscale = xscale - 1;
% X_m(:,:,1) = (X_m(:,:,1) - X_m(:,iteration,1)) ./ X_m(:,iteration,1);
% X_f(:,:,1) = (X_f(:,:,1) - X_f(:,iteration,1)) ./ X_f(:,iteration,1);
% X_m(:,:,1) = X_m(:,:,1) ./ xscale;
% X_f(:,:,1) = X_f(:,:,1) ./ xscale;

% y-axis limits
% X_f = X_m;
% X_m = X_f;
ylower = zeros(length(X_m(:,1,1)),1); yupper = ylower; 
for i = 1:length(ylower)
    ylower(i) = 0.9*min(min(X_m(i,:,1)), min(X_f(i,:,1)));
    yupper(i) = 1.1*max(max(X_m(i,:,1)), max(X_f(i,:,1)));
    if ylower(i) == yupper(i)
        ylower(i) = ylower(i) - 10^(-5); yupper(i) = yupper(i) + 10^(-5);
    end
    if ylower(i) == inf || yupper(i) == inf || isnan(ylower(i)) || isnan(yupper(i))
        ylower(i) = -1; yupper(i) = 1;
    end
end

% X_f = zeros(size(X_f));
% X_m = zeros(size(X_m));

% Plot all variables vs sodium intake. ------------------------------------

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

        plot(s(i,j), xscale,X_m((i-1)*15 + j,:,fixed_ss),'b', ...
                     xscale,X_f((i-1)*15 + j,:,fixed_ss),'r');
        
        xlim([lower, upper])
        ylim([ylower((i-1)*15 + j), yupper((i-1)*15 + j)])
        
        xlabel('$\Phi_{sodin}$', 'Interpreter','latex', 'FontSize',15)
        title(names((i-1)*15 + j), 'Interpreter','latex', 'FontSize',15)
%         legend('Male', 'Female')
    end
end

% Plot Sodium Intake vs Mean Arterial Pressure. ---------------------------

g = figure('DefaultAxesFontSize',14);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 3.5]);
plot(X_m(42,:,fixed_ss),xscale,'-', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3);
% xlim([80, 120])
ylim([lower, upper])
ax = gca;
% ax.XTick = (80 : 10 : 120);
xlabel('MAP (mmHg)')
ylabel({'Fold change in'; 'sodium excretion'})
hold on
plot(X_f(42,:,fixed_ss),xscale,'-', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3)
legend('Male','Female', 'Location','Northwest')
hold off

% Plot with different scenarios. ------------------------------------------

h = figure('pos',[100 100 675 450]);
plot(X_m(42,:,1),xscale,'b-' , 'LineWidth',3, 'DisplayName','M Normal')
% xlim([80, 160])
ylim([lower, upper])
set(gca,'FontSize',14)
xlabel(names(42)       , 'Interpreter','latex', 'FontSize',22, 'FontWeight','bold')
ylabel('$\Phi_{sodin}$', 'Interpreter','latex', 'FontSize',22, 'FontWeight','bold')
legend('-DynamicLegend');
hold all
plot(X_f(42,:,1),xscale,'r-', 'LineWidth',3, 'DisplayName','F Normal')
legend('-DynamicLegend');

scen_linestyle_m = {'b-x', 'b-o', 'b-+', 'b-*',};
scen_linestyle_f = {'r-x', 'r-o', 'r-+', 'r-*',};
scen_legend = {'RSNA', 'AT2R', 'RAS', 'Reab'};
for ss = 2:num_scen-2
    hold all
    plot(X_m(42,:,ss),xscale,scen_linestyle_m{ss-1}, 'LineWidth',3, 'DisplayName',scen_legend{ss-1})
    legend('-DynamicLegend');
    hold all
    plot(X_f(42,:,ss),xscale,scen_linestyle_f{ss-1}, 'LineWidth',3, 'DisplayName',scen_legend{ss-1})
    legend('-DynamicLegend');
end

i = figure('pos',[100 100 675 450]);
plot(X_m(42,:,1),xscale,'b-' , 'LineWidth',3, 'DisplayName','M Normal')
% xlim([80, 160])
ylim([lower, upper])
set(gca,'FontSize',14)
xlabel(names(42)       , 'Interpreter','latex', 'FontSize',22, 'FontWeight','bold')
ylabel('$\Phi_{sodin}$', 'Interpreter','latex', 'FontSize',22, 'FontWeight','bold')
legend('-DynamicLegend');
hold all
plot(X_f(42,:,1),xscale,'r-', 'LineWidth',3, 'DisplayName','F Normal')
legend('-DynamicLegend');

scen_linestyle_m = {'b--', 'b:'};
scen_linestyle_f = {'r--', 'r:'};
scen_legend = {'ACEi', 'Ang II'};
for ss = num_scen-1:num_scen
    hold all
    plot(X_m(42,:,ss),xscale,scen_linestyle_m{ss-(num_scen-2)}, 'LineWidth',3, 'DisplayName',scen_legend{ss-(num_scen-2)})
    legend('-DynamicLegend');
    hold all
    plot(X_f(42,:,ss),xscale,scen_linestyle_f{ss-(num_scen-2)}, 'LineWidth',3, 'DisplayName',scen_legend{ss-(num_scen-2)})
    legend('-DynamicLegend');
end

% Plot male - female bar graph for each scenario. -------------------------

% X_m/f = (variable, iteration, scenario)
deltaMAP_m = reshape(X_m(42,end,1:end-2) - X_m(42,iteration,1:end-2), [1,num_scen-2]);
deltaMAP_f = reshape(X_f(42,end,1:end-2) - X_f(42,iteration,1:end-2), [1,num_scen-2]);
MAP_comp = deltaMAP_m(1) - [deltaMAP_m(1), deltaMAP_f];
scen_comp = categorical({'M - M'       , 'M - F'       , ...
                         'M - F M RSNA', 'M - F M AT2R', ...
                         'M - F M RAS' , 'M - F M Reab'});
scen_comp = reordercats(scen_comp,{'M - M'       , 'M - F'       , ...
                                   'M - F M RSNA', 'M - F M AT2R', ...
                                   'M - F M RAS' , 'M - F M Reab'});

j = figure('DefaultAxesFontSize',10);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.15, 3.5]);
s1(1) = subplot(1,2,1); 
s1(2) = subplot(1,2,2); 

plot(s1(1), X_m(42,:,fixed_ss),xscale,'-', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3);
% xlim(s1(1), [80, 120]); set(s1(1),'XTick', [80,100,120]);
ylim(s1(1), [lower, upper])
xlabel(s1(1), 'MAP (mmHg)', 'FontSize',14*1.1); ylabel(s1(1), {'Fold change in'; 'sodium excretion'}, 'FontSize',14);
hold(s1(1), 'on')
plot(s1(1), X_f(42,:,fixed_ss),xscale,'-', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3)
legend(s1(1), {'Male','Female'}, 'Location','Northwest', 'FontSize',14)
hold(s1(1), 'off')
title(s1(1), 'A', 'FontSize',14)

bar(s1(2), scen_comp,MAP_comp,'k');
% set(gca,'xticklabel',scen_comp_text);
% xtickget = get(gca,'xticklabel');  
% set(gca,'xticklabel',xtickget,'fontsize',6)
% xtickangle(s1(2),90)
% xlim(s1(2), [1-1,6+1])
ylim(s1(2), [-1,7])
xlabel(s1(2), 'Scenario', 'FontSize',14); ylabel(s1(2), '\DeltaMAP (mmHg)', 'FontSize',14);
% hAxes.XAxis.FontSize = 6;
title(s1(2), 'B', 'FontSize',14)

% Save figures. -----------------------------------------------------------

% if     strcmp(win,  'fixed')
%     save_data_name = sprintf('all_vars_vs_Phisodin_fixed_Phiwin.fig' );
% elseif strcmp(win, 'varied')
%     save_data_name = sprintf('all_vars_vs_Phisodin_varied_Phiwin.fig');
% end
% save_data_name = strcat('Figures/', save_data_name);
% savefig([f;g;h;i;j], save_data_name)

end


























