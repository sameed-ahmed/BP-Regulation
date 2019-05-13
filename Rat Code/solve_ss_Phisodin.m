% This script varies the sodium intake and then calculates the steady state
% solution of the system for each inputted value. To avoid the solver not
% converging, the initial guess for solution to the system is taken as the
% previous solution value. That is, IG_i = SOL_i-1.
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

% Boolean to fix/vary water intake.
% win =  'fixed';
win = 'varied';

% Number of iterations below/above baseline.
iteration = 51;
% Fold decrease/increase.
lower = 1/5; upper = 5;

% Scenarios
% Normal - normal conditions
% ACEi   - Angiotensin convernting enzyme inhibitor
% AngII  - Ang II infusion
scenario = {'Normal', 'ACEi', 'AngII'};
num_scen = length(scenario);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of variables changes based on whether water intake is
% fixed/varied.
if     strcmp(win,  'fixed')
    num_vars = 91-1;
elseif strcmp(win, 'varied')
    num_vars = 91;
end

% Baseline of water intake for each scenario if it is fixed.
Phi_win_bl_m = zeros(num_scen,2*iteration-1);
Phi_win_bl_f = zeros(num_scen,2*iteration-1);
% Load data for baseline water intake and renal perfusion pressure for each
% scenario.
load(  'male_ss_data_scenario_Normal.mat', 'SSdata');
Phi_win_bl_m(1,1) = SSdata(27);
clear SSdata;
% load('female_ss_data_scenario_Normal.mat', 'SSdata');
% Phi_win_bl_f(1,1) = SSdata(27);
% clear SSdata;
load(  'male_ss_data_scenario_AngII.mat', 'SSdata');
Phi_win_bl_m(2,1) = SSdata(27);
clear SSdata;
% load('female_ss_data_scenario_AngII.mat', 'SSdata');
% Phi_win_bl_f(2,1) = SSdata(27);
% clear SSdata;
load(  'male_ss_data_scenario_ACEi.mat', 'SSdata');
Phi_win_bl_m(3,1) = SSdata(27);
clear SSdata;
% load('female_ss_data_scenario_ACEi.mat', 'SSdata');
% Phi_win_bl_f(3,1) = SSdata(27);
% clear SSdata;

% Range for fold decrease/increase.
iter_range_l = linspace(lower, 1, iteration);
iter_range_u = linspace(1, upper, iteration);
% Delete overlapping baseline value.
iter_range_l(end) = '';
% Combine decrease/increase ranges.
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

gender   = {'male',     'female'  };
change   = {'decrease', 'increase'};

for ss = 1:num_scen % scenario
for gg = 1:1        % gender
for cc = 1:2        % change

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

% Load data for steady state initial value. 
if strcmp(scenario{ss}, 'Normal')
    if     strcmp(gender{gg}, 'male')
        load(  'male_ss_data_scenario_Normal.mat', 'SSdata');
    elseif strcmp(gender{gg}, 'female')
        load('female_ss_data_scenario_Normal.mat', 'SSdata');
    end
elseif strcmp(scenario{ss}, 'ACEi')
    if     strcmp(gender{gg}, 'male')
        load(  'male_ss_data_scenario_ACEi.mat', 'SSdata');
    elseif strcmp(gender{gg}, 'female')
        load('female_ss_data_scenario_ACEi.mat', 'SSdata');
    end
elseif strcmp(scenario{ss}, 'AngII')
    if     strcmp(gender{gg}, 'male')
        load(  'male_ss_data_scenario_AngII.mat', 'SSdata');
    elseif strcmp(gender{gg}, 'female')
        load('female_ss_data_scenario_AngII.mat', 'SSdata');
    end
end

% Retrieve and replace parameters in fixed variable equations.
fixed_ind = [2, 10, 14, 20, 24, 43, 48, 61, 65, 70, 87];
fixed_var_pars = SSdata(fixed_ind);
SF = 4.5*10^(-3)*10^(3);
a = 0.2787 * exp(SSdata(32) * 0.2281 / SF);
fixed_var_pars = [fixed_var_pars; a];
SSdata(fixed_ind) = 1;
SSdataIG = SSdata;
clear SSdata

% Delete Phi_win if it is fixed.
if     strcmp(win,  'fixed')
    SSdataIG(27) = '';
elseif strcmp(win, 'varied')
end

for iter = 1:iteration % range

%% Parameters

if     strcmp(gender{gg}, 'male')
    gen = 1;
elseif strcmp(gender{gg}, 'female')
    gen = 0;
end

% Scaling factor
% Rat flow = Human flow x SF
if     strcmp(gender{gg}, 'male')
    SF = 4.5*10^(-3)*10^(3);
elseif strcmp(gender{gg}, 'female')
    SF = 2/3 * 4.5*10^(-3)*10^(3);
end

N_rsna    = 1;
R_aass    = 31.67 / SF;   % mmHg min / l
R_eass    = 51.66 / SF;   % mmHg min / l
P_B       = 18;           % mmHg
P_go      = 28;           % mmHg
C_gcf     = 0.00781 * SF;
if     strcmp(gender{gg}, 'male')
    eta_ptsodreab_eq = 0.93; 
    eta_dtsodreab_eq = 0.77; 
    eta_cdsodreab_eq = 0.15;
elseif strcmp(gender{gg}, 'female')
    eta_ptsodreab_eq = 0.5;
    eta_dtsodreab_eq = 0.5; 
    eta_cdsodreab_eq = 0.972;
end
if     strcmp(gender{gg}, 'male')
    eta_ptwreab_eq = 0.86; 
    eta_dtwreab_eq = 0.60; 
    eta_cdwreab_eq = 0.78;
elseif strcmp(gender{gg}, 'female')
    eta_ptwreab_eq = 0.5;
    eta_dtwreab_eq = 0.5; 
    eta_cdwreab_eq = 0.972;
end
K_vd      = 0.00001;
K_bar     = 16.6 / SF;    % mmHg min / l
R_bv      = 3.4 / SF;     % mmHg min / l
T_adh     = 6;            % min
% Phi_sodin = 1.2278;   % mEq / min
C_K       = 5;               % mEq / l 
T_al      = 30;           % min LISTED AS 30 IN TABLE %listed as 60 in text will only change dN_al
N_rs      = 1;            % ng / ml / min

% Baseline/range of sodium intake.
Phi_sodin_bl_m = 1.2278;
Phi_sodin_bl_f = 1.2278; % CHECK IF DIFFERENT LATER.
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

% RAS
h_renin   = 12;      % min
h_AGT     = 10*60;   % min
h_AngI    = 0.5;     % min
h_AngII   = 0.66;    % min
h_Ang17   = 30;      % min
h_AngIV   = 0.5;     % min
h_AT1R    = 12;      % min
h_AT2R    = 12;      % min

% Male and female different parameters for RAS
if     strcmp(gender{gg}, 'male')
    X_PRCPRA = 135.59/17.312;
    k_AGT    = 801.02;
    c_ACE    = 0.096833;
    c_Chym   = 0.010833;
    c_NEP    = 0.012667;
    c_ACE2   = 0.0026667;
    c_IIIV   = 0.29800;
    c_AT1R   = 0.19700;
    c_AT2R   = 0.065667;
    AT1R_eq  = 20.4807902818665;
    AT2R_eq  = 6.82696474842298;
elseif strcmp(gender{gg}, 'female')
    X_PRCPRA = 114.22/17.312;
    k_AGT    = 779.63;
    c_ACE    = 0.11600;
    c_Chym   = 0.012833;
    c_NEP    = 0.0076667;
    c_ACE2   = 0.00043333;
    c_IIIV   = 0.29800;
    c_AT1R   = 0.19700;
    c_AT2R   = 0.065667;
    AT1R_eq  = 20.4538920068419;
    AT2R_eq  = 6.81799861123497;
end

% Input Phi_win if it is fixed.
if     strcmp(gender{gg}, 'male')
    Phi_win_input = Phi_win_bl_m(ss,1);
elseif strcmp(gender{gg}, 'female')
    Phi_win_input = Phi_win_bl_f(ss,1);
end

pars = [N_rsna; R_aass; R_eass; P_B; P_go; C_gcf; eta_ptsodreab_eq; ...
        eta_dtsodreab_eq; eta_cdsodreab_eq; eta_ptwreab_eq; ...
        eta_dtwreab_eq; eta_cdwreab_eq; K_vd; K_bar; R_bv; T_adh; ...
        Phi_sodin; C_K; T_al; N_rs; X_PRCPRA; h_renin; h_AGT; h_AngI; ...
        h_AngII; h_Ang17; h_AngIV; h_AT1R; h_AT2R; k_AGT; c_ACE; ...
        c_Chym; c_NEP; c_ACE2; c_IIIV; c_AT1R; c_AT2R; AT1R_eq; ...
        AT2R_eq; gen; SF];

%% Drugs

% drugs = [Ang II inf rate fmol/(ml min), ACEi target level]
if     strcmp(scenario{ss}, 'Normal')
    drugs = [0   , 0];
elseif strcmp(scenario{ss}, 'ACEi')
    drugs = [0   , 1]; % Hall 2018
elseif strcmp(scenario{ss}, 'AngII')
    if     strcmp(gender{gg}, 'male')
        drugs = [(3/3)*10984, 0]; % Sampson 2008
    elseif strcmp(gender{gg}, 'female')
        drugs = [(2/3)*10984, 0]; % Sampson 2008
    end
end

%% Variables initial guess

names  = {'$rsna$'; '$\alpha_{map}$'; '$\alpha_{rap}$'; '$R_{r}$'; ...
          '$\beta_{rsna}$'; '$\Phi_{rb}$'; '$\Phi_{gfilt}$'; '$P_{f}$'; ...
          '$P_{gh}$'; '$\Sigma_{tgf}$'; '$\Phi_{filsod}$'; ...
          '$\Phi_{pt-sodreab}$'; '$\eta_{pt-sodreab}$'; ...
          '$\gamma_{filsod}$'; '$\gamma_{at}$'; '$\gamma_{rsna}$'; ...
          '$\Phi_{md-sod}$'; '$\Phi_{dt-sodreab}$'; ...
          '$\eta_{dt-sodreab}$'; '$\psi_{al}$'; '$\Phi_{dt-sod}$'; ...
          '$\Phi_{cd-sodreab}$'; '$\eta_{cd-sodreab}$'; ...
          '$\lambda_{dt}$'; '$\lambda_{anp}$'; '$\Phi_{u-sod}$'; ...
          '$\Phi_{win}$'; '$V_{ecf}$'; '$V_{b}$'; '$P_{mf}$'; ...
          '$\Phi_{vr}$'; '$\Phi_{co}$'; '$P_{ra}$'; '$vas$'; ...
          '$vas_{f}$'; '$vas_{d}$'; '$R_{a}$'; '$R_{ba}$'; '$R_{vr}$'; ...
          '$R_{tp}$'; '$P_{ma}$'; '$\epsilon_{aum}$'; '$a_{auto}$'; ...
          '$a_{chemo}$'; '$a_{baro}$'; '$C_{adh}$'; '$N_{adh}$'; ...
          '$N_{adhs}$'; '$\delta_{ra}$'; '$\Phi_{pt-wreab}$'; ...
          '$\eta_{pt-wreab}$'; '$\mu_{pt-sodreab}$'; '$\Phi_{md-u}$'; ...
          '$\Phi_{dt-wreab}$'; '$\eta_{dt-wreab}$'; ...
          '$\mu_{dt-sodreab}$'; '$\Phi_{dt-u}$'; '$\Phi_{cd-wreab}$'; ...
          '$\eta_{cd-wreab}$'; '$\mu_{cd-sodreab}$'; '$\mu_{adh}$'; ...
          '$\Phi_{u}$'; '$M_{sod}$'; '$C_{sod}$'; '$\nu_{md-sod}$'; ...
          '$\nu_{rsna}$'; '$C_{al}$'; '$N_{al}$'; '$N_{als}$'; ...
          '$\xi_{k/sod}$'; '$\xi_{map}$'; '$\xi_{at}$'; ...
          '$\hat{C}_{anp}$'; '$AGT$'; '$\nu_{AT1}$'; '$R_{sec}$'; ...
          '$PRC$'; '$PRA$'; '$Ang I$'; '$Ang II$'; ...
          '$Ang II_{AT1R-bound}$'; '$Ang II_{AT2R-bound}$'; ...
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
exitflag, output] = fsolve(@(x) bp_reg_solve_Phisodin(t,x,x_p0,pars, ...
                                                      fixed_var_pars, ...
                                                      drugs,win, ...
                                                      Phi_win_input), ...
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
% Update next initial guess as current solution
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
    Phi_win_bl_m(ss,:) = Phi_win_bl_m(ss,1) * ones(1,2*iteration-1);
    Phi_win_bl_f(ss,:) = Phi_win_bl_f(ss,1) * ones(1,2*iteration-1);
    X_m(:,:,ss) = [X(1:26,:,1,ss); Phi_win_bl_m(ss,:); X(27:end,:,1,ss)];
    X_f(:,:,ss) = [X(1:26,:,2,ss); Phi_win_bl_f(ss,:); X(27:end,:,2,ss)];
elseif strcmp(win, 'varied')
    X_m(:,:,ss) = X(:,:,1,ss); X_f(:,:,ss) = X(:,:,2,ss);
end

end % scenario

% x-axis
xscale = iter_range;

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
end
% X_f = zeros(size(X_f));
% X_m = zeros(size(X_m));

f = gobjects(7,1);
s = gobjects(7,15);
% Loop through each set of subplots.
for i = 1:7
    f(i) = figure;
%     f(i) = figure('pos',[750 500 650 450]);
    % This is to avoid the empty plots in the last subplot set.
    if i == 7
        last_plot = 1;
    else
        last_plot = 15;
    end
    % Loop through each subplot within a set of subplots.
    for j = 1:last_plot
        s(i,j) = subplot(3,5,j);

        plot(s(i,j), xscale,X_m((i-1)*15 + j,:,1),'b', ...
                     xscale,X_f((i-1)*15 + j,:,1),'r');
        
        xlim([lower, upper])
        ylim([ylower((i-1)*15 + j), yupper((i-1)*15 + j)])
        
        xlabel('$\Phi_{sodin}$', 'Interpreter','latex', 'FontSize',15)
        title(names((i-1)*15 + j), 'Interpreter','latex', 'FontSize',15)
%         legend('Male', 'Female')
    end
end

% Data
% Relative MAP change, Relative sodium excretion change

pn_m = [1, 1.12; 1, 4.9]; % 
pn_m(1,:) = pn_m(1,:)*101;
pn_f = [1, 1.17; 1, 8.0]; % 
pn_f(1,:) = pn_f(1,:)*100;

% Plot Sodium Intake vs Mean Arterial Pressure

g = figure('pos',[100 100 675 450]);
plot(X_m(41,:,1),xscale,'b-', X_f(41,:,1),xscale,'r-', 'LineWidth',3)
% xlim([90, 120])
ylim([lower, upper])
legend('Male', 'Female')
set(gca,'FontSize',14)
xlabel(names(41)       , 'Interpreter','latex', 'FontSize',22, 'FontWeight','bold')
ylabel('$\Phi_{sodin}$', 'Interpreter','latex', 'FontSize',22, 'FontWeight','bold')
hold all
plot(pn_m(1,:),pn_m(2,:),'bx', pn_f(1,:),pn_f(2,:),'rx', ...
     'MarkerSize',10, 'LineWidth',3)

%

h = figure('pos',[100 100 675 450]);
plot(X_m(41,:,1),xscale,'b-' , 'LineWidth',3, 'DisplayName','M Normal')
% xlim([80, 160])
ylim([lower, upper])
set(gca,'FontSize',14)
xlabel(names(41)       , 'Interpreter','latex', 'FontSize',22, 'FontWeight','bold')
ylabel('$\Phi_{sodin}$', 'Interpreter','latex', 'FontSize',22, 'FontWeight','bold')
legend('-DynamicLegend');
hold all
plot(X_f(41,:,1),xscale,'r-', 'LineWidth',3, 'DisplayName','F Normal')
legend('-DynamicLegend');

hold all
plot(X_m(41,:,2),xscale,'b--' , 'LineWidth',3, 'DisplayName','M ACEi')
legend('-DynamicLegend');
hold all
plot(X_f(41,:,2),xscale,'r--', 'LineWidth',3, 'DisplayName','F ACEi')
legend('-DynamicLegend');

hold all
plot(X_m(41,:,3),xscale,'b:' , 'LineWidth',3, 'DisplayName','M AngII')
legend('-DynamicLegend');
hold all
plot(X_f(41,:,3),xscale,'r:', 'LineWidth',3, 'DisplayName','F AngII')
legend('-DynamicLegend');

% % Save figures.
% 
% % if     strcmp(win,  'fixed')
% %     savefig(f, 'all_vars_vs_Phisodin_fixed_Phiwin.fig' )
% %     savefig(g, 'Phisodin_vs_Pma_fixed_Phiwin.fig'      )
% % elseif strcmp(win, 'varied')
% %     savefig(f, 'all_vars_vs_Phisodin_varied_Phiwin.fig')
% %     savefig(g, 'Phisodin_vs_Pma_varied_Phiwin.fig'     )
% % end
% 
% if     strcmp(win,  'fixed')
%     savefig(f, 'all_vars_vs_Phisodin_fixed_Phiwin_new_Phitwreab.fig'  )
%     savefig(g, 'Phisodin_vs_Pma_fixed_Phiwin_new_Phitwreab.fig'       )
%     savefig(h, 'Phisodin_vs_Pma_fixed_Phiwin_new_Phitwreab_AngII.fig' )
% elseif strcmp(win, 'varied')
%     savefig(f, 'all_vars_vs_Phisodin_varied_Phiwin_new_Phitwreab.fig' )
%     savefig(g, 'Phisodin_vs_Pma_varied_Phiwin_new_Phitwreab.fig'      )
%     savefig(h, 'Phisodin_vs_Pma_varied_Phiwin_new_Phitwreab_AngII.fig')
% end

end


























