% This simulates the blood pressure regulation model bp_reg_mod.m for 
% various scenarios for a given virtual individual.
% 
% Parameters and steady state data for the virtual population are 
% calculated by create_par_bs_rep.m and create_vp.m, respectively.

% Input
% sim_scenario: scenario to simulate
% sample_num  : index of virtual individual
% Output
% plots and saves figures corresponding to the scenario.

function run_sim_hyp_scen

close all

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

%% Variable names for plotting.
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
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Begin user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation scenarios
% Baseline - no experiment
% AngII    - Ang II infusion
% Sodin    - sodium loading
% RPP      - manipulate renal perfusion pressure
sim_scenario = {'Baseline', 'AngII', 'Sodin', 'RPP'};
exact_sim_scen = 1;

% Species
spe_ind = 2;

% Bootstrap replicate sample number
sample_num = random('Discrete Uniform',1000)
% sample_num = 119

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

species = {'human', 'rat'   };
sex     = {'male' , 'female'};

PARS = cell(1,2); SSDATA = cell(1,2); 

for sex_ind = 1:2 % sex

%% Load bootstrap replicate parameters & variables created by create_par_bs_rep.m.

% Parameters
load_data_name_pars = sprintf('%s_%s_pars_scenario_Pri_Hyp_bs_rep1000.mat', ...
                              species{spe_ind},sex{sex_ind});
load(load_data_name_pars, 'pars_rep');
num_pars   = size(pars_rep,1);
num_sample = size(pars_rep,2);
PARS{sex_ind} = pars_rep;

% pars0_est = [24.5664; 11.9122; 3.2875; 1.9027; 1.9448; 1.4909; 1.4893; 4.8474]; % diverges for j = 076
% par_ind   = [13     ; 14     ; 4     ; 21    ; 18    ; 3     ; 15    ; 41    ];
% PARS{sex_ind}(par_ind,sample_num) = pars0_est;

% Variables
load_data_name_vars = sprintf('%s_%s_ss_data_scenario_Pri_Hyp_bs_rep1000.mat', ...
                              species{spe_ind},sex{sex_ind});
load(load_data_name_vars, 'SSdata_rep');
num_vars   = size(SSdata_rep,1);
SSDATA{sex_ind} = SSdata_rep;

end % sex

%% Run some simulation.

if     strcmp(sim_scenario{exact_sim_scen}, 'Baseline')
    fig = run_sim(PARS,SSDATA,sample_num);
elseif strcmp(sim_scenario{exact_sim_scen}, 'AngII')
    fig = run_sim_AngII(PARS,SSDATA,sample_num);
elseif strcmp(sim_scenario{exact_sim_scen}, 'Sodin')
    fig = solve_ss_Phisodin(PARS,SSDATA,sample_num);
elseif strcmp(sim_scenario{exact_sim_scen}, 'RPP')
    fig = run_sim_RPP(PARS,SSDATA,sample_num);
end

%% Save figures.

save_data_name = sprintf('Pri_hyp_sim_%s_VI%s.fig', ...
                         sim_scenario{exact_sim_scen},num2str(sample_num));
save_data_name = strcat('Figures/', save_data_name);
savefig(fig, save_data_name)
% ---
save_data_name = sprintf('Pri_hyp_sim_%s_VI%s.png', ...
                         sim_scenario{exact_sim_scen},num2str(sample_num));
save_data_name = strcat('Figures/', save_data_name);
exportgraphics(fig(end), save_data_name)

%% ---------------------------------
%% Subfunctions
%% ---------------------------------

% -------------------------------------------------------------------------
% Steady state simulation
% -------------------------------------------------------------------------

function f = run_sim(PARS,SSDATA,sample_num)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Begin user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scenario
scenario = {'Normal'};
fixed_ss = 1;

% Number of days to run simulation after change; Day at which to induce change;
days = 10; day_change = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize variables and time.
X = cell(1,2); T = cell(1,2);

for sex_ind_run_sim = 1:2 % sex
% for sex_ind_run_sim = fixed_sex_ind:fixed_sex_ind % sex 

varargin_input = {scenario{fixed_ss},true};

% Retreive parameters and steady state data.
pars = PARS{sex_ind_run_sim}; SSdata = SSDATA{sex_ind_run_sim}; 
% % % Note: SSdata is cleared above. Therefore it is not a scoping variable. % % % 

pars = pars(:,sample_num); SSdata = SSdata(:,sample_num);

%% Solve DAE

% Initial condition for the variables and their derivatives. 
% System is initially at steady state, so the derivative is 0.
x0 = SSdata; x_p0 = zeros(num_vars,1);

% Time at which to keep steady state, change a parameter, etc.
tchange = day_change*1440;

% Initial time (min); Final time (min);
t0 = 0*1440; tend = tchange + days*1440;

% Time vector
tspan = [t0, tend];

% ode options
options = odeset('MaxStep',100); % default is 0.1*abs(t0-tf)

% Solve dae
[t,x] = ode15i(@(t,x,x_p) ...
               bp_reg_mod(t,x,x_p,pars,tchange,varargin_input{:}), ...
               tspan, x0, x_p0, options);

T{sex_ind_run_sim} = t';
X{sex_ind_run_sim} = x';

end % sex

%% Plot

% Retrieve male and female.
t_m = T{1}; X_m = X{1}; 
t_f = T{2}; X_f = X{2}; 

% x-axis limits
xlower = t0; xupper = tend; 

% Convert minutes to days for longer simulations.
t_m = t_m/1440; t_f = t_f/1440; tchange = tchange/1440; 
xlower = xlower/1440; xupper = xupper/1440; 

% y-axis limits
ylower = zeros(length(X_m(:,1)),1); yupper = ylower; 
for i = 1:length(ylower)
    ylower(i) = 0.95*min( min(X_m(i,:)), min(X_f(i,:)) );
    yupper(i) = 1.05*max( max(X_m(i,:)), max(X_f(i,:)) );
    if abs(yupper(i)) < eps*100
        ylower(i) = -10^(-5); yupper(i) = 10^(-5);
    end
end

f1 = gobjects(7,1);
s1 = gobjects(7,15);
% Loop through each set of subplots.
for i = 1:7
%     f(i) = figure();
    f1(i) = figure('pos',[750 500 650 450]);
    % This is to avoid the empty plots in the last subplot set.
    if i == 7
        last_plot = mod(num_vars, 15);
    else
        last_plot = 15;
    end
    % Loop through each subplot within a set of subplots.
    for j = 1:last_plot
        s1(i,j) = subplot(3,5,j);
        s1(i,j).Position = s1(i,j).Position + [0 0 0.01 0];
        
        plot(s1(i,j), t_m,X_m((i-1)*15 + j,:),'b', t_f,X_f((i-1)*15 + j,:),'r');
        
        xlim([xlower, xupper])
        ylim([ylower((i-1)*15 + j), yupper((i-1)*15 + j)])
        
%         Days
        ax = gca;
        ax.XTick = (tchange+0*(1) : 1 : tchange+days*(1));
        ax.XTickLabel = {'0' ,'1' ,'2' ,'3' ,'4' ,'5' ,'6' ,...
                         '7' ,'8' ,'9' ,'10','11','12','13',...
                         '14','15','16','17','18','19','20',...
                         '21','22','23','24','25','26'};
        
%         legend('Male', 'Female')
        title(names((i-1)*15 + j), 'Interpreter','latex', 'FontSize',15)
    end
end

f = [f1];

end % run_sim

% -------------------------------------------------------------------------
% Ang II infusion
% -------------------------------------------------------------------------

function f = run_sim_AngII(PARS,SSDATA,sample_num)

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
% m_RSNA_&_m_Reab - male RSNA & fractional sodium and water reabsorption
scenario = {'Normal', 'm_RSNA', 'm_AT2R', 'm_RAS', 'm_Reab', ...
            'm_RAS_m_Reab', 'm_RSNA_m_Reab'};
% scenario = {'Normal', 'm_RSNA', 'm_AT2R', 'm_RAS', 'm_Reab', ...
%             'm_RAS_m_Reab', 'm_RSNA_m_Reab', ...
%             'Pri_Hyp'};
num_scen = length(scenario);
fixed_ss = 1;

% Number of days to run simulation after change; Day at which to induce change;
days = 14; day_change = 1;
% Number of points for plotting resolution
% N = ((days+1)*1440) / 2;
N = (days+1)*10 + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize variables.
% X = (variables, points, sex, scenario)
X = zeros(num_vars,N,2,num_scen);

for sce_ind               = 1:num_scen % scenario
for sex_ind_AngII = 1:2        % sex
% for sex_ind_AngII = fixed_sex_ind:fixed_sex_ind        % sex

varargin_input = {scenario{sce_ind},true};

% Retreive parameters and steady state data.
pars = PARS{sex_ind_AngII}; SSdata = SSDATA{sex_ind_AngII}; 
% % % Note: SSdata is cleared above. Therefore it is not a scoping variable. % % % 

pars = pars(:,sample_num); SSdata = SSdata(:,sample_num);

%% Drugs

% Ang II inf rate fmol/(ml min)
if     strcmp(sex{sex_ind_AngII}, 'male')
%     kappa_AngII = 2022; % Sampson 2008
    kappa_AngII = 910; % Sullivan 2010
%     kappa_AngII = 630; % Sullivan 2010
%     kappa_AngII = 0;
elseif strcmp(sex{sex_ind_AngII}, 'female')
%     kappa_AngII = 2060; % Sampson 2008
    kappa_AngII = 505; % Sullivan 2010
%     kappa_AngII = 630; % Sullivan 2010
%     kappa_AngII = 0;
end

varargin_input = [varargin_input, 'AngII',kappa_AngII];

%% Solve DAE

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
options = odeset('MaxStep',1000); % default MaxStep is 0.1*abs(t0-tf)

% Solve dae
[t,x] = ode15i(@(t,x,x_p) ...
               bp_reg_mod(t,x,x_p,pars,tchange,varargin_input{:}), ...
               tspan, x0, x_p0, options);

% X = (variables, points, sex, scenario)
X(:,:,sex_ind_AngII,sce_ind) = x';

end % sex
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

% Interesting variables to plot.
var_ind = [33;41;42;9;73;74;6;7;27;92;93;29]; sub_var_num = length(var_ind);

% Plot all vars vs time. --------------------------------------------------

f1 = gobjects(7,1);
s1 = gobjects(7,15);
% Loop through each set of subplots.
for i = 1:7
    f1(i) = figure('pos',[750 500 650 450]);
    % This is to avoid the empty plots in the last subplot set.
    if i == 7
        last_plot = mod(num_vars, 15);
    else
        last_plot = 15;
    end
    % Loop through each subplot within a set of subplots.
    for j = 1:last_plot
        s1(i,j) = subplot(3,5,j);
        s1(i,j).Position = s1(i,j).Position + [0 0 0.01 0];
        
        plot(s1(i,j), t,X_m((i-1)*15 + j,:,fixed_ss),'b', ...
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

% h(1) = figure('DefaultAxesFontSize',14);
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 2.5]);
% plot(t,RSNA_m(:,fixed_ss),'-', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3);
% xlim([xlower, xupper]); ylim([0.50,1.05]);
% ax = gca;
% ax.XTick = (tchange+0*(1) : 2 : tchange+days*(1));
% ax.XTickLabel = {'0','2','4','6','8','10','12','14'};
% xlabel('Time (days)'); ylabel('RSNA');
% hold on
% plot(t,RSNA_f(:,fixed_ss),'-', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3)
% [~, hobj, ~, ~] = legend({'Male','Female'}, 'FontSize',7,'Location','Northeast');
% hl = findobj(hobj,'type','line');
% set(hl,'LineWidth',1.5);
% % ---
f3 = figure('DefaultAxesFontSize',14);
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

% Plot Mean Arterial Pressure vs Time. ------------------------------------

% Data from Sullivan 2010. MAP is in difference from baseline.
load_data_name_AngII = sprintf('%s_male_AngII_data_bs_rep.mat'  , species{spe_ind});
load(load_data_name_AngII, 'AngII_data_rep');
MAPdata_m = AngII_data_rep(sample_num,:);
load_data_name_AngII = sprintf('%s_female_AngII_data_bs_rep.mat', species{spe_ind});
load(load_data_name_AngII, 'AngII_data_rep');
MAPdata_f = AngII_data_rep(sample_num,:);

tdata = [0+1 , 1+1 , 2+1 , 3+1 , 4+1 , 5+1 , 6+1 ,...
         7+1 , 8+1 , 9+1 , 10+1, 11+1, 12+1, 13+1, 14+1];

% Substract MAP by baseline for each sex and all scenarios.
% X_m/f = (variable, points, scenario)
MAP_m = zeros(N,num_scen); MAP_f = zeros(N,num_scen);
for i = 1:num_scen
    MAP_m(:,i) = X_m(42,:,i) - X_m(42,1,i);
    MAP_f(:,i) = X_f(42,:,i) - X_f(42,1,i);
end
% MAP_m = reshape(X_m(42,:,i) - X_m(42,1,i), [N,num_scen]);
% MAP_f = reshape(X_f(42,:,i) - X_f(42,1,i), [N,num_scen]);

f4 = figure('DefaultAxesFontSize',16*1.5);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7, 5]);
plot(t,MAP_m(:,fixed_ss),'-', 'Color',[0.203, 0.592, 0.835], 'LineWidth',4.5);
xlim([xlower, xupper]); ylim([-1, 60]);
ax = gca;
ax.XTick = (tchange+0*(1) : 2 : tchange+days*(1));
ax.XTickLabel = {'0','2','4','6','8','10','12','14'};
xlabel('Time (days)'); ylabel('\DeltaMAP (mmHg)');
hold on
plot(t,MAP_f(:,fixed_ss),'-', 'Color',[0.835, 0.203, 0.576], 'LineWidth',4.5)
plot(tdata,MAPdata_m,'o', 'DisplayName',  'Male data', 'Color',[0.203, 0.592, 0.835], 'MarkerSize',9, 'LineWidth',3)
plot(tdata,MAPdata_f,'o', 'DisplayName','Female data', 'Color',[0.835, 0.203, 0.576], 'MarkerSize',9, 'LineWidth',3)
[~, hobj, ~, ~] = legend({'Male sim','Female sim','Male data','Female data'}, 'FontSize',8*1.5,'Location','Northwest');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',2.25);
title('A', 'FontWeight','normal')

% % Plot male - female bar graph for each scenario --------------------------
% 
% % Retrieve last MAP value for each sex and all scenarios.
% % X_m/f = (variable, points, scenario)
% deltaMAP_m = reshape(X_m(42,end,1:end) - X_m(42,1,1:end), [1,num_scen]);
% deltaMAP_f = reshape(X_f(42,end,1:end) - X_f(42,1,1:end), [1,num_scen]);
% % Substract all MAP values from all scenarios from normal male MAP.
% MAP_comp = deltaMAP_m(1) - [deltaMAP_m(1), deltaMAP_f];
% % String for bar graph labels.
% scen_comp = categorical({'M - M       ', 'M - F       ', ...
%                          'M - F M RSNA', 'M - F M AT2R', ...
%                          'M - F M RAS' , 'M - F M Reab', ...
%                          'M - F M RAS \newline& Reab'  , ...
%                          'M - F M RSNA\newline& Reab'  });
% % Horizontal bar graph
% scen_comp = reordercats(scen_comp,{'M - F M RSNA\newline& Reab'  , ...
%                                    'M - F M RAS \newline& Reab'  , ...
%                                    'M - F M Reab', 'M - F M RAS ', ...
%                                    'M - F M AT2R', 'M - F M RSNA', ...
%                                    'M - F       ', 'M - M       '});
% 
% k = figure('DefaultAxesFontSize',10);
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 3.5]);
% barh(scen_comp,MAP_comp,'k');
% % set(gca,'yticklabel',scen_comp_text);
% % ytickget = get(gca,'yticklabel');  
% % set(gca,'yticklabel',ytickget,'fontsize',6)
% % ylim([1-1,6+1])
% xlim([0,20])
% ylabel('Scenario', 'FontSize',14); xlabel('\DeltaMAP (mmHg)', 'FontSize',14);
% ytickangle(00)
% % hAxes.YAxis.FontSize = 6;

% % % BAR GRAPH IS PRODUCING NO CHANG BETWEEN M AND F SCENARIOS BECAUSE
% PARS ARE NOT BEING RETREIVED FROM 'get_pars.m'. % % %

f = [f1;f2;f3;f4];

end % run_sim_AngII

% -------------------------------------------------------------------------
% Sodium intake
% -------------------------------------------------------------------------

function f = solve_ss_Phisodin(PARS,SSDATA,sample_num)

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
scenario = {'Normal', 'm_RSNA', 'm_AT2R', 'm_RAS', 'm_Reab', ...
            'ACEi', 'AngII'};
% scenario = {'Normal', 'm_RSNA', 'm_AT2R', 'm_RAS', 'm_Reab', ...
%             'ACEi', 'AngII', ...
%             'Pri_Hyp'};
num_scen = length(scenario);
% Index of scenario to plot for all variables
fixed_ss = 1;

% Number of iterations below/above baseline.
iteration = 11; % must be odd number for symmetry
% Fold decrease/increase.
lower = 1/4; upper = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Range for fold decrease/increase.
iter_range_l = linspace(lower, 1, iteration);
iter_range_u = linspace(1, upper, iteration);
% Delete overlapping baseline value.
iter_range_l(end) = '';
% Combine decrease/increase ranges into single vector.
iter_range   = [iter_range_l, iter_range_u];

% Initialize variables.
% X = (variable, iteration, sex, scenario)
X = zeros(num_vars,2*iteration-1,2,num_scen);
% Retrieve male/female. 
% X_m/f = (variable, iteration, scenario)
X_m = zeros(num_vars,2*iteration-1,num_scen);
X_f = zeros(num_vars,2*iteration-1,num_scen);

change  = {'decrease', 'increase'};

for sce_ind = fixed_ss:fixed_ss % scenario
for sex_ind_Phisodin = 1:2        % sex
for cha_ind = 1:2        % change

varargin_input = {scenario{sce_ind},true};

% Retreive parameters and steady state data.
pars = PARS{sex_ind_Phisodin}; SSdata = SSDATA{sex_ind_Phisodin}; 
% % % Note: SSdata is cleared above. Therefore it is not a scoping variable. % % % 

pars = pars(:,sample_num); SSdata = SSdata(:,sample_num);
SSdataIG = SSdata; clear SSdata;

% Baseline/range of sodium intake.
Phi_sodin_bl = pars(17);
Phi_sodin_range = Phi_sodin_bl * iter_range;

for iter = 1:iteration % range

%% Parameters

% Vary sodium intake.
if     strcmp(change{cha_ind}, 'decrease')
    Phi_sodin = Phi_sodin_range(iteration-iter+1);
elseif strcmp(change{cha_ind}, 'increase')
    Phi_sodin = Phi_sodin_range(iteration+iter-1);
end
pars(17) = Phi_sodin;

%% Drugs

% drugs = [ACEi target level, Ang II inf rate fmol/(ml min)]
if     strcmp(scenario{sce_ind}, 'ACEi' )
        varargin_input = {'ACEi' ,1   }; % Hall 1980
elseif strcmp(scenario{sce_ind}, 'AngII')
    if     strcmp(sex{sex_ind_Phisodin}, 'male'  )
        varargin_input = {'AngII',2022}; % Sampson 2008
    elseif strcmp(sex{sex_ind_Phisodin}, 'female')
        varargin_input = {'AngII',2060}; % Sampson 2008
    end
end

%% Variables initial guess

% Initial guess for the variables.
% Find the steady state solution, so the derivative is 0.
% Arbitrary value for time to input, greater than tchange + deltat.
x0 = SSdataIG; x_p0 = zeros(num_vars,1); t = 30;

% Time at which to change place holder.
tchange = 0;

%% Find steady state solution

options = optimset('Display','off');
[SSdata, residual, ...
 exitflag, output] = fsolve(@(x) ...
                            bp_reg_mod(t,x,x_p0,pars,tchange,varargin_input{:}), ...
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
if     strcmp(change{cha_ind}, 'decrease')
    X(:,iteration-iter+1,sex_ind_Phisodin,sce_ind) = SSdata;
elseif strcmp(change{cha_ind}, 'increase')
    X(:,iteration+iter-1,sex_ind_Phisodin,sce_ind) = SSdata;
end

% To avoid the solver not converging, the initial guess for solution to the
% system is taken as the previous solution value. That is, IG_i = SOL_i-1.
% Update next initial guess as current solution.
SSdataIG = SSdata;

% Sanity check to see script's progress. Also a check for where to
% troubleshoot in case the solver does not converge.
fprintf('%s %s %s iteration = %s out of %s \n', ...
        scenario{sce_ind},sex{sex_ind_Phisodin},change{cha_ind},num2str(iter),num2str(iteration))

end % range
end % change
end % sex

% Retrieve data
X_m(:,:,sce_ind) = X(:,:,1,sce_ind); X_f(:,:,sce_ind) = X(:,:,2,sce_ind);

end % scenario

%% Plot

% x-axis
xscale = iter_range;

% y-axis limits
ylower = zeros(length(X_m(:,1,fixed_ss)),1); yupper = ylower; 
for i = 1:length(ylower)
    ylower(i) = 0.9*min(min(X_m(i,:,fixed_ss)), min(X_f(i,:,fixed_ss)));
    yupper(i) = 1.1*max(max(X_m(i,:,fixed_ss)), max(X_f(i,:,fixed_ss)));
    if ylower(i) == yupper(i)
        ylower(i) = ylower(i) - 10^(-5); yupper(i) = yupper(i) + 10^(-5);
    end
    if ylower(i) == inf || yupper(i) == inf || isnan(ylower(i)) || isnan(yupper(i))
        ylower(i) = -1; yupper(i) = 1;
    end
end

% Interesting variables to plot.
var_ind = [33;41;42;9;73;74;6;7;27;92;93;29]; sub_var_num = length(var_ind);

% Plot all variables vs sodium intake. ------------------------------------

f1 = gobjects(7,1);
s1 = gobjects(7,15);
% Loop through each set of subplots.
for i = 1:7
%     f(i) = figure();
    f1(i) = figure('pos',[750 500 650 450]);
    % This is to avoid the empty plots in the last subplot set.
    if i == 7
        last_plot = mod(num_vars, 15);
    else
        last_plot = 15;
    end
    % Loop through each subplot within a set of subplots.
    for j = 1:last_plot
        s1(i,j) = subplot(3,5,j);

        plot(s1(i,j), xscale,X_m((i-1)*15 + j,:,fixed_ss),'b', ...
                     xscale,X_f((i-1)*15 + j,:,fixed_ss),'r');
        
        xlim([lower, upper])
        ylim([ylower((i-1)*15 + j), yupper((i-1)*15 + j)])
        
        xlabel('$\Phi_{sodin}$', 'Interpreter','latex', 'FontSize',15)
        title(names((i-1)*15 + j), 'Interpreter','latex', 'FontSize',15)
%         legend('Male', 'Female')
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

    plot(s2(j), xscale,X_m(var_ind(j),:,fixed_ss), 'Color',[0.203, 0.592, 0.835], 'LineWidth',2.5);
    hold(s2(j), 'on')
    plot(s2(j), xscale,X_f(var_ind(j),:,fixed_ss), 'Color',[0.835, 0.203, 0.576], 'LineWidth',2.5);
    hold(s2(j), 'off')

    xlim([lower, upper])
    ylim([ylower(var_ind(j)), yupper(var_ind(j))])

    ylabel(names(var_ind(j)), 'Interpreter','latex', 'FontSize',16)
end
legend(s2(1),'Male','Female', 'Location','east')
xlh = xlabel(s2(11),'$\Phi_{sodin}$', 'Interpreter','latex', 'FontSize',16);
xlh.Position(2) = xlh.Position(2) - 0.005;

% Plot all other quantities of interest. ----------------------------------

% CSOD; CADH; BV; for each sex and all scenarios.
% X_m/f = (variable, iteration, scenario)
CSOD_m = reshape(X_m(52,:,:), [2*iteration-1,num_scen]);
CSOD_f = reshape(X_f(52,:,:), [2*iteration-1,num_scen]);
CADH_m = reshape(X_m(47,:,:), [2*iteration-1,num_scen]);
CADH_f = reshape(X_f(47,:,:), [2*iteration-1,num_scen]);
BV_m   = reshape(X_m(30,:,:), [2*iteration-1,num_scen]);
BV_f   = reshape(X_f(30,:,:), [2*iteration-1,num_scen]);
% Plot as relative change in order to compare male and female.
CSOD_m_bl = CSOD_m(iteration,:);
CSOD_f_bl = CSOD_f(iteration,:);
CADH_m_bl = CADH_m(iteration,:);
CADH_f_bl = CADH_f(iteration,:);
BV_m_bl   = BV_m  (iteration,:);
BV_f_bl   = BV_f  (iteration,:);
for i = 1:2*iteration-1
    CSOD_m(i,:) = CSOD_m(i,:) ./ CSOD_m_bl;
    CSOD_f(i,:) = CSOD_f(i,:) ./ CSOD_f_bl;
    CADH_m(i,:) = CADH_m(i,:) ./ CADH_m_bl;
    CADH_f(i,:) = CADH_f(i,:) ./ CADH_f_bl;
    BV_m  (i,:) = BV_m  (i,:) ./ BV_m_bl  ;
    BV_f  (i,:) = BV_f  (i,:) ./ BV_f_bl  ;
end

% Filtration fraction for sodium and urine for each sex and all scenarios.
FRNA_m = reshape((X_m(11,:,:) - X_m(27,:,:)) ./ X_m(11,:,:), [2*iteration-1,num_scen]) * 100;
FRNA_f = reshape((X_f(11,:,:) - X_f(27,:,:)) ./ X_f(11,:,:), [2*iteration-1,num_scen]) * 100;
FRW_m  = reshape((X_m( 7,:,:) - X_m(92,:,:)) ./ X_m( 7,:,:), [2*iteration-1,num_scen]) * 100;
FRW_f  = reshape((X_f( 7,:,:) - X_f(92,:,:)) ./ X_f( 7,:,:), [2*iteration-1,num_scen]) * 100;
% Plot as relative change in order to compare male and female.
FRNA_m_bl = FRNA_m(iteration,:);
FRNA_f_bl = FRNA_f(iteration,:);
FRW_m_bl  = FRW_m (iteration,:);
FRW_f_bl  = FRW_f (iteration,:);
for i = 1:2*iteration-1
    FRNA_m(i,:) = FRNA_m(i,:) ./ FRNA_m_bl;
    FRNA_f(i,:) = FRNA_f(i,:) ./ FRNA_f_bl;
    FRW_m(i,:)  = FRW_m(i,:)  ./ FRW_m_bl ;
    FRW_f(i,:)  = FRW_f(i,:)  ./ FRW_f_bl ;
end

f3 = figure('DefaultAxesFontSize',14);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.15, 5]);
s_main(1) = subplot(2,2,1); 
s_main(2) = subplot(2,2,2); 
s_main(3) = subplot(2,2,3);
s_main(4) = subplot(2,2,4); 

plot(s_main(1), xscale,CSOD_m(:,fixed_ss), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim(s_main(1), [lower, upper]);
set(s_main(1), 'XTick', [1/5, 1, 2, 3, 4, 5]);
set(s_main(1), 'XTickLabel', {'^{1}/_{5}','1','2','3','4','5'});
% ylim(s_main(1), [0.99,1.03])
xlabel(s_main(1), 'Na^+ Intake (relative)'); ylabel(s_main(1), 'C_{Na^+} (relative)');
hold(s_main(1), 'on')
plot(s_main(1), xscale,CSOD_f(:,fixed_ss), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold(s_main(1), 'off')
[~, hobj, ~, ~] = legend(s_main(1), {'Male','Female'}, 'FontSize',7,'Location','Southeast');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
title(s_main(1), 'A')

plot(s_main(2), xscale,CADH_m(:,fixed_ss), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim(s_main(2), [lower, upper]);
set(s_main(2), 'XTick', [1/5, 1, 2, 3, 4, 5]);
set(s_main(2), 'XTickLabel', {'^{1}/_{5}','1','2','3','4','5'});
% ylim(s_main(2), [0.6,2.2])
xlabel(s_main(2), 'Na^+ Intake (relative)'); ylabel(s_main(2), 'C_{ADH} (relative)');
hold(s_main(2), 'on')
plot(s_main(2), xscale,CADH_f(:,fixed_ss), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold(s_main(2), 'off')
title(s_main(2), 'B')

plot(s_main(3), xscale,BV_m  (:,fixed_ss), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim(s_main(3), [lower, upper]);
set(s_main(3), 'XTick', [1/5, 1, 2, 3, 4, 5]);
set(s_main(3), 'XTickLabel', {'^{1}/_{5}','1','2','3','4','5'});
% ylim(s_main(3), [0.985,1.015])
xlabel(s_main(3), 'Na^+ Intake (relative)'); ylabel(s_main(3), 'BV (relative)');
hold(s_main(3), 'on')
plot(s_main(3), xscale,BV_f  (:,fixed_ss), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold(s_main(3), 'off')
title(s_main(3), 'C')

plot(s_main(4), xscale,FRNA_m(:,fixed_ss) ,'-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim(s_main(4), [lower, upper]);
set(s_main(4), 'XTick', [1/5, 1, 2, 3, 4, 5]);
set(s_main(4), 'XTickLabel', {'^{1}/_{5}','1','2','3','4','5'});
% ylim(s_main(4), [0.95,1.02])
xlabel(s_main(4), 'Na^+ Intake (relative)'); ylabel(s_main(4), 'FR (relative)');
hold(s_main(4), 'on')
plot(s_main(4), xscale,FRW_m (:,fixed_ss), '--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(s_main(4), xscale,FRNA_f(:,fixed_ss), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(s_main(4), xscale,FRW_f (:,fixed_ss), '--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
fakeplot = zeros(2, 1);
fakeplot(1) = plot(s_main(4), NaN,NaN, 'k-' );
fakeplot(2) = plot(s_main(4), NaN,NaN, 'k--');
[~, hobj, ~, ~] = legend(fakeplot, {'FR_{Na^+}','FR_{U}'}, 'FontSize',7,'Location','Southwest');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
hold(s_main(4), 'off')
title(s_main(4), 'D')

% Plot Sodium Intake vs Mean Arterial Pressure. ---------------------------

MAP_bl_m = X_m(42,iteration,fixed_ss);
MAP_bl_f = X_f(42,iteration,fixed_ss);
point_m = [4,15+MAP_bl_m; 4,25+MAP_bl_m];
point_f = [4, 5+MAP_bl_f; 4,10+MAP_bl_f];

f4 = figure('DefaultAxesFontSize',16*1.5);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5*2, 2.5*2]);
plot(xscale,X_m(42,:,fixed_ss),'-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3*1.5);
ylim([125, 175])
xlim([lower, upper+0.1])
ax = gca;
ax.YTick = (130 : 10 : 170);
ylabel('MAP (mmHg)')
% xlabel({'Fold change in'; 'sodium excretion'})
xlabel({'Fold change in sodium excretion'})
hold on
plot(point_m(:,1),point_m(:,2),'o-', 'Color',[0.203, 0.592, 0.835], 'LineWidth',2*1.5, 'MarkerSize',6*1.5)
plot(xscale,X_f(42,:,fixed_ss),'-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3*1.5)
plot(point_f(:,1),point_f(:,2),'o-', 'Color',[0.835, 0.203, 0.576], 'LineWidth',2*1.5, 'MarkerSize',6*1.5)
legend('Male sim','Female sim','Male data','Female data', 'Location','Northwest')
[~, hobj, ~, ~] = legend({'Male sim','Male data','Female sim','Female data'}, 'FontSize',8*1.5,'Location','Northwest');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5*1.5);
hold off
title('B', 'FontWeight','normal')

% % Plot with different scenarios. ------------------------------------------
% 
% f5 = figure('pos',[100 100 675 450]);
% plot(X_m(42,:,1),xscale,'b-' , 'LineWidth',3, 'DisplayName','M Normal')
% % xlim([80, 160])
% ylim([lower, upper])
% set(gca,'FontSize',14)
% xlabel(names(42)       , 'Interpreter','latex', 'FontSize',22, 'FontWeight','bold')
% ylabel('$\Phi_{sodin}$', 'Interpreter','latex', 'FontSize',22, 'FontWeight','bold')
% legend('-DynamicLegend');
% hold all
% plot(X_f(42,:,1),xscale,'r-', 'LineWidth',3, 'DisplayName','F Normal')
% legend('-DynamicLegend');
% 
% scen_linestyle_m = {'b-x', 'b-o', 'b-+', 'b-*',};
% scen_linestyle_f = {'r-x', 'r-o', 'r-+', 'r-*',};
% scen_legend = {'RSNA', 'AT2R', 'RAS', 'Reab'};
% for sce_ind = 2:num_scen-2
%     hold all
%     plot(X_m(42,:,sce_ind),xscale,scen_linestyle_m{sce_ind-1}, 'LineWidth',3, 'DisplayName',scen_legend{sce_ind-1})
%     legend('-DynamicLegend');
%     hold all
%     plot(X_f(42,:,sce_ind),xscale,scen_linestyle_f{sce_ind-1}, 'LineWidth',3, 'DisplayName',scen_legend{sce_ind-1})
%     legend('-DynamicLegend');
% end
% 
% f6 = figure('pos',[100 100 675 450]);
% plot(X_m(42,:,1),xscale,'b-' , 'LineWidth',3, 'DisplayName','M Normal')
% % xlim([80, 160])
% ylim([lower, upper])
% set(gca,'FontSize',14)
% xlabel(names(42)       , 'Interpreter','latex', 'FontSize',22, 'FontWeight','bold')
% ylabel('$\Phi_{sodin}$', 'Interpreter','latex', 'FontSize',22, 'FontWeight','bold')
% legend('-DynamicLegend');
% hold all
% plot(X_f(42,:,1),xscale,'r-', 'LineWidth',3, 'DisplayName','F Normal')
% legend('-DynamicLegend');
% 
% scen_linestyle_m = {'b--', 'b:'};
% scen_linestyle_f = {'r--', 'r:'};
% scen_legend = {'ACEi', 'Ang II'};
% for sce_ind = num_scen-1:num_scen
%     hold all
%     plot(X_m(42,:,sce_ind),xscale,scen_linestyle_m{sce_ind-(num_scen-2)}, 'LineWidth',3, 'DisplayName',scen_legend{sce_ind-(num_scen-2)})
%     legend('-DynamicLegend');
%     hold all
%     plot(X_f(42,:,sce_ind),xscale,scen_linestyle_f{sce_ind-(num_scen-2)}, 'LineWidth',3, 'DisplayName',scen_legend{sce_ind-(num_scen-2)})
%     legend('-DynamicLegend');
% end
% 
% % Plot male - female bar graph for each scenario. -------------------------
% 
% % X_m/f = (variable, iteration, scenario)
% deltaMAP_m = reshape(X_m(42,end,1:end-2) - X_m(42,iteration,1:end-2), [1,num_scen-2]);
% deltaMAP_f = reshape(X_f(42,end,1:end-2) - X_f(42,iteration,1:end-2), [1,num_scen-2]);
% MAP_comp = deltaMAP_m(1) - [deltaMAP_m(1), deltaMAP_f];
% scen_comp = categorical({'M - M'       , 'M - F'       , ...
%                          'M - F M RSNA', 'M - F M AT2R', ...
%                          'M - F M RAS' , 'M - F M Reab'});
% scen_comp = reordercats(scen_comp,{'M - M'       , 'M - F'       , ...
%                                    'M - F M RSNA', 'M - F M AT2R', ...
%                                    'M - F M RAS' , 'M - F M Reab'});
% 
% f7 = figure('DefaultAxesFontSize',10);
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.15, 3.5]);
% s1(1) = subplot(1,2,1); 
% s1(2) = subplot(1,2,2); 
% 
% plot(s1(1), X_m(42,:,fixed_ss),xscale,'-', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3);
% xlim(s1(1), [90, 120]); 
% % set(s1(1),'XTick', [80,100,120]);
% ylim(s1(1), [lower, upper])
% xlabel(s1(1), 'MAP (mmHg)', 'FontSize',14*1.1); ylabel(s1(1), {'Fold change in'; 'sodium excretion'}, 'FontSize',14);
% hold(s1(1), 'on')
% plot(s1(1), X_f(42,:,fixed_ss),xscale,'-', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3)
% legend(s1(1), {'Male','Female'}, 'Location','Southeast', 'FontSize',14)
% hold(s1(1), 'off')
% title(s1(1), 'A', 'FontSize',14)
% 
% bar(s1(2), scen_comp,MAP_comp,'k');
% % set(gca,'xticklabel',scen_comp_text);
% % xtickget = get(gca,'xticklabel');  
% % set(gca,'xticklabel',xtickget,'fontsize',6)
% % xtickangle(s1(2),90)
% % xlim(s1(2), [1-1,6+1])
% ylim(s1(2), [-2,5])
% xlabel(s1(2), 'Scenario', 'FontSize',14); ylabel(s1(2), '\DeltaMAP (mmHg)', 'FontSize',14);
% % hAxes.XAxis.FontSize = 6;
% title(s1(2), 'B', 'FontSize',14)

% % % BAR GRAPH IS PRODUCING NO CHANG BETWEEN M AND F SCENARIOS BECAUSE
% PARS ARE NOT BEING RETREIVED FROM 'get_pars.m'. % % %

f = [f1;f2;f3;f4];

end % solve_ss_Phisodin

% -------------------------------------------------------------------------
% Renal perfusion pressure
% -------------------------------------------------------------------------

function f = run_sim_RPP(PARS,SSDATA,sample_num)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Begin user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Renal perfusion pressure perturbation
% Enter postive for increase or negative for decrease.
RPP_per = [-20; 0; 20];
num_per = length(RPP_per);
% Index of RPP to plot for all variables
exact_per = 3;

% Scenarios
scenario = {'Normal'};
% scenario = {'Pri_Hyp'};
num_scen = length(scenario);
% Index of scenario to plot for all variables
exact_scen = 1;

% Number of points for plotting resolution
num_points = 121;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize variables.
% X = (variables, points, sex, perturbation, scenario)
X = zeros(num_vars,num_points,2,num_per,num_scen);
% Retrieve male/female. 
% X_m/f = (variables, points, perturbation, scenario)
X_m = zeros(num_vars,num_points,num_per,num_scen);
X_f = zeros(num_vars,num_points,num_per,num_scen);

% Need to store male and female RPP for plotting later.
% RPP = (sex, scenario)
RPP = zeros(2,num_scen);

for per_ind     = 1:num_per  % perturbation
for sce_ind     = 1:num_scen % scenario
for sex_ind_RPP = 1:2        % sex

% Retreive parameters and steady state data.
pars = PARS{sex_ind_RPP}; SSdata = SSDATA{sex_ind_RPP}; 
% % % Note: SSdata is cleared above. Therefore it is not a scoping variable. % % % 

pars = pars(:,sample_num); SSdata = SSdata(:,sample_num);

varargin_input = {'RPP',{RPP_per(per_ind), SSdata}, scenario{sce_ind},true, ...
                  'Denerve',{true, SSdata}, 'Fixed Water Intake',{true, SSdata}};

%% Solve DAE

% Renal Perfusion Pressure.
RPP(sex_ind_RPP,sce_ind) = SSdata(42);

% Initial condition for the variables and their derivatives. 
% System is initially at steady state, so the derivative is 0.
x0 = SSdata; x_p0 = zeros(num_vars,1);

% Time at which to keep steady state, change a parameter, etc.
tchange = 10;

% Initial time (min); Final time (min); Points per minute;
t0 = 0; tend = tchange + 50; ppm = (num_points-1)/(tend-t0);

% Time vector
tspan = linspace(t0,tend,num_points);

% ode options
options = odeset('RelTol',1e-1, 'AbsTol',1e-2, 'MaxStep',1e-2);

% Solve dae
[t,x] = ode15i(@(t,x,x_p) ...
               bp_reg_mod(t,x,x_p,pars,tchange,varargin_input{:}), ...
               tspan, x0, x_p0, options);
t = t'; x = x';

% Store solution.
% X = (variables, points, sex, perturbation, scenario)
X(:,:,sex_ind_RPP,per_ind,sce_ind) = x;

end % sex
end % scenario
end % perturbation

%% Retrieve data and visualize

% X_m/f = (variables, points, perturbation, scenario)
X_m(:,:,:,:) = X(:,:,1,:,:);
X_f(:,:,:,:) = X(:,:,2,:,:);

% x-axis limits
xlower = t0; xupper = tend; 

% y-axis limits
ylower = zeros(num_vars); yupper = ylower; 
for i = 1:num_vars
    ylower(i) = 0.9*min(min(X_m(i,:,exact_per,exact_scen)), min(X_f(i,:,exact_per,exact_scen)));
    yupper(i) = 1.1*max(max(X_m(i,:,exact_per,exact_scen)), max(X_f(i,:,exact_per,exact_scen)));
    if abs(yupper(i)) < eps*100
        ylower(i) = -10^(-5); yupper(i) = 10^(-5);
    end
end

% % Plot all variables vs time. ---------------------------------------------
% 
% f1  = gobjects(7,1);
% s1 = gobjects(7,15);
% % Loop through each set of subplots.
% for i = 1:7
%     f1(i) = figure('pos',[750 500 650 450]);
%     % This is to avoid the empty plots in the last subplot set.
%     if i == 7
%         last_plot = mod(num_vars, 15);
%     else
%         last_plot = 15;
%     end
%     % Loop through each subplot within a set of subplots.
%     for j = 1:last_plot
%         s1(i,j) = subplot(3,5,j);
%         s1(i,j).Position = s1(i,j).Position + [0 0 0.01 0];
%         
%         plot(s1(i,j), t,X_m((i-1)*15+j,:,exact_per,exact_scen),'b' , ...
%                       t,X_f((i-1)*15+j,:,exact_per,exact_scen),'r');
%         
% %         xlim([xlower, xupper])
%         ylim([ylower((i-1)*15 + j), yupper((i-1)*15 + j)])
% 
% %         legend('Male', 'Female')
%         title(names((i-1)*15 + j), 'Interpreter','latex', 'FontSize',15)
%     end
% end
% 
% % Plot renal perfusion pressure input vs time. ----------------------------
% 
% tplot   = [t0:1:tend];
% RPPplot = zeros(1,length(tplot));
% RPPplot(1        :tchange) = RPP(1,exact_scen);
% RPPplot(tchange+1:tend+1 ) = RPP(1,exact_scen) + RPP_per(exact_per);
% g = figure('pos',[100 100 675 450]);
% plot(tplot,RPPplot, 'LineWidth',3)
% xlabel('$t$ (min)', 'Interpreter','latex')
% ylabel('$RPP$'    , 'Interpreter','latex')

% Plot data as in Hilliard 2011. ------------------------------------------

% Time average quantity from 10-30 minutes after perturbation in RPP.
% RPP at 80, 100, 120.
% Phi_rb = var(6), Phi_gfilt = var(7), Phi_u = var(92), Phi_usod = var(26)
% rel = value divided by baseline value; act = actual value

% X_m/f = (variables, points, perturbation, scenario)

time_int    = (tchange+10)*ppm+1:(tchange+30)*ppm+1;
time_points = length(time_int);
RBF_rel_m  = zeros(num_per,num_scen); RBF_rel_f  = zeros(num_per,num_scen);  
GFR_rel_m  = zeros(num_per,num_scen); GFR_rel_f  = zeros(num_per,num_scen); 
UF_rel_m   = zeros(num_per,num_scen); UF_rel_f   = zeros(num_per,num_scen); 
USOD_rel_m = zeros(num_per,num_scen); USOD_rel_f = zeros(num_per,num_scen); 

for sce_ind = 1:num_scen
    for per_ind = 1:num_per
        RBF_rel_m (per_ind,sce_ind) = (sum(X_m(6 , time_int, per_ind, sce_ind)) / time_points) ...
                                    / (sum(X_m(6 , time_int, 2 , sce_ind)) / time_points);
        GFR_rel_m (per_ind,sce_ind) = (sum(X_m(7 , time_int, per_ind, sce_ind)) / time_points) ...
                                    / (sum(X_m(7 , time_int, 2 , sce_ind)) / time_points);
        UF_rel_m  (per_ind,sce_ind) = (sum(X_m(92, time_int, per_ind, sce_ind)) / time_points) ...
                                    / (sum(X_m(92, time_int, 2 , sce_ind)) / time_points);
        USOD_rel_m(per_ind,sce_ind) = (sum(X_m(27, time_int, per_ind, sce_ind)) / time_points) ...
                                    / (sum(X_m(27, time_int, 2 , sce_ind)) / time_points);
        
        RBF_rel_f (per_ind,sce_ind) = (sum(X_f(6 , time_int, per_ind, sce_ind)) / time_points) ...
                                    / (sum(X_f(6 , time_int, 2 , sce_ind)) / time_points);
        GFR_rel_f (per_ind,sce_ind) = (sum(X_f(7 , time_int, per_ind, sce_ind)) / time_points) ...
                                    / (sum(X_f(7 , time_int, 2 , sce_ind)) / time_points);
        UF_rel_f  (per_ind,sce_ind) = (sum(X_f(92, time_int, per_ind, sce_ind)) / time_points) ...
                                    / (sum(X_f(92, time_int, 2 , sce_ind)) / time_points);
        USOD_rel_f(per_ind,sce_ind) = (sum(X_f(27, time_int, per_ind, sce_ind)) / time_points) ...
                                    / (sum(X_f(27, time_int, 2 , sce_ind)) / time_points);
    end
end

% RPP
RPP_m = RPP(1,exact_scen) + RPP_per; RPP_f = RPP(2,exact_scen) + RPP_per; 

% Data --------------------------------------------------------------------

% Yes AT2R
RBFdata_yes_at2r_rel_m  = [0.8894; 1.0000; 1.0609]; RBFdata_yes_at2r_rel_f  = [0.8562; 1.0000; 1.1045]; 
GFRdata_yes_at2r_rel_m  = [0.7816; 1.0000; 1.0083]; GFRdata_yes_at2r_rel_f  = [0.7395; 1.0000; 1.1357]; 
UFdata_yes_at2r_rel_m   = [0.6163; 1.0000; 1.4535]; UFdata_yes_at2r_rel_f   = [0.8097; 1.0000; 2.2327]; 
USODdata_yes_at2r_rel_m = [0.4000; 1.0000; 1.8744]; USODdata_yes_at2r_rel_f = [0.5979; 1.0000; 3.0815]; 
% Block AT2R
RBFdata_blk_at2r_rel_m  = [0.8652; 1.0000; 1.2381]; RBFdata_blk_at2r_rel_f  = [0.5404; 1.0000; 1.0961]; 
GFRdata_blk_at2r_rel_m  = [0.7476; 1.0000; 1.2907]; GFRdata_blk_at2r_rel_f  = [0.2642; 1.0000; 1.3678]; 
UFdata_blk_at2r_rel_m   = [0.7106; 1.0000; 1.7368]; UFdata_blk_at2r_rel_f   = [0.4597; 1.0000; 2.7119]; 
USODdata_blk_at2r_rel_m = [0.0000; 1.0000; 3.5712]; USODdata_blk_at2r_rel_f = [0.5161; 1.0000; 4.2295]; 

% Male combined array for ease of plotting
RBFdata_rel_m  = [zeros(3,1), RBFdata_yes_at2r_rel_m , RBFdata_blk_at2r_rel_m ];
GFRdata_rel_m  = [zeros(3,1), GFRdata_yes_at2r_rel_m , GFRdata_blk_at2r_rel_m ];
UFdata_rel_m   = [zeros(3,1), UFdata_yes_at2r_rel_m  , UFdata_blk_at2r_rel_m  ];
USODdata_rel_m = [zeros(3,1), USODdata_yes_at2r_rel_m, USODdata_blk_at2r_rel_m];
% Female combined array for ease of plotting
RBFdata_rel_f  = [zeros(3,1), RBFdata_yes_at2r_rel_f , RBFdata_blk_at2r_rel_f ];
GFRdata_rel_f  = [zeros(3,1), GFRdata_yes_at2r_rel_f , GFRdata_blk_at2r_rel_f ];
UFdata_rel_f   = [zeros(3,1), UFdata_yes_at2r_rel_f  , UFdata_blk_at2r_rel_f  ];
USODdata_rel_f = [zeros(3,1), USODdata_yes_at2r_rel_f, USODdata_blk_at2r_rel_f];

% Subplot -----------------------------------------------------------------

f3 = figure('DefaultAxesFontSize',14);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.15, 5]);
s_rel1(1) = subplot(2,2,1); 
s_rel1(2) = subplot(2,2,2); 
s_rel1(3) = subplot(2,2,3);
s_rel1(4) = subplot(2,2,4); 

plot(s_rel1(1), RPP_m,RBF_rel_m     (:,exact_scen) ,'x-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
% xlim(s_rel1(1), [75,125]); set(s_rel1(1),'XTick', [80,100,120]);
% ylim(s_rel1(1), [0.6,1.2])
xlabel(s_rel1(1), 'RPP (mmHg)'); ylabel(s_rel1(1), 'RBF (relative)');
hold(s_rel1(1), 'on')
plot(s_rel1(1), RPP_m,RBFdata_rel_m (:,2         ) ,'o--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(s_rel1(1), RPP_f,RBF_rel_f     (:,exact_scen) ,'x-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(s_rel1(1), RPP_f,RBFdata_rel_f (:,2         ) ,'o--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8); 
hold(s_rel1(1), 'off')
[~, hobj, ~, ~] = legend(s_rel1(1), {'Male sim','Male data','Female sim','Female data'}, 'FontSize',7,'Location','Southeast');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
title(s_rel1(1), 'A')

plot(s_rel1(2), RPP_m,GFR_rel_m     (:,exact_scen) ,'x-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
% xlim(s_rel1(2), [75,125]); set(s_rel1(2),'XTick', [80,100,120]);
% ylim(s_rel1(2), [0.6,1.2])
xlabel(s_rel1(2), 'RPP (mmHg)'); ylabel(s_rel1(2), 'GFR (relative)');
hold(s_rel1(2), 'on')
plot(s_rel1(2), RPP_m,GFRdata_rel_m (:,2         ) ,'o--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(s_rel1(2), RPP_f,GFR_rel_f     (:,exact_scen) ,'x-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(s_rel1(2), RPP_f,GFRdata_rel_f (:,2         ) ,'o--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8); 
hold(s_rel1(2), 'off')
title(s_rel1(2), 'B')

plot(s_rel1(3), RPP_m,UF_rel_m      (:,exact_scen) ,'x-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
% xlim(s_rel1(3), [75,125]); set(s_rel1(3),'XTick', [80,100,120]);
% ylim(s_rel1(3), [0.0,3.5])
xlabel(s_rel1(3), 'RPP (mmHg)'); ylabel(s_rel1(3), 'UF (relative)');
hold(s_rel1(3), 'on')
plot(s_rel1(3), RPP_m,UFdata_rel_m  (:,2         ) ,'o--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(s_rel1(3), RPP_f,UF_rel_f      (:,exact_scen) ,'x-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(s_rel1(3), RPP_f,UFdata_rel_f  (:,2         ) ,'o--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8); 
hold(s_rel1(3), 'off')
title(s_rel1(3), 'C')

plot(s_rel1(4), RPP_m,USOD_rel_m    (:,exact_scen) ,'x-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
% xlim(s_rel1(4), [75,125]); set(s_rel1(4),'XTick', [80,100,120]);
% ylim(s_rel1(4), [0.0,3.5])
xlabel(s_rel1(4), 'RPP (mmHg)'); ylabel(s_rel1(4), 'USOD (relative)');
hold(s_rel1(4), 'on')
plot(s_rel1(4), RPP_m,USODdata_rel_m(:,2         ) ,'o--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(s_rel1(4), RPP_f,USOD_rel_f    (:,exact_scen) ,'x-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(s_rel1(4), RPP_f,USODdata_rel_f(:,2         ) ,'o--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8); 
hold(s_rel1(4), 'off')
title(s_rel1(4), 'D')

f = [f3];

end

end % sim_hyp_vp





























