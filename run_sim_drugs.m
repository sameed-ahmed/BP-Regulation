% This simulates the blood pressure regulation model bp_reg.m for drug administration.
% 
% Steady state data is calculated by solve_ss_scenario.m.

function run_sim_drugs

close all

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Begin user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Physiological scenarios
% Normal  - Normal conditions
% m_RAS   - male RAS pars
% m_Reab  - male fractional sodium and water reabsorption
% Pri_Hyp - essential/primary hypertension
scenario1 = {'Normal', 'm_RSNA', 'm_AT2R', 'm_RAS', 'm_Reab', ...
             'm_RAS_m_Reab', 'm_RSNA_m_Reab'};
% scenario1 = {'Normal', 'm_RSNA', 'm_AT2R', 'm_RAS', 'm_Reab', ...
%              'm_RAS_m_Reab', 'm_RSNA_m_Reab', ...
%              'Pri_Hyp'};
fixed_ss1 = 1;
% Drug scenarios
% Normal - Normal conditions
% AngII  - Ang II infusion fmol/(ml min)
% ACEi   - Angiotensin converting enzyme inhibitor %
% ARB1   - Angiotensin receptor 1 blocker %
% ARB2   - Angiotensin receptor 2 blocker %
% DRI    - Direct renin inhibitor %
% MRB    - Aldosterone blocker (MR?) %
% RSS    - Renin secretion stimulator (thiazide?) % % NOT COMPLETE
scenario2 = {'Normal', 'AngII', 'ACEi', 'ARB1', 'ARB2', 'DRI', 'MRB', 'RSS'};
fixed_ss2 = [4];
num_scen = length(scenario1);

% Species
spe_ind = 2;

% Number of days to run simulation after change; Day at which to induce change;
days = 14; day_change = 1;
% Number of points for plotting resolution
% N = ((days+1)*1440) / 2;
N = (days+1)*100 + 1;

% Bootstrap replicate sample number
sample_num = random('Discrete Uniform',1000)
% sample_num = 42 % male and female MAP similar
% sample_num = 208
% sample_num = 655

% Drug dose
drug_dose = 0.95

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

species = {'human', 'rat'   };
sex     = {'male' , 'female'};

% Number of variables
num_vars = 93;

% Initialize variables.
% X = (variables, points, sex, scenario)
X_dy = zeros(num_vars,N,2,num_scen);
% X = (variables, sex, scenario)
X_ss = zeros(num_vars,2,num_scen);

for sce_ind = fixed_ss1:fixed_ss1 % scenario
for sex_ind = 1:2        % sex

varargin_input = {scenario1{sce_ind},true};

%% Load bootstrap replicate parameters & variables created by create_par_bs_rep.m.

% Parameters
load_data_name_pars = sprintf('%s_%s_pars_scenario_Pri_Hyp_bs_rep1000.mat', ...
                              species{spe_ind},sex{sex_ind});
load(load_data_name_pars, 'pars_rep');
num_pars   = size(pars_rep,1);
num_sample = size(pars_rep,2);
pars_rep = pars_rep(:,sample_num);

% Variables
load_data_name_vars = sprintf('%s_%s_ss_data_scenario_Pri_Hyp_bs_rep1000.mat', ...
                              species{spe_ind},sex{sex_ind});
load(load_data_name_vars, 'SSdata_rep');
num_vars   = size(SSdata_rep,1);
SSdata_rep = SSdata_rep(:,sample_num);

%% Drugs

for i = 1:length(fixed_ss2)
    if     strcmp(scenario2{fixed_ss2(i)}, 'AngII')
        if     strcmp(sex{sex_ind}, 'male')
            varargin_input = [varargin_input, 'AngII',910]; % Sullivan 2010
        elseif strcmp(sex{sex_ind}, 'female')
            varargin_input = [varargin_input, 'AngII',505]; % Sullivan 2010
        end
    elseif strcmp(scenario2{fixed_ss2(i)}, 'ACEi' )
            varargin_input = [varargin_input, 'ACEi' ,drug_dose]; % 
    elseif strcmp(scenario2{fixed_ss2(i)}, 'ARB1' )
            varargin_input = [varargin_input, 'ARB1' ,drug_dose]; % 
    elseif strcmp(scenario2{fixed_ss2(i)}, 'ARB2' )
            varargin_input = [varargin_input, 'ARB2' ,drug_dose]; % 
    elseif strcmp(scenario2{fixed_ss2(i)}, 'DRI'  )
            varargin_input = [varargin_input, 'DRI'  ,drug_dose]; % 
    elseif strcmp(scenario2{fixed_ss2(i)}, 'MRB'  )
            varargin_input = [varargin_input, 'MRB'  ,drug_dose]; % 
    elseif strcmp(scenario2{fixed_ss2(i)}, 'RSS'  )
            varargin_input = [varargin_input, 'RSS'  ,drug_dose]; % 
    end
end

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

%% Solve DAE dynamic

% Initial condition for the variables and their derivatives. 
% System is initially at steady state, so the derivative is 0.
x0_dy = SSdata_rep; x_p0_dy = zeros(num_vars,1);

% Time at which to keep steady state, change a parameter, etc.
tchange_dy = day_change*1440;

% Initial time (min); Final time (min);
t0 = 0*1440; tend = tchange_dy + days*1440;

% Time vector
tspan = linspace(t0,tend,N);

% ODE options
options_dy = odeset('MaxStep',1000); % default is 0.1*abs(t0-tf)
% Solve dae
[t_dy,x] = ode15i(@(t,x,x_p) ...
               bp_reg_mod(t,x,x_p,pars_rep,tchange_dy,varargin_input{:}), ...
               tspan, x0_dy, x_p0_dy, options_dy);

% X = (variables, points, sex, scenario)
X_dy(:,:,sex_ind,sce_ind) = x';

% %% Solve system steady state
% 
% % Initial guess for the variables.
% % Find the steady state solution, so the derivative is 0.
% % Arbitrary value for time to input, greater than tchange + deltat.
% x0_ss = SSdata_rep; x_p0_ss = zeros(num_vars,1); t_ss = 30;
% 
% % Time at which to change place holder.
% tchange_ss = 0;
% 
% % Solver options
% options_ss = optimset();
% % Solve system
% [SSdata, residual, ...
%  exitflag, output] = fsolve(@(x) ...
%                             bp_reg_mod(t_ss,x,x_p0_ss,pars_rep,tchange_ss,varargin_input{:}), ...
%                             x0_ss, options_ss);
% 
% % Check for solver convergence.
% if exitflag == 0
%     disp('Solver did not converge.')
%     disp(output)
% end
% 
% % Check for imaginary solution.
% if not (isreal(SSdata))
%     disp('Imaginary number returned.')
% end
% 
% % Set any values that are within machine precision of 0 equal to 0.
% for i = 1:length(SSdata)
%     if abs(SSdata(i)) < eps*100
%         SSdata(i) = 0;
%     end
% end
% 
% % X = (variables, sex, scenario)
% X_ss(:,sex_ind,sce_ind) = SSdata;

end % sex
end % scenario

%% Plot

% Retrieve male and female.
% X_m/f = (variables, points, scenario)
t_dy = t_dy';
X_dy_m = reshape(X_dy(:,:,1,:), [num_vars,N,num_scen]); 
X_dy_f = reshape(X_dy(:,:,2,:), [num_vars,N,num_scen]); 
% % X_m/f = (variables, scenario)
% X_ss_m = reshape(X_ss(:,  1,:), [num_vars,  num_scen]); 
% X_ss_f = reshape(X_ss(:,  2,:), [num_vars,  num_scen]); 

% x-axis limits
xlower = t0; xupper = tend; 

% Convert minutes to days for longer simulations.
t_dy = t_dy/1440; tchange_dy = tchange_dy/1440; 
xlower = xlower/1440; xupper = xupper/1440; 

% y-axis limits
ylower = zeros(num_vars,1); yupper = ylower; 
for i = 1:length(ylower)
    ylower(i) = 0.95*min( min(X_dy_m(i,:,fixed_ss1)), min(X_dy_f(i,:,fixed_ss1)) );
    yupper(i) = 1.05*max( max(X_dy_m(i,:,fixed_ss1)), max(X_dy_f(i,:,fixed_ss1)) );
    if abs(yupper(i)) < eps*100
        ylower(i) = -10^(-5); yupper(i) = 10^(-5);
    end
end

%% Plot all vars vs time. -------------------------------------------------

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
        
        plot(s(i,j), t_dy,X_dy_m((i-1)*15 + j,:,fixed_ss1),'b', ...
                     t_dy,X_dy_f((i-1)*15 + j,:,fixed_ss1),'r');
        
        xlim([xlower, xupper])
        ylim([ylower((i-1)*15 + j), yupper((i-1)*15 + j)])
        
%         Days
        ax.XTick = (tchange_dy+0*(1) : 2 : tchange_dy+days*(1));
        ax.XTickLabel = {'0','2','4','6','8','10','12','14'};

%         legend('Male', 'Female')
        title(names((i-1)*15 + j), 'Interpreter','latex', 'FontSize',15)
    end
end

%% Plot interesting variables. --------------------------------------------

% Interesting variables to plot.
var_ind = [33;41;42;9;73;74;6;7;27;92;93;29]; sub_var_num = length(var_ind);

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

    plot(s2(j), t_dy,X_dy_m(var_ind(j),:,fixed_ss1), 'Color',[0.203, 0.592, 0.835], 'LineWidth',2.5);
    hold(s2(j), 'on')
    plot(s2(j), t_dy,X_dy_f(var_ind(j),:,fixed_ss1), 'Color',[0.835, 0.203, 0.576], 'LineWidth',2.5);
    hold(s2(j), 'off')

    xlim([xlower, xupper])
    ylim([ylower(var_ind(j)), yupper(var_ind(j))])
    
    set(s2(j), 'XTick', [tchange_dy+0*(1) : 2 : tchange_dy+days*(1)]);
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

%% Plot Mean Arterial Pressure vs Time. -----------------------------------

% Substract MAP by baseline for each sex and all scenarios.
% X_m/f = (variable, points, scenario)
MAP_m = zeros(N,num_scen); MAP_f = zeros(N,num_scen);
for i = 1:num_scen
    MAP_m(:,i) = (X_dy_m(42,:,i) - X_dy_m(42,1,i)) ./ X_dy_m(42,1,i) * 100;
    MAP_f(:,i) = (X_dy_f(42,:,i) - X_dy_f(42,1,i)) ./ X_dy_f(42,1,i) * 100;
end
% MAP_m = reshape(X_m(42,:,i) - X_m(42,1,i), [N,num_scen]);
% MAP_f = reshape(X_f(42,:,i) - X_f(42,1,i), [N,num_scen]);

g1 = figure('DefaultAxesFontSize',14);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 2.5]);
plot(t_dy,MAP_m(:,fixed_ss1),'-', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3);
xlim([xlower, xupper]); %ylim([-20, 5]);
ax = gca;
ax.XTick = (tchange_dy+0*(1) : 2 : tchange_dy+days*(1));
ax.XTickLabel = {'0','2','4','6','8','10','12','14'};
xlabel('Time (days)'); ylabel('\DeltaMAP (mmHg)');
hold on
plot(t_dy,MAP_f(:,fixed_ss1),'-', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3)
[~, hobj, ~, ~] = legend({'Male sim','Female sim'}, 'FontSize',7,'Location','Northwest');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);

%% Plot all other quantities of interest. ---------------------------------

% GFR; BV; RSNA; REA/RR for each sex and all scenarios.
% X_m/f = (variable, points, scenario)
GFR_m  = reshape(X_dy_m( 7,:,:), [N,num_scen]);
GFR_f  = reshape(X_dy_f( 7,:,:), [N,num_scen]);
BV_m   = reshape(X_dy_m(30,:,:), [N,num_scen]);
BV_f   = reshape(X_dy_f(30,:,:), [N,num_scen]);
R_m    = reshape(X_dy_m(74,:,:) ./ X_dy_m( 4,:,:), [N,num_scen]);
R_f    = reshape(X_dy_f(74,:,:) ./ X_dy_f( 4,:,:), [N,num_scen]);
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
FRNA_m = reshape((X_dy_m(11,:,:) - X_dy_m(27,:,:)) ./ X_dy_m(11,:,:), [N,num_scen]) * 100;
FRNA_f = reshape((X_dy_f(11,:,:) - X_dy_f(27,:,:)) ./ X_dy_f(11,:,:), [N,num_scen]) * 100;
FRW_m  = reshape((X_dy_m( 7,:,:) - X_dy_m(92,:,:)) ./ X_dy_m( 7,:,:), [N,num_scen]) * 100;
FRW_f  = reshape((X_dy_f( 7,:,:) - X_dy_f(92,:,:)) ./ X_dy_f( 7,:,:), [N,num_scen]) * 100;
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

g2 = figure('DefaultAxesFontSize',14);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.15, 5]);
s_2(1) = subplot(2,2,1); 
s_2(2) = subplot(2,2,2); 
s_2(3) = subplot(2,2,3);
s_2(4) = subplot(2,2,4); 

plot(s_2(1), t_dy,R_m   (:,fixed_ss1), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim(s_2(1), [xlower, xupper]);
set(s_2(1), 'XTick', [tchange_dy+0*(1) : 2 : tchange_dy+days*(1)]);
set(s_2(1), 'XTickLabel', {'0','2','4','6','8','10','12','14'});
% ylim(s_main(1), [0.75,1.05])
xlabel(s_2(1), 'Time (days)'); ylabel(s_2(1), 'R_{EA}/R_R (relative)');
hold(s_2(1), 'on')
plot(s_2(1), t_dy,R_f   (:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold(s_2(1), 'off')
[~, hobj, ~, ~] = legend(s_2(1), {'Male','Female'}, 'FontSize',7,'Location','Northeast');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
title(s_2(1), 'A')

plot(s_2(2), t_dy,FRNA_m(:,fixed_ss1) ,'-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim(s_2(2), [xlower, xupper]);
set(s_2(2), 'XTick', [tchange_dy+0*(1) : 2 : tchange_dy+days*(1)]);
set(s_2(2), 'XTickLabel', {'0','2','4','6','8','10','12','14'});
% ylim(s_main(2), [97,100])
xlabel(s_2(2), 'Time (days)'); ylabel(s_2(2), 'FR (relative)');
hold(s_2(2), 'on')
plot(s_2(2), t_dy,FRW_m (:,fixed_ss1), '--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(s_2(2), t_dy,FRNA_f(:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(s_2(2), t_dy,FRW_f (:,fixed_ss1), '--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
fakeplot = zeros(2, 1);
fakeplot(1) = plot(s_2(2), NaN,NaN, 'k-' );
fakeplot(2) = plot(s_2(2), NaN,NaN, 'k--');
[~, hobj, ~, ~] = legend(fakeplot, {'FR_{Na^+}','FR_{U}'}, 'FontSize',7,'Location','Northeast');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
hold(s_2(2), 'off')
title(s_2(2), 'B')

plot(s_2(3), t_dy,GFR_m (:,fixed_ss1), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim(s_2(3), [xlower, xupper]);
set(s_2(3), 'XTick', [tchange_dy+0*(1) : 2 : tchange_dy+days*(1)]);
set(s_2(3), 'XTickLabel', {'0','2','4','6','8','10','12','14'});
% ylim(s_main(3), [0.75,1.35])
xlabel(s_2(3), 'Time (days)'); ylabel(s_2(3), 'GFR (relative)');
hold(s_2(3), 'on')
plot(s_2(3), t_dy,GFR_f (:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold(s_2(3), 'off')
title(s_2(3), 'C')

plot(s_2(4), t_dy,BV_m  (:,fixed_ss1), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim(s_2(4), [xlower, xupper]);
set(s_2(4), 'XTick', [tchange_dy+0*(1) : 2 : tchange_dy+days*(1)]);
set(s_2(4), 'XTickLabel', {'0','2','4','6','8','10','12','14'});
% ylim(s_main(4), [1,1.25])
xlabel(s_2(4), 'Time (days)'); ylabel(s_2(4), 'BV (relative)');
hold(s_2(4), 'on')
plot(s_2(4), t_dy,BV_f  (:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold(s_2(4), 'off')
title(s_2(4), 'D')

%% Plot all other quantities of interest. ---------------------------------

% GFR; BV; RSNA; REA/RR for each sex and all scenarios.
% X_m/f = (variable, points, scenario)
MAP_m   = reshape(X_dy_m(42,:,:), [N,num_scen]);
MAP_f   = reshape(X_dy_f(42,:,:), [N,num_scen]);
PRA_m   = reshape(X_dy_m(66,:,:), [N,num_scen]);
PRA_f   = reshape(X_dy_f(66,:,:), [N,num_scen]);
ANGI_m  = reshape(X_dy_m(67,:,:), [N,num_scen]);
ANGI_f  = reshape(X_dy_f(67,:,:), [N,num_scen]);
ANGII_m = reshape(X_dy_m(68,:,:), [N,num_scen]);
ANGII_f = reshape(X_dy_f(68,:,:), [N,num_scen]);
% Plot as relative change in order to compare male and female.
MAP_m_bl   = MAP_m  (1,:);
MAP_f_bl   = MAP_f  (1,:);
PRA_m_bl   = PRA_m  (1,:);
PRA_f_bl   = PRA_f  (1,:);
ANGI_m_bl  = ANGI_m (1,:);
ANGI_f_bl  = ANGI_f (1,:);
ANGII_m_bl = ANGII_m(1,:);
ANGII_f_bl = ANGII_f(1,:);
for i = 1:N
    MAP_m  (i,:) = MAP_m  (i,:) ./ MAP_m_bl  ;
    MAP_f  (i,:) = MAP_f  (i,:) ./ MAP_f_bl  ;
    PRA_m  (i,:) = PRA_m  (i,:) ./ PRA_m_bl  ;
    PRA_f  (i,:) = PRA_f  (i,:) ./ PRA_f_bl  ;
    ANGI_m (i,:) = ANGI_m (i,:) ./ ANGI_m_bl ;
    ANGI_f (i,:) = ANGI_f (i,:) ./ ANGI_f_bl ;
    ANGII_m(i,:) = ANGII_m(i,:) ./ ANGII_m_bl;
    ANGII_f(i,:) = ANGII_f(i,:) ./ ANGII_f_bl;
end

g3 = figure('DefaultAxesFontSize',14);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.15, 5]);
s_3(1) = subplot(2,2,1); 
s_3(2) = subplot(2,2,2); 
s_3(3) = subplot(2,2,3);
s_3(4) = subplot(2,2,4); 

plot(s_3(1), t_dy,MAP_m   (:,fixed_ss1), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim(s_3(1), [xlower, xupper]);
set(s_3(1), 'XTick', [tchange_dy+0*(1) : 2 : tchange_dy+days*(1)]);
set(s_3(1), 'XTickLabel', {'0','2','4','6','8','10','12','14'});
% ylim(s_main(1), [0.75,1.05])
xlabel(s_3(1), 'Time (days)'); ylabel(s_3(1), 'MAP (relative)');
hold(s_3(1), 'on')
plot(s_3(1), t_dy,MAP_f   (:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold(s_3(1), 'off')
[~, hobj, ~, ~] = legend(s_3(1), {'Male','Female'}, 'FontSize',7,'Location','Northeast');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
title(s_3(1), 'A')

plot(s_3(2), t_dy,PRA_m(:,fixed_ss1) ,'-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim(s_3(2), [xlower, xupper]);
set(s_3(2), 'XTick', [tchange_dy+0*(1) : 2 : tchange_dy+days*(1)]);
set(s_3(2), 'XTickLabel', {'0','2','4','6','8','10','12','14'});
% ylim(s_main(2), [97,100])
xlabel(s_3(2), 'Time (days)'); ylabel(s_3(2), 'PRA (relative)');
hold(s_3(2), 'on')
plot(s_3(2), t_dy,PRA_m (:,fixed_ss1), '--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
hold(s_3(2), 'off')
title(s_3(2), 'B')

plot(s_3(3), t_dy,ANGI_m (:,fixed_ss1), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim(s_3(3), [xlower, xupper]);
set(s_3(3), 'XTick', [tchange_dy+0*(1) : 2 : tchange_dy+days*(1)]);
set(s_3(3), 'XTickLabel', {'0','2','4','6','8','10','12','14'});
% ylim(s_main(3), [0.75,1.35])
xlabel(s_3(3), 'Time (days)'); ylabel(s_3(3), 'Ang I (relative)');
hold(s_3(3), 'on')
plot(s_3(3), t_dy,ANGI_f (:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold(s_3(3), 'off')
title(s_3(3), 'C')

plot(s_3(4), t_dy,ANGII_m  (:,fixed_ss1), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim(s_3(4), [xlower, xupper]);
set(s_3(4), 'XTick', [tchange_dy+0*(1) : 2 : tchange_dy+days*(1)]);
set(s_3(4), 'XTickLabel', {'0','2','4','6','8','10','12','14'});
% ylim(s_main(4), [1,1.25])
xlabel(s_3(4), 'Time (days)'); ylabel(s_3(4), 'Ang II (relative)');
hold(s_3(4), 'on')
plot(s_3(4), t_dy,ANGII_f  (:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold(s_3(4), 'off')
title(s_3(4), 'D')

% % % Save figures. -----------------------------------------------------------
% % 
% save_data_name = sprintf('all_vars_%s.fig', scenario2{fixed_ss2});
% save_data_name = strcat('Figures/', save_data_name);
% savefig([f;f2;g;h';k], save_data_name)

end






























