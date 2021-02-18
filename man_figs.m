% This script to output the specific figures in Ahmed 2021.

function man_figs

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

% Physiological scenarios
% Normal  - Normal conditions
scenario1 = {'Normal'};
fixed_ss1 = 1;
num_scen = length(scenario1);
% Drug scenarios
% Normal - Normal conditions
% ACEi   - Angiotensin converting enzyme inhibitor % 96
% ARB1   - Angiotensin receptor 1 blocker % 94
% CCB    - Calcium channel blocker % 85
% TZD    - Thiazide diuretic % 50 100
scenario2 = {'Normal', 'ACEi', 'ARB1', 'CCB', 'TZD'};
fixed_ss2 = [5];

% Species
spe_ind = 2;

% Number of days to run simulation after change; Day at which to induce change;
days = 7; day_change = 1;
% Number of points for plotting resolution
% N = ((days+1)*1440) / 2;
N = (days+1)*100 + 1;

% Bootstrap replicate sample number
% sample_num = random('Discrete Uniform',1000)
% sample_num = 42 % male and female MAP similar
sample_num = 422 % manuscript figures
% sample_num = 655
% sample_num = 958 % manuscript figures CCB

% Drug dose
% Inibition level for primary effect
drug_dose = 0.50
% TZD effect on vasodilation
drug_dose_vaso = 0
% TZD effect on renin secretion
a = 11/9; b = 1/9;
% drug_dose_rsec = drug_dose + 0.5 % TZD
% drug_dose_rsec = 2*drug_dose
drug_dose_rsec = a * drug_dose ./ (b + drug_dose)
% drug_dose_rsec = 0
% drug_dose_rsec = 1

% Drug dose range
num_dose = 100;
drug_dose_range = linspace(0,0.99,num_dose);

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

for sce_ind = fixed_ss1:fixed_ss1 % scenario
for sex_ind = 1:2        % sex

% Set optional parameters.
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
    if     strcmp(scenario2{fixed_ss2(i)}, 'ACEi' )
            varargin_input = [varargin_input, 'ACEi' ,drug_dose]; % 
    elseif strcmp(scenario2{fixed_ss2(i)}, 'ARB1' )
            varargin_input = [varargin_input, 'ARB1' ,drug_dose]; % 
    elseif strcmp(scenario2{fixed_ss2(i)}, 'CCB'  )
            varargin_input = [varargin_input, 'CCB'  ,[drug_dose,2/3]]; % 
    elseif strcmp(scenario2{fixed_ss2(i)}, 'TZD'  )
        if     strcmp(sex{sex_ind}, 'male')
            varargin_input = [varargin_input, 'TZD'  ,[drug_dose/1.0,drug_dose_vaso,drug_dose_rsec]]; % 
        elseif strcmp(sex{sex_ind}, 'female')
            varargin_input = [varargin_input, 'TZD'  ,[drug_dose/1.0,drug_dose_vaso,drug_dose_rsec]]; % 
        end
    end
end

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

end % sex
end % scenario

%% Load bootstrap replicate parameters & variables before and after drug dose.

% Parameters
load_data_name_pars = sprintf(...
      '%s_male_pars_scenario_Pri_Hyp_bs_rep1000.mat', species{spe_ind});
load(load_data_name_pars, 'pars_rep');
load_data_name_pars = sprintf(...
    '%s_female_pars_scenario_Pri_Hyp_bs_rep1000.mat', species{spe_ind});
load(load_data_name_pars, 'pars_rep');

% Variables after drug dose and hypertensive baseline before drug dose
% X_m/f = (variable, sample, scenario)
load_data_name_vars = sprintf('%s_male_ss_data_scenario_Pri_Hyp_%s.mat'  , ...
                              species{spe_ind},scenario2{fixed_ss2})
load(load_data_name_vars, ...
     'X_ss_m' , 'X_ss_mean_m' , 'X_ss_std_m' , ...
     'X_bl_m' , 'X_bl_mean_m' , 'X_bl_std_m' , ...
     'X_rel_m', 'X_rel_mean_m', 'X_rel_std_m');
load_data_name_vars = sprintf('%s_female_ss_data_scenario_Pri_Hyp_%s.mat', ...
                              species{spe_ind},scenario2{fixed_ss2})
load(load_data_name_vars, ...
     'X_ss_f' , 'X_ss_mean_f' , 'X_ss_std_f' , ...
     'X_bl_f' , 'X_bl_mean_f' , 'X_bl_std_f' , ...
     'X_rel_f', 'X_rel_mean_f', 'X_rel_std_f');

% Compute 95% confidence interval
X_bl_lower_m  = X_bl_mean_m  - 2*X_bl_std_m ; X_bl_upper_m  = X_bl_mean_m  + 2*X_bl_std_m ;
X_ss_lower_m  = X_ss_mean_m  - 2*X_ss_std_m ; X_ss_upper_m  = X_ss_mean_m  + 2*X_ss_std_m ;
X_rel_lower_m = X_rel_mean_m - 2*X_rel_std_m; X_rel_upper_m = X_rel_mean_m + 2*X_rel_std_m;
% ---
X_bl_lower_f  = X_bl_mean_f  - 2*X_bl_std_f ; X_bl_upper_f  = X_bl_mean_f  + 2*X_bl_std_f ;
X_ss_lower_f  = X_ss_mean_f  - 2*X_ss_std_f ; X_ss_upper_f  = X_ss_mean_f  + 2*X_ss_std_f ;
X_rel_lower_f = X_rel_mean_f - 2*X_rel_std_f; X_rel_upper_f = X_rel_mean_f + 2*X_rel_std_f;

%% Post processing

% Retrieve male and female.
% X_m/f = (variables, points, scenario)
t_dy = t_dy';
X_dy_m = reshape(X_dy(:,:,1,:), [num_vars,N,num_scen]); 
X_dy_f = reshape(X_dy(:,:,2,:), [num_vars,N,num_scen]); 

% x-axis limits
xlower = t0; xupper = tend; 

% x-axis
xscale = drug_dose_range * 100;
xlabel_name = strcat(scenario2{fixed_ss2}, ' %');


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

%% Quantities used for calibration/validation for each sex and all scenarios.

% X_m/f = (variable, points, scenario)
BV_m   = reshape(X_dy_m(30,:,:), [N,num_scen]);
BV_f   = reshape(X_dy_f(30,:,:), [N,num_scen]);
R_m    = reshape(X_dy_m(74,:,:) ./ X_dy_m( 4,:,:), [N,num_scen]);
R_f    = reshape(X_dy_f(74,:,:) ./ X_dy_f( 4,:,:), [N,num_scen]);
CO_m   = reshape(X_dy_m(33,:,:), [N,num_scen]);
CO_f   = reshape(X_dy_f(33,:,:), [N,num_scen]);
TPR_m  = reshape(X_dy_m(41,:,:), [N,num_scen]);
TPR_f  = reshape(X_dy_f(41,:,:), [N,num_scen]);
UNA_m  = reshape(X_dy_m(27,:,:), [N,num_scen]);
UNA_f  = reshape(X_dy_f(27,:,:), [N,num_scen]);
UW_m   = reshape(X_dy_m(92,:,:), [N,num_scen]);
UW_f   = reshape(X_dy_f(92,:,:), [N,num_scen]);
CALD_m = reshape(X_dy_m(55,:,:), [N,num_scen]);
CALD_f = reshape(X_dy_f(55,:,:), [N,num_scen]);
RSNA_m = reshape(X_dy_m( 1,:,:), [N,num_scen]);
RSNA_f = reshape(X_dy_f( 1,:,:), [N,num_scen]);
% 
FRNA_m = reshape((X_dy_m(11,:,:) - X_dy_m(27,:,:)) ./ X_dy_m(11,:,:), [N,num_scen]) * 100;
FRNA_f = reshape((X_dy_f(11,:,:) - X_dy_f(27,:,:)) ./ X_dy_f(11,:,:), [N,num_scen]) * 100;
FRW_m  = reshape((X_dy_m( 7,:,:) - X_dy_m(92,:,:)) ./ X_dy_m( 7,:,:), [N,num_scen]) * 100;
FRW_f  = reshape((X_dy_f( 7,:,:) - X_dy_f(92,:,:)) ./ X_dy_f( 7,:,:), [N,num_scen]) * 100;
% 
MAP_m   = reshape(X_dy_m(42,:,:), [N,num_scen]);
MAP_f   = reshape(X_dy_f(42,:,:), [N,num_scen]);
PRA_m   = reshape(X_dy_m(66,:,:), [N,num_scen]);
PRA_f   = reshape(X_dy_f(66,:,:), [N,num_scen]);
ANGI_m  = reshape(X_dy_m(67,:,:), [N,num_scen]);
ANGI_f  = reshape(X_dy_f(67,:,:), [N,num_scen]);
ANGII_m = reshape(X_dy_m(68,:,:), [N,num_scen]);
ANGII_f = reshape(X_dy_f(68,:,:), [N,num_scen]);
CSOD_m  = reshape(X_dy_m(52,:,:), [N,num_scen]);
CSOD_f  = reshape(X_dy_f(52,:,:), [N,num_scen]);
% Plot as relative change in order to compare male and female.
BV_m_bl   = BV_m  (1,:);
BV_f_bl   = BV_f  (1,:);
R_m_bl    = R_m   (1,:);
R_f_bl    = R_f   (1,:);
CO_m_bl   = CO_m  (1,:);
CO_f_bl   = CO_f  (1,:);
TPR_m_bl  = TPR_m (1,:);
TPR_f_bl  = TPR_f (1,:);
UNA_m_bl  = UNA_m (1,:);
UNA_f_bl  = UNA_f (1,:);
UW_m_bl   = UW_m  (1,:);
UW_f_bl   = UW_f  (1,:);
CALD_m_bl = CALD_m(1,:);
CALD_f_bl = CALD_f(1,:);
RSNA_m_bl = RSNA_m(1,:);
RSNA_f_bl = RSNA_f(1,:);
% 
FRNA_m_bl = FRNA_m(1,:);
FRNA_f_bl = FRNA_f(1,:);
FRW_m_bl  = FRW_m (1,:);
FRW_f_bl  = FRW_f (1,:);
% 
MAP_m_bl   = MAP_m  (1,:);
MAP_f_bl   = MAP_f  (1,:);
PRA_m_bl   = PRA_m  (1,:);
PRA_f_bl   = PRA_f  (1,:);
ANGI_m_bl  = ANGI_m (1,:);
ANGI_f_bl  = ANGI_f (1,:);
ANGII_m_bl = ANGII_m(1,:);
ANGII_f_bl = ANGII_f(1,:);
CSOD_m_bl  = CSOD_m (1,:);
CSOD_f_bl  = CSOD_f (1,:);
for i = 1:N
    BV_m  (i,:) = BV_m  (i,:) ./ BV_m_bl  ;
    BV_f  (i,:) = BV_f  (i,:) ./ BV_f_bl  ;
    R_m   (i,:) = R_m   (i,:) ./ R_m_bl   ;
    R_f   (i,:) = R_f   (i,:) ./ R_f_bl   ;
    CO_m  (i,:) = CO_m  (i,:) ./ CO_m_bl  ;
    CO_f  (i,:) = CO_f  (i,:) ./ CO_f_bl  ;
    TPR_m (i,:) = TPR_m (i,:) ./ TPR_m_bl ;
    TPR_f (i,:) = TPR_f (i,:) ./ TPR_f_bl ;
    UNA_m (i,:) = UNA_m (i,:) ./ UNA_m_bl ;
    UNA_f (i,:) = UNA_f (i,:) ./ UNA_f_bl ;
    UW_m  (i,:) = UW_m  (i,:) ./ UW_m_bl  ;
    UW_f  (i,:) = UW_f  (i,:) ./ UW_f_bl  ;
    CALD_m(i,:) = CALD_m(i,:) ./ CALD_m_bl;
    CALD_f(i,:) = CALD_f(i,:) ./ CALD_f_bl;
    RSNA_m(i,:) = RSNA_m(i,:) ./ RSNA_m_bl;
    RSNA_f(i,:) = RSNA_f(i,:) ./ RSNA_f_bl;
%     
    FRNA_m(i,:) = FRNA_m(i,:) ./ FRNA_m_bl;
    FRNA_f(i,:) = FRNA_f(i,:) ./ FRNA_f_bl;
    FRW_m (i,:) = FRW_m (i,:) ./ FRW_m_bl ;
    FRW_f (i,:) = FRW_f (i,:) ./ FRW_f_bl ;
% 
    MAP_m  (i,:) = MAP_m  (i,:) ./ MAP_m_bl  ;
    MAP_f  (i,:) = MAP_f  (i,:) ./ MAP_f_bl  ;
    PRA_m  (i,:) = PRA_m  (i,:) ./ PRA_m_bl  ;
    PRA_f  (i,:) = PRA_f  (i,:) ./ PRA_f_bl  ;
    ANGI_m (i,:) = ANGI_m (i,:) ./ ANGI_m_bl ;
    ANGI_f (i,:) = ANGI_f (i,:) ./ ANGI_f_bl ;
    ANGII_m(i,:) = ANGII_m(i,:) ./ ANGII_m_bl;
    ANGII_f(i,:) = ANGII_f(i,:) ./ ANGII_f_bl;
    CSOD_m(i,:)  = CSOD_m (i,:) ./ CSOD_m_bl ;
    CSOD_f(i,:)  = CSOD_f (i,:) ./ CSOD_f_bl ;
end
%%

% Actual and % change in MAP. 
MAP_ss_mean_m   = X_ss_mean_m (42,:);
MAP_ss_lower_m  = X_ss_lower_m(42,:);
MAP_ss_upper_m  = X_ss_upper_m(42,:);
MAP_rel_mean_m  = X_rel_mean_m (42,:);
MAP_rel_lower_m = X_rel_lower_m(42,:);
MAP_rel_upper_m = X_rel_upper_m(42,:);
% ---
MAP_ss_mean_f   = X_ss_mean_f (42,:);
MAP_ss_lower_f  = X_ss_lower_f(42,:);
MAP_ss_upper_f  = X_ss_upper_f(42,:);
MAP_rel_mean_f  = X_rel_mean_f (42,:);
MAP_rel_lower_f = X_rel_lower_f(42,:);
MAP_rel_upper_f = X_rel_upper_f(42,:);

% y-axis limits
ylower_MAP = min( min(MAP_rel_lower_m),min(MAP_rel_lower_f) );

%% Plot quantities for manuscript. ----------------------------------------

man_f(1) = figure('DefaultAxesFontSize',18);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5.0, 3.5]);
plot(t_dy,R_m   (:,fixed_ss1), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim([xlower, xupper]);
xticks([tchange_dy+0*(1) : 1 : tchange_dy+days*(1)]);
xticklabels({'0','1','2','3','4','5','6','7'});
% ylim(s_main(1), [0.75,1.05])
xlabel('Time (days)'); ylabel('EAR/RVR (relative)');
hold on
plot(t_dy,R_f   (:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold off
legend('Male','Female', 'FontSize',9,'Location','Southeast');
if     strcmp(scenario2{fixed_ss2}, 'ACEi')
    title('C', 'FontWeight','normal')
elseif strcmp(scenario2{fixed_ss2}, 'ARB1')
    title('C', 'FontWeight','normal')
elseif strcmp(scenario2{fixed_ss2}, 'CCB' )
    title('B', 'FontWeight','normal')
elseif strcmp(scenario2{fixed_ss2}, 'TZD' )
    title('-', 'FontWeight','normal')
end

man_f(2) = figure('DefaultAxesFontSize',18);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5.0, 3.5]);
plot(t_dy,FRNA_m(:,fixed_ss1) ,'-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
xlim([xlower, xupper]);
xticks([tchange_dy+0*(1) : 1 : tchange_dy+days*(1)]);
xticklabels({'0','1','2','3','4','5','6','7','8'});
% ylim(s_main(2), [97,100])
xlabel('Time (days)'); ylabel('FR (relative)');
hold on
plot(t_dy,FRW_m (:,fixed_ss1), '--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(t_dy,FRNA_f(:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(t_dy,FRW_f (:,fixed_ss1), '--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
fakeplot = zeros(2, 1);
fakeplot(1) = plot(NaN,NaN, 'k-' );
fakeplot(2) = plot(NaN,NaN, 'k--');
[~, hobj, ~, ~] = legend(fakeplot, {'FR_{Na^+}','FR_{W}'}, 'FontSize',9,'Location','Southeast');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
hold off
if     strcmp(scenario2{fixed_ss2}, 'ACEi')
    title('E', 'FontWeight','normal')
elseif strcmp(scenario2{fixed_ss2}, 'ARB1')
    title('E', 'FontWeight','normal')
elseif strcmp(scenario2{fixed_ss2}, 'CCB' )
    title('D', 'FontWeight','normal')
elseif strcmp(scenario2{fixed_ss2}, 'TZD' )
    title('B', 'FontWeight','normal')
end

man_f(3) = figure('DefaultAxesFontSize',18);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5.0, 3.5]);
plot(t_dy,BV_m   (:,fixed_ss1), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim([xlower, xupper]);
xticks([tchange_dy+0*(1) : 1 : tchange_dy+days*(1)]);
xticklabels({'0','1','2','3','4','5','6','7'});
% ylim(s_main(1), [0.75,1.05])
xlabel('Time (days)'); ylabel('BV (relative)');
hold on
plot(t_dy,BV_f   (:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold off
legend('Male','Female', 'FontSize',9,'Location','Northeast');
if     strcmp(scenario2{fixed_ss2}, 'ACEi')
    title('F', 'FontWeight','normal')
elseif strcmp(scenario2{fixed_ss2}, 'ARB1')
    title('F', 'FontWeight','normal')
elseif strcmp(scenario2{fixed_ss2}, 'CCB' )
    title('-', 'FontWeight','normal')
elseif strcmp(scenario2{fixed_ss2}, 'TZD' )
    title('-', 'FontWeight','normal')
end

man_f(4) = figure('DefaultAxesFontSize',18);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5.0, 3.5]);
plot(t_dy,TPR_m   (:,fixed_ss1), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim([xlower, xupper]);
xticks([tchange_dy+0*(1) : 1 : tchange_dy+days*(1)]);
xticklabels({'0','1','2','3','4','5','6','7'});
% ylim(s_main(1), [0.75,1.05])
xlabel('Time (days)'); ylabel('TPR (relative)');
hold on
plot(t_dy,TPR_f   (:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold off
legend('Male','Female', 'FontSize',9,'Location','Northeast');
if     strcmp(scenario2{fixed_ss2}, 'ACEi')
    title('-', 'FontWeight','normal')
elseif strcmp(scenario2{fixed_ss2}, 'ARB1')
    title('-', 'FontWeight','normal')
elseif strcmp(scenario2{fixed_ss2}, 'CCB' )
    title('F', 'FontWeight','normal')
elseif strcmp(scenario2{fixed_ss2}, 'TZD' )
    title('E', 'FontWeight','normal')
end

man_f(5) = figure('DefaultAxesFontSize',18);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5.0, 3.5]);
plot(t_dy,UNA_m (:,fixed_ss1) ,'-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
xlim([xlower, xupper]);
xticks([tchange_dy+0*(1) : 1 : tchange_dy+days*(1)]);
xticklabels({'0','1','2','3','4','5','6','7','8'});
% ylim(s_main(2), [97,100])
xlabel('Time (days)'); ylabel('UF (relative)');
hold on
plot(t_dy,UW_m  (:,fixed_ss1), '--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(t_dy,UNA_f (:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(t_dy,UW_f  (:,fixed_ss1), '--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
fakeplot = zeros(2, 1);
fakeplot(1) = plot(NaN,NaN, 'k-' );
fakeplot(2) = plot(NaN,NaN, 'k--');
[~, hobj, ~, ~] = legend(fakeplot, {'UF_{Na^+}','UF_{W}'}, 'FontSize',9,'Location','Northeast');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
hold off
if     strcmp(scenario2{fixed_ss2}, 'ACEi')
    title('-', 'FontWeight','normal')
elseif strcmp(scenario2{fixed_ss2}, 'ARB1')
    title('-', 'FontWeight','normal')
elseif strcmp(scenario2{fixed_ss2}, 'CCB' )
    title('-', 'FontWeight','normal')
elseif strcmp(scenario2{fixed_ss2}, 'TZD' )
    title('C', 'FontWeight','normal')
end

man_f(6) = figure('DefaultAxesFontSize',18);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5.0, 3.5]);
plot(t_dy,CALD_m   (:,fixed_ss1), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim([xlower, xupper]);
xticks([tchange_dy+0*(1) : 1 : tchange_dy+days*(1)]);
xticklabels({'0','1','2','3','4','5','6','7'});
% ylim(s_main(1), [0.75,1.05])
xlabel('Time (days)'); ylabel('ALD (relative)');
hold on
plot(t_dy,CALD_f   (:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold off
legend('Male','Female', 'FontSize',9,'Location','Northeast');
if     strcmp(scenario2{fixed_ss2}, 'ACEi')
    title('D', 'FontWeight','normal')
elseif strcmp(scenario2{fixed_ss2}, 'ARB1')
    title('D', 'FontWeight','normal')
elseif strcmp(scenario2{fixed_ss2}, 'CCB' )
    title('-', 'FontWeight','normal')
elseif strcmp(scenario2{fixed_ss2}, 'TZD' )
    title('-', 'FontWeight','normal')
end

man_f(7) = figure('DefaultAxesFontSize',18);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5.0, 3.5]);
plot(t_dy,RSNA_m   (:,fixed_ss1), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim([xlower, xupper]);
xticks([tchange_dy+0*(1) : 1 : tchange_dy+days*(1)]);
xticklabels({'0','1','2','3','4','5','6','7'});
% ylim(s_main(1), [0.75,1.05])
xlabel('Time (days)'); ylabel('RSNA (relative)');
hold on
plot(t_dy,RSNA_f   (:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold off
legend('Male','Female', 'FontSize',9,'Location','Northeast');
if     strcmp(scenario2{fixed_ss2}, 'ACEi')
    title('-', 'FontWeight','normal')
elseif strcmp(scenario2{fixed_ss2}, 'ARB1')
    title('-', 'FontWeight','normal')
elseif strcmp(scenario2{fixed_ss2}, 'CCB' )
    title('C', 'FontWeight','normal')
elseif strcmp(scenario2{fixed_ss2}, 'TZD' )
    title('-', 'FontWeight','normal')
end

man_f(8) = figure('DefaultAxesFontSize',18);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5.0, 3.5]);
plot(t_dy,CO_m   (:,fixed_ss1), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim([xlower, xupper]);
xticks([tchange_dy+0*(1) : 1 : tchange_dy+days*(1)]);
xticklabels({'0','1','2','3','4','5','6','7'});
% ylim(s_main(1), [0.75,1.05])
xlabel('Time (days)'); ylabel('CO (relative)');
hold on
plot(t_dy,CO_f   (:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold off
legend('Male','Female', 'FontSize',9,'Location','Northeast');
if     strcmp(scenario2{fixed_ss2}, 'ACEi')
    title('-', 'FontWeight','normal')
elseif strcmp(scenario2{fixed_ss2}, 'ARB1')
    title('-', 'FontWeight','normal')
elseif strcmp(scenario2{fixed_ss2}, 'CCB' )
    title('E', 'FontWeight','normal')
elseif strcmp(scenario2{fixed_ss2}, 'TZD' )
    title('D', 'FontWeight','normal')
end

% figure('DefaultAxesFontSize',18);
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5.0, 3.5]);
% plot(t_dy,MAP_m   (:,fixed_ss1), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
% xlim([xlower, xupper]);
% xticks([tchange_dy+0*(1) : 1 : tchange_dy+days*(1)]);
% xticklabels({'0','1','2','3','4','5','6','7'});
% % ylim(s_main(1), [0.75,1.05])
% xlabel('Time (days)'); ylabel('MAP (relative)');
% hold on
% plot(t_dy,MAP_f   (:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
% hold off
% legend('Male','Female', 'FontSize',9,'Location','Northeast');
% if     strcmp(scenario2{fixed_ss2}, 'ACEi')
%     title('-', 'FontWeight','normal')
% elseif strcmp(scenario2{fixed_ss2}, 'ARB1')
%     title('-', 'FontWeight','normal')
% elseif strcmp(scenario2{fixed_ss2}, 'CCB' )
%     title('-', 'FontWeight','normal')
% elseif strcmp(scenario2{fixed_ss2}, 'TZD' )
%     title('-', 'FontWeight','normal')
% end

man_f(9) = figure('DefaultAxesFontSize',18);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5.0, 3.5]);
plot(xscale,MAP_rel_mean_m , '-', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
% xlim([xlower, xupper]);
xticks(0:20:100);
% ylim([ylower_MAP, 0]);
xlabel(xlabel_name); ylabel('% \DeltaMAP');
hold on
plot(xscale,MAP_rel_mean_f , '-', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(xscale,MAP_rel_lower_m, ':', 'Color',[0.203, 0.592, 0.835], 'LineWidth',1, 'MarkerSize',8);
plot(xscale,MAP_rel_lower_f, ':', 'Color',[0.835, 0.203, 0.576], 'LineWidth',1, 'MarkerSize',8);
plot(xscale,MAP_rel_upper_m, ':', 'Color',[0.203, 0.592, 0.835], 'LineWidth',1, 'MarkerSize',8);
plot(xscale,MAP_rel_upper_f, ':', 'Color',[0.835, 0.203, 0.576], 'LineWidth',1, 'MarkerSize',8);
hold off
legend('Male','Female', 'FontSize',9,'Location','Southwest');
if     strcmp(scenario2{fixed_ss2}, 'ACEi')
    title('A', 'FontWeight','normal')
elseif strcmp(scenario2{fixed_ss2}, 'ARB1')
    title('B', 'FontWeight','normal')
elseif strcmp(scenario2{fixed_ss2}, 'CCB' )
    title('A', 'FontWeight','normal')
elseif strcmp(scenario2{fixed_ss2}, 'TZD' )
    title('A', 'FontWeight','normal')
end

%% Save figures. ----------------------------------------------------------

save_data_name = sprintf('man_figs_%s%s%%_VI%s.fig', ...
                         scenario2{fixed_ss2},num2str(drug_dose*100),num2str(sample_num));
save_data_name = strcat('Figures/', save_data_name);
savefig([man_f], save_data_name)
% ---
fig_names = {'R','FR','BV','TPR','UF','ALD','RSNA','CO','MAP'};
for i = 1:length(man_f)
save_data_name = sprintf('%s_%s.png', ...
                         scenario2{fixed_ss2},fig_names{i});
save_data_name = strcat('Figures/', save_data_name);
exportgraphics([man_f(i)], save_data_name)
end

% save_data_name = sprintf('man_figs_%s%s%%_VI%s.png', ...
%                          scenario2{fixed_ss2},num2str(drug_dose*100),num2str(sample_num));
% save_data_name = strcat('Figures/', save_data_name);
% exportgraphics([man_f(5)], save_data_name)


end






























