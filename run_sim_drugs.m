% This simulates the blood pressure regulation model bp_reg.m for drug administration.
% 
% Steady state data is calculated by solve_ss_scenario.m.

function run_sim_drugs

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
% m_RAS   - male RAS pars
% m_Reab  - male fractional sodium and water reabsorption
% Pri_Hyp - essential/primary hypertension
scenario1 = {'Normal', 'm_RSNA', 'm_AT2R', 'm_RAS', 'm_Reab', ...
             'm_RAS_m_Reab', 'm_RSNA_m_Reab'};
% scenario1 = {'Normal', 'm_RSNA', 'm_AT2R', 'm_RAS', 'm_Reab', ...
%              'm_RAS_m_Reab', 'm_RSNA_m_Reab', ...
%              'Pri_Hyp'};
fixed_ss1 = 1;
num_scen = length(scenario1);
% Drug scenarios
% Normal - Normal conditions
% ACEi   - Angiotensin converting enzyme inhibitor % 95
% ARB1   - Angiotensin receptor 1 blocker % 94
% CCB    - Calcium channel blocker % 84
% DIU    - Thiazide diuretic % 0.5 1?
% ARB2   - Angiotensin receptor 2 blocker %
% DRI    - Direct renin inhibitor %
% MRB    - Aldosterone blocker (MR?) %
% RSS    - Renin secretion stimulator (thiazide?) % % NOT COMPLETE
% AngII  - Ang II infusion fmol/(ml min)
scenario2 = {'Normal', 'ACEi', 'ARB1', 'CCB', 'DIU', ...
             'ARB2'  , 'DRI' , 'MRB' , 'RSS', 'AngII'};
fixed_ss2 = [4];

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
% sample_num = 208
sample_num = 655

% Drug dose
drug_dose = 0.60
drug_dose_vaso = 0           % DIU
% a = 3; b = 1;
a = 11/9; b = 1/9;
% drug_dose_rsec = drug_dose + 0.5 % DIU
% drug_dose_rsec = 2*drug_dose
drug_dose_rsec = a * drug_dose ./ (b + drug_dose)
% drug_dose_rsec = 0
% drug_dose_rsec = 1

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
load_data_name_pars = sprintf('%s_%s_pars_scenario_Pri_Hyp_bs_rep1000OLD.mat', ...
                              species{spe_ind},sex{sex_ind});
load(load_data_name_pars, 'pars_rep');
num_pars   = size(pars_rep,1);
num_sample = size(pars_rep,2);
pars_rep = pars_rep(:,sample_num);

% Variables
load_data_name_vars = sprintf('%s_%s_ss_data_scenario_Pri_Hyp_bs_rep1000OLD.mat', ...
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
    elseif strcmp(scenario2{fixed_ss2(i)}, 'DIU'  )
        if     strcmp(sex{sex_ind}, 'male')
            varargin_input = [varargin_input, 'DIU'  ,[drug_dose/1.0,drug_dose_vaso,drug_dose_rsec]]; % 
        elseif strcmp(sex{sex_ind}, 'female')
            varargin_input = [varargin_input, 'DIU'  ,[drug_dose/1.0,drug_dose_vaso,drug_dose_rsec]]; % 
        end
    elseif strcmp(scenario2{fixed_ss2(i)}, 'ARB2' )
            varargin_input = [varargin_input, 'ARB2' ,drug_dose]; % 
    elseif strcmp(scenario2{fixed_ss2(i)}, 'DRI'  )
            varargin_input = [varargin_input, 'DRI'  ,drug_dose]; % 
    elseif strcmp(scenario2{fixed_ss2(i)}, 'MRB'  )
            varargin_input = [varargin_input, 'MRB'  ,drug_dose]; % 
    elseif strcmp(scenario2{fixed_ss2(i)}, 'RSS'  )
            varargin_input = [varargin_input, 'RSS'  ,drug_dose]; % 

    elseif strcmp(scenario2{fixed_ss2(i)}, 'AngII')
        if     strcmp(sex{sex_ind}, 'male')
            varargin_input = [varargin_input, 'AngII',910]; % Sullivan 2010
        elseif strcmp(sex{sex_ind}, 'female')
            varargin_input = [varargin_input, 'AngII',505]; % Sullivan 2010
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

%% Post processing

% Retrieve male and female.
% X_m/f = (variables, points, scenario)
t_dy = t_dy';
X_dy_m = reshape(X_dy(:,:,1,:), [num_vars,N,num_scen]); 
X_dy_f = reshape(X_dy(:,:,2,:), [num_vars,N,num_scen]); 
% % X_m/f = (variables, scenario)
% X_ss_m = reshape(X_ss(:,  1,:), [num_vars,  num_scen]); 
% X_ss_f = reshape(X_ss(:,  2,:), [num_vars,  num_scen]); 

% Total urine and sodium excretion
% 27, 92
% t_dy
% (t_dy(2:102) - t_dy(1:101))'
% X_dy_m(92,6,fixed_ss1)
% Baseline and drug 4 and 24 hour excretion
% U12_b_m  = sum( X_dy_m(92,  1:  1+10,fixed_ss1) .* (t_dy(  2:  2+10) - t_dy(  1:  1+10)) );
% U12_d_m  = sum( X_dy_m(92,101:101+10,fixed_ss1) .* (t_dy(102:102+10) - t_dy(101:101+10)) );
% U12_d_m/U12_b_m
% S12_b_m  = sum( X_dy_m(27,  1:  1+10,fixed_ss1) .* (t_dy(  2:  2+10) - t_dy(  1:  1+10)) );
% S12_d_m  = sum( X_dy_m(27,101:101+10,fixed_ss1) .* (t_dy(102:102+10) - t_dy(101:101+10)) );
% S12_d_m/S12_b_m
% U12_b_f  = sum( X_dy_f(92,  1:  1+10,fixed_ss1) .* (t_dy(  2:  2+10) - t_dy(  1:  1+10)) );
% U12_d_f  = sum( X_dy_f(92,101:101+10,fixed_ss1) .* (t_dy(102:102+10) - t_dy(101:101+10)) );
% U12_d_f/U12_b_f
% S12_b_f  = sum( X_dy_f(27,  1:  1+10,fixed_ss1) .* (t_dy(  2:  2+10) - t_dy(  1:  1+10)) );
% S12_d_f  = sum( X_dy_f(27,101:101+10,fixed_ss1) .* (t_dy(102:102+10) - t_dy(101:101+10)) );
% S12_d_f/S12_b_f
% ---
% U4_b_m  = sum( X_dy_m(92,  1: 18,fixed_ss1) .* (t_dy(  2: 19) - t_dy(  1: 18)) )
% U4_d_m  = sum( X_dy_m(92,101:118,fixed_ss1) .* (t_dy(102:119) - t_dy(101:118)) )
% U4_d_m/U4_b_m
% S4_b_m  = sum( X_dy_m(27,  1: 18,fixed_ss1) .* (t_dy(  2: 19) - t_dy(  1: 18)) )
% S4_d_m  = sum( X_dy_m(27,101:118,fixed_ss1) .* (t_dy(102:119) - t_dy(101:118)) )
% S4_d_m/S4_b_m
% ---
% U5_b_m  = sum( X_dy_m(92,  1:  1+21,fixed_ss1) .* (t_dy(  2:  2+21) - t_dy(  1:  1+21)) )
% U5_d_m  = sum( X_dy_m(92,101:101+21,fixed_ss1) .* (t_dy(102:102+21) - t_dy(101:101+21)) )
% U5_d_m/U5_b_m
% S5_b_m  = sum( X_dy_m(27,  1:  1+21,fixed_ss1) .* (t_dy(  2:  2+21) - t_dy(  1:  1+21)) )
% S5_d_m  = sum( X_dy_m(27,101:101+21,fixed_ss1) .* (t_dy(102:102+21) - t_dy(101:101+21)) )
% S5_d_m/S5_b_m
% ---
% U24_b_m = sum( X_dy_m(92,  1:  1+100,fixed_ss1) .* (t_dy(  2:  2+100) - t_dy(  1:  1+100)) )
% U24_d_m = sum( X_dy_m(92,101:101+100,fixed_ss1) .* (t_dy(102:102+100) - t_dy(101:101+100)) )
% U24_d_m/U24_b_m
% S24_b_m = sum( X_dy_m(27,  1:  1+100,fixed_ss1) .* (t_dy(  2:  2+100) - t_dy(  1:  1+100)) )
% S24_d_m = sum( X_dy_m(27,101:101+100,fixed_ss1) .* (t_dy(102:102+100) - t_dy(101:101+100)) )
% S24_d_m/S24_b_m
% U24_b_f = sum( X_dy_f(92,  1:  1+100,fixed_ss1) .* (t_dy(  2:  2+100) - t_dy(  1:  1+100)) )
% U24_d_f = sum( X_dy_f(92,101:101+100,fixed_ss1) .* (t_dy(102:102+100) - t_dy(101:101+100)) )
% U24_d_f/U24_b_f
% S24_b_f = sum( X_dy_f(27,  1:  1+100,fixed_ss1) .* (t_dy(  2:  2+100) - t_dy(  1:  1+100)) )
% S24_d_f = sum( X_dy_f(27,101:101+100,fixed_ss1) .* (t_dy(102:102+100) - t_dy(101:101+100)) )
% S24_d_f/S24_b_f

%% Plot all vars vs time. -------------------------------------------------

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
        
        plot(s1(i,j), t_dy,X_dy_m((i-1)*15 + j,:,fixed_ss1),'b', ...
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

f2 = figure('pos',[000 000 600 600], 'DefaultAxesFontSize',14);
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

%% Plot quantities that explain mechanisms. -------------------------------

% BV; RSNA; REA/RR for each sex and all scenarios.
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
% Plot as relative change in order to compare male and female.
BV_m_bl  = BV_m (1,:);
BV_f_bl  = BV_f (1,:);
R_m_bl   = R_m  (1,:);
R_f_bl   = R_f  (1,:);
CO_m_bl  = CO_m (1,:);
CO_f_bl  = CO_f (1,:);
TPR_m_bl = TPR_m(1,:);
TPR_f_bl = TPR_f(1,:);
UNA_m_bl = UNA_m(1,:);
UNA_f_bl = UNA_f(1,:);
UW_m_bl  = UW_m (1,:);
UW_f_bl  = UW_f (1,:);
for i = 1:N
    BV_m (i,:) = BV_m (i,:) ./ BV_m_bl ;
    BV_f (i,:) = BV_f (i,:) ./ BV_f_bl ;
    R_m  (i,:) = R_m  (i,:) ./ R_m_bl  ;
    R_f  (i,:) = R_f  (i,:) ./ R_f_bl  ;
    CO_m (i,:) = CO_m (i,:) ./ CO_m_bl ;
    CO_f (i,:) = CO_f (i,:) ./ CO_f_bl ;
    TPR_m(i,:) = TPR_m(i,:) ./ TPR_m_bl;
    TPR_f(i,:) = TPR_f(i,:) ./ TPR_f_bl;
    UNA_m(i,:) = UNA_m(i,:) ./ UNA_m_bl;
    UNA_f(i,:) = UNA_f(i,:) ./ UNA_f_bl;
    UW_m (i,:) = UW_m (i,:) ./ UW_m_bl ;
    UW_f (i,:) = UW_f (i,:) ./ UW_f_bl ;
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

% ACEi, ARB
g1 = figure('DefaultAxesFontSize',14);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 3]);
t1 = tiledlayout(1,3,'TileSpacing','Normal','Padding','Compact');

nexttile
plot(t_dy,R_m   (:,fixed_ss1), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
xlim([xlower+0.75, xupper]);
xticks([tchange_dy+0*(1) : 1 : tchange_dy+days*(1)]);
xticklabels({'0','1','2','3','4','5','6','7','8'});
xlabel('Time (days)'); ylabel('EAR/RVR (relative)');
hold on
plot(t_dy,R_f   (:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold off
[~, hobj, ~, ~] = legend('Male','Female', 'FontSize',7, 'Location','Southeast');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
title('A')

nexttile
plot(t_dy,FRNA_m(:,fixed_ss1) ,'-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
xlim([xlower+0.75, xupper]);
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
[~, hobj, ~, ~] = legend(fakeplot, {'FR_{Na^+}','FR_{W}'}, 'FontSize',7,'Location','Southeast');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
hold off
title('B')

nexttile
plot(t_dy,BV_m  (:,fixed_ss1), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
xlim([xlower+0.75, xupper]);
xticks([tchange_dy+0*(1) : 1 : tchange_dy+days*(1)]);
xticklabels({'0','1','2','3','4','5','6','7','8'});
% ylim(s_main(4), [1,1.25])
xlabel('Time (days)'); ylabel('BV (relative)');
hold on
plot(t_dy,BV_f  (:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold off
title('C')

% CCB
g2 = figure('DefaultAxesFontSize',14);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 3]);
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10*2/3, 3]);
t2 = tiledlayout(1,3,'TileSpacing','Normal','Padding','Compact');

nexttile
plot(t_dy,R_m   (:,fixed_ss1), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
xlim([xlower+0.75, xupper]);
xticks([tchange_dy+0*(1) : 1 : tchange_dy+days*(1)]);
xticklabels({'0','1','2','3','4','5','6','7','8'});
xlabel('Time (days)'); ylabel('EAR/RVR (relative)');
hold on
plot(t_dy,R_f   (:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold off
[~, hobj, ~, ~] = legend('Male','Female', 'FontSize',7, 'Location','Southeast');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
title('A')

nexttile
plot(t_dy,FRNA_m(:,fixed_ss1) ,'-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
xlim([xlower+0.75, xupper]);
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
[~, hobj, ~, ~] = legend(fakeplot, {'FR_{Na^+}','FR_{W}'}, 'FontSize',7,'Location','Southeast');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
hold off
title('B')

nexttile
plot(t_dy,TPR_m  (:,fixed_ss1), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
xlim([xlower+0.75, xupper]);
xticks([tchange_dy+0*(1) : 1 : tchange_dy+days*(1)]);
xticklabels({'0','1','2','3','4','5','6','7','8'});
% ylim(s_main(4), [1,1.25])
xlabel('Time (days)'); ylabel('TPR (relative)');
hold on
plot(t_dy,TPR_f  (:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold off
title('C')

% DIU
g3 = figure('DefaultAxesFontSize',14);
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 3]);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10*3/3, 3*2]);
t3 = tiledlayout(2,2,'TileSpacing','Normal','Padding','Compact');

nexttile
plot(t_dy,FRNA_m(:,fixed_ss1) ,'-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
xlim([xlower+0.75, xupper]);
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
[~, hobj, ~, ~] = legend(fakeplot, {'FR_{Na^+}','FR_{W}'}, 'FontSize',7,'Location','Southeast');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
hold off
title('A')

nexttile
plot(t_dy,UNA_m (:,fixed_ss1) ,'-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
xlim([xlower+0.75, xupper]);
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
[~, hobj, ~, ~] = legend(fakeplot, {'UF_{Na^+}','UF_{W}'}, 'FontSize',7,'Location','Northeast');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
hold off
title('B')

nexttile
plot(t_dy,CO_m   (:,fixed_ss1), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
xlim([xlower+0.75, xupper]);
xticks([tchange_dy+0*(1) : 1 : tchange_dy+days*(1)]);
xticklabels({'0','1','2','3','4','5','6','7','8'});
% ylim(s_main(4), [1,1.25])
xlabel('Time (days)'); ylabel('CO (relative)');
hold on
plot(t_dy,CO_f  (:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold off
title('C')

nexttile
plot(t_dy,TPR_m (:,fixed_ss1), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
xlim([xlower+0.75, xupper]);
xticks([tchange_dy+0*(1) : 1 : tchange_dy+days*(1)]);
xticklabels({'0','1','2','3','4','5','6','7','8'});
% ylim(s_main(4), [1,1.25])
xlabel('Time (days)'); ylabel('TPR (relative)');
hold on
plot(t_dy,TPR_f (:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold off
title('D')

%% Plot quantities for calibration/validation. ----------------------------

% Calibration/validation quantities for each sex and all scenarios.
% X_m/f = (variable, points, scenario)
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

h1 = figure('DefaultAxesFontSize',14);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.15, 5]);
t = tiledlayout(2,3,'TileSpacing','Normal','Padding','Compact');

nexttile
plot(t_dy,MAP_m   (:,fixed_ss1), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim([xlower, xupper]);
xticks([tchange_dy+0*(1) : 2 : tchange_dy+days*(1)]);
xticklabels({'0','2','4','6','8','10','12','14'});
% ylim(s_main(1), [0.75,1.05])
xlabel('Time (days)'); ylabel('MAP (relative)');
hold on
plot(t_dy,MAP_f   (:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold off
legend('Male','Female', 'FontSize',7,'Location','Northeast');
title('A')

nexttile
plot(t_dy,PRA_m(:,fixed_ss1) ,'-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim([xlower, xupper]);
xticks([tchange_dy+0*(1) : 2 : tchange_dy+days*(1)]);
xticklabels({'0','2','4','6','8','10','12','14'});
% ylim(s_main(2), [97,100])
xlabel('Time (days)'); ylabel('PRA (relative)');
hold on
plot(t_dy,PRA_f (:,fixed_ss1), '-', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold off
title('B')

nexttile
plot(t_dy,ANGI_m (:,fixed_ss1), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim([xlower, xupper]);
xticks([tchange_dy+0*(1) : 2 : tchange_dy+days*(1)]);
xticklabels({'0','2','4','6','8','10','12','14'});
% ylim(s_main(3), [0.75,1.35])
xlabel('Time (days)'); ylabel('Ang I (relative)');
hold on
plot(t_dy,ANGI_f (:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold off
title('C')

nexttile
plot(t_dy,ANGII_m  (:,fixed_ss1), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim([xlower, xupper]);
xticks([tchange_dy+0*(1) : 2 : tchange_dy+days*(1)]);
xticklabels({'0','2','4','6','8','10','12','14'});
% ylim(s_main(4), [1,1.25])
xlabel('Time (days)'); ylabel('Ang II (relative)');
hold on
plot(t_dy,ANGII_f  (:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold off
title('D')

nexttile
plot(t_dy,CSOD_m   (:,fixed_ss1), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim([xlower, xupper]);
xticks([tchange_dy+0*(1) : 2 : tchange_dy+days*(1)]);
xticklabels({'0','2','4','6','8','10','12','14'});
% ylim(s_main(4), [1,1.25])
xlabel('Time (days)'); ylabel('[Na^{+}] (relative)');
hold on
plot(t_dy,CSOD_f   (:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold off
title('E')

%% Save figures. ----------------------------------------------------------

% save_data_name = sprintf('Pri_hyp_sim_%s%s%%_VI%s.fig', ...
%                          scenario2{fixed_ss2},num2str(drug_dose*100),num2str(sample_num));
% save_data_name = strcat('Figures/', save_data_name);
% savefig([f1;f2;g1;g2;g3;h1], save_data_name)
% ---
% save_data_name = sprintf('Pri_hyp_sim_all_vars_%s%s%%_VI%s.fig', ...
%                          scenario2{fixed_ss2},num2str(drug_dose*100),num2str(sample_num));
% save_data_name = strcat('Figures/', save_data_name);
% savefig([f1], save_data_name)

end






























