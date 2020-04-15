% 

function run_sim_AngII_test_hyp_pars

close all

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Begin user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters to perturb
% K_bar            - par 13; - 0%, +200%
% R_bv             - par 14; - 0%, +200%
% % C_gcf            - par 8 ; -20%
% R_aass           - par 4 ; - 0%, +100%
% N_rs             - par 21; - 0%, +100%
% N_als_eq         - par 18; - 0%, +100%
% N_rsna           - par 3 ; - 0%, +100%
% Indices
par_ind = [13;14;4;21;18;3];
par_num = length(par_ind);
% Range for parameters
par_range_lower = [0  ;0  ;0  ;0  ;0  ;0  ]/100;
par_range_upper = [100;400;100;100;100;100]/100;

% Scenarios
% Normal - Normal conditions
% m_RSNA - male RSNA
% m_AT2R - male AT2R
% m_RAS  - male RAS pars
% m_Reab - male fractional sodium and water reabsorption
% m_RAS_&_m_Reab - male RAS pars & fractional sodium and water reabsorption
% scenario = {'Normal', 'm_RSNA', 'm_AT2R', 'm_RAS', 'm_Reab', ...
%             'm_RSNA_m_Reab'};
scenario = {'Normal', 'm_RSNA', 'm_AT2R', 'm_RAS', 'm_Reab', ...
            'm_RSNA_m_Reab', ...
            'Pri_Hyp'};
num_scen = length(scenario);
fixed_ss = 7;

% Species
spe_ind = 2;

% Number of days to run simulation after change; Day at which to induce change;
days = 14; day_change = 1;
% Number of points for plotting resolution
% N = ((days+1)*1440) / 2;
N = (days+1)*10 + 1;


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

for sce_ind = fixed_ss:fixed_ss % scenario
for sex_ind = 1:2        % sex

varargin_input = {scenario{1},true};

% Parameter input
pars = get_pars(species{spe_ind}, sex{sex_ind}, varargin_input{:});

%% Uniformly randomly sample perturbed parameters from range.

% Set interval bounds for parameters.
lower = pars(par_ind) - par_range_lower .* pars(par_ind);
upper = pars(par_ind) + par_range_upper .* pars(par_ind);

% % Sample between 0 and 1.
% % rng('default')
% ran_vec = random('unif',0,1,par_num,1);
% % Sample within interval.
% ran_vec = lower + ran_vec .* (upper - lower);
% % Replace input parameters with newly sampled parameter.
% pars(par_ind) = ran_vec;
% % pars(par_ind) = [30.2589; 12.5101; 3.6061; 1.19  ; 1.50  ; 1.12  ]; % diverges for j = 15
% % pars(par_ind) = [17.9423; 13.9705; 5.6973; 1.5338; 1.1092; 1.8258]; % mycon get stuck for j = 15
% % pars(par_ind) = [17.9780; 7.5624 ; 4.6954; 1.0288; 1.1770; 1.8975]; % ode15i get stuck for j = 15

%% Solve steady state solution with new parameters.

% Initial value
% This initial condition is the steady state data value taken from
% solve_ss_scenario.m.

% Set name for data file to be loaded based upon sex and scenario.    
load_data_name = sprintf('%s_%s_ss_data_scenario_%s.mat', ...
                         species{spe_ind},sex{sex_ind},scenario{1});
% Load data for steady state initial value. 
load(load_data_name, 'SSdata');
SSdataIG = SSdata; clear SSdata;

% Initial guess for the variables.
% Find the steady state solution, so the derivative is 0.
% Arbitrary value for time to input, greater than tchange + deltat.
x0_ss = SSdataIG; x_p0_ss = zeros(num_vars,1); t_ss = 30;

% Time at which to change place holder.
tchange_ss = 0;

exitflag = 0;
iter = 0;

% options_ss = optimset();
% [SSdata, ~, exitflag, ~] = ...
%     fsolve(@(x) ...
%            bp_reg_mod(t_ss,x,x_p0_ss,pars,tchange_ss,varargin_input{:}), ...
%            x0_ss, options_ss);
% 
% if exitflag == 0
%     error('%s fail', sex{sex_ind})
% end

while exitflag == 0 
    fprintf('****** iteration = %s ******', num2str(iter))

    ran_vec = random('unif',0,1,par_num,1);
    ran_vec = lower + ran_vec .* (upper - lower);
    pars(par_ind) = ran_vec;
    
    options_ss = optimset();
    [SSdata, ~, exitflag, ~] = ...
        fsolve(@(x) ...
               bp_reg_mod(t_ss,x,x_p0_ss,pars,tchange_ss,varargin_input{:}), ...
               x0_ss, options_ss);
    
    iter = iter + 1;
end

%% Drugs

% Ang II inf rate fmol/(ml min)
if     strcmp(sex{sex_ind}, 'male')
%     kappa_AngII = 2022; % Sampson 2008
    kappa_AngII = 910; % Sullivan 2010
%     kappa_AngII = 630; % Sullivan 2010
elseif strcmp(sex{sex_ind}, 'female')
%     kappa_AngII = 2060; % Sampson 2008
    kappa_AngII = 505; % Sullivan 2010
%     kappa_AngII = 630; % Sullivan 2010
end

varargin_input = [varargin_input, 'AngII',kappa_AngII];

%% Solve DAE

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
options = odeset('MaxStep',1000); % default is 0.1*abs(t0-tf)

% tic
% Solve dae
[t,x] = ode15i(@(t,x,x_p) ...
               bp_reg_mod(t,x,x_p,pars,tchange,varargin_input{:}), ...
               tspan, x0, x_p0, options);
% toc

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

% Plot Mean Arterial Pressure vs Time. ------------------------------------

% % Data from Sampson 2008. MAP is in difference from baseline.
% Data from Sullivan 2010. MAP is in difference from baseline.
tdata     = [0+1 , 1+1 , 2+1 , 3+1 , 4+1 , 5+1 , 6+1 ,...
             7+1 , 8+1 , 9+1 , 10+1, 11+1, 12+1, 13+1, 14+1];
MAPdata_m = [0   , -1.3, 2.3 , 8.9 , 15.5, 18.3, 22.7, 22.6, ...
             28.6, 31.2, 30.9, 32.8, 37.4, 41.4, 40.3];
MAPdata_f = [0   , 5.2 ,  5.3,  3.9,  3.6,  5.9,    8,   13, ...
             15.7, 17.4, 19.8, 23.7, 25.8,  23.5,  24];

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

end

















