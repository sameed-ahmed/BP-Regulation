% This simulates the blood pressure regulation model bp_reg.m.
% It is adopted with modifications from Karaaslan 2005 and Leete 2018.
% 
% Steady state data is calculated by solve_ss_baseline.m or solve_ss_scenario.m.

% function [SSdata, f] = run_sim
function run_sim

close all

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Begin user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scenarios
% Normal  - Normal conditions
% m_RAS   - male RAS pars
% m_Reab  - male fractional sodium and water reabsorption
% Pri_Hyp - essential/primary hypertension
scenario1 = {'Normal', 'm_RSNA', 'm_AT2R', 'm_RAS', 'm_Reab', ...
             'm_RAS_m_Reab', 'm_RSNA_m_Reab', ...
             'Pri_Hyp'};
scenario2 = {'Normal', 'AngII', 'ACEi', 'ARB'};
fixed_ss1 = 8;
fixed_ss2 = 1;

% Species
spe_ind = 2;

% Number of days to run simulation after change; Day at which to induce change;
days = 1; day_change = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

species = {'human', 'rat'   };
sex     = {'male' , 'female'};

% Number of variables
num_vars = 93;

% Initialize variables and time.
X = cell(1,2);
T = cell(1,2);

for sex_ind = 1:2 % sex

varargin_input = {scenario1{fixed_ss1},true};

%% Parameters

% Parameter input
pars = get_pars(species{spe_ind}, sex{sex_ind}, varargin_input{:});
% pars([1;7;20;29;55]+58)

%% Drugs

% drugs = [Ang II inf rate fmol/(ml min), ACEi target level, ARB target level, AT2R decay rate]
if     strcmp(scenario2{fixed_ss2}, 'AngII')
    if     strcmp(sex{sex_ind}, 'male')
        varargin_input = {varargin_input{:}, 'AngII',2022}; % Sampson 2008
    elseif strcmp(sex{sex_ind}, 'female')
        varargin_input = {varargin_input{:}, 'AngII',2060}; % Sampson 2008
    end
elseif strcmp(scenario2{fixed_ss2}, 'ACEi')
        varargin_input = {varargin_input{:}, 'ACEi',0.78}; % Leete 2018
elseif strcmp(scenario2{fixed_ss2}, 'ARB')
        varargin_input = {varargin_input{:}, 'ARB',0.67}; % Leete 2018
end

%% Solve DAE

% Initial value
% This initial condition is the steady state data value taken from
% solve_ss_scenario.m.

% Set name for data file to be loaded based upon sex and scenario.    
load_data_name = sprintf('%s_%s_ss_data_scenario_%s.mat', ...
                         species{spe_ind},sex{sex_ind},scenario1{fixed_ss1});
% load_data_name = sprintf('%s_%s_ss_data_scenario_%s.mat', ...
%                          species{sp},sex{sex_ind},scenario1{1});
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
tspan = [t0, tend];

% ode options
options = odeset('MaxStep',1); % default is 0.1*abs(t0-tf)

% Solve dae
[t,x] = ode15i(@(t,x,x_p) ...
               bp_reg_mod(t,x,x_p,pars,tchange,varargin_input{:}), ...
               tspan, x0, x_p0, options);

T{sex_ind} = t';
X{sex_ind} = x';

end % sex

%% Plot

% Retrieve male and female.
t_m = T{1}; t_f = T{2};
X_m = X{1}; X_f = X{2};

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
        
        plot(s(i,j), t_m,X_m((i-1)*15 + j,:),'b', t_f,X_f((i-1)*15 + j,:),'r');
        
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

% Save figures.

if strcmp(scenario1{fixed_ss1}, 'Normal') && strcmp(scenario2{fixed_ss2}, 'Normal')
    save_data_name = sprintf('all_vars_baseline.fig');
    save_data_name = strcat('Figures/', save_data_name);
    savefig(f, save_data_name)
end

end






























