% This simulates the blood pressure regulation model bp_reg_mod.m for 
% Ang II infusion across the entire virtual population. The result is
% plotted along with the experimental data used to calibrate the
% pathophysiological parameters related to hypertension.
% 
% Parameters and steady state data for the virtual population are 
% calculated by create_par_bs_rep.m and create_vp.m, respectively.
% 
% This produces Figure 1 in the manuscript.
% 
% Input
% none
% Output
% plots and saves figures corresponding to the scenario.

function man_figs_AngII

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
scenario = {'Normal'};

% Simulation scenarios
% AngII    - Ang II infusion
sim_scenario = {'AngII'};
exact_sim_scen = 1;

% Species
spe_ind = 2;

% Number of variables
num_vars = 93;
% Number of samples
num_sample = 1000;

% Ang II infusion initialization ------------------------------------------
% Number of days to run simulation after change; Day at which to induce change;
days = 14; day_change = 1;
% Number of points for plotting resolution
N = (days+1)*10 + 1;
% Time at which to keep steady state, change a parameter, etc.
tchange = day_change*1440;
% Initial time (min); Final time (min);
t0 = 0*1440; tend = tchange + days*1440;
% Time vector
tspan = linspace(t0,tend,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

species = {'human', 'rat'   };
sex     = {'male' , 'female'};

% Initialize variables.
% MAP = (points, samples, sex)
MAP = zeros(N,num_sample,2);

% Create parallel pool on cluster. 
parpool

for sex_ind = 1:2 % sex

varargin_input = {scenario{1},true};

%% Drugs

% Ang II inf rate fmol/(ml min)
if     strcmp(sex{sex_ind}, 'male')
%     kappa_AngII = 2022; % Sampson 2008
    kappa_AngII = 910; % Sullivan 2010
%     kappa_AngII = 630; % Sullivan 2010
%     kappa_AngII = 0;
elseif strcmp(sex{sex_ind}, 'female')
%     kappa_AngII = 2060; % Sampson 2008
    kappa_AngII = 505; % Sullivan 2010
%     kappa_AngII = 630; % Sullivan 2010
%     kappa_AngII = 0;
end

varargin_input = [varargin_input, 'AngII',kappa_AngII];

%% Load bootstrap replicate parameters & variables created by create_par_bs_rep.m.

% Parameters
load_data_name_pars = sprintf('%s_%s_pars_scenario_Pri_Hyp_bs_rep1000.mat', ...
                              species{spe_ind},sex{sex_ind});
load(load_data_name_pars, 'pars_rep');
num_pars   = size(pars_rep,1);

% Variables
load_data_name_vars = sprintf('%s_%s_ss_data_scenario_Pri_Hyp_bs_rep1000.mat', ...
                              species{spe_ind},sex{sex_ind});
load(load_data_name_vars, 'SSdata_rep');
num_vars   = size(SSdata_rep,1);

%% Run simulation across whole VP.

parfor sam_ind = 1:num_sample % VP
% for sam_ind = 3:3

% Retreive parameters and steady state data.
pars = pars_rep(:,sam_ind); SSdata = SSdata_rep(:,sam_ind);

%% Solve DAE

% Initial condition for the variables and their derivatives. 
% System is initially at steady state, so the derivative is 0.
x0 = SSdata; x_p0 = zeros(num_vars,1);

% ode options
options = odeset('MaxStep',1000); % default MaxStep is 0.1*abs(t0-tf)
% options = odeset();

% Solve dae
[~,x] = ode15i(@(t,x,x_p) ...
               bp_reg_mod(t,x,x_p,pars,tchange,varargin_input{:}), ...
               tspan, x0, x_p0, options);

% This portion is because for some virtual individual parametrizations, the
% DAE solver is crashing due to need a smaller MaxStep. So if the DAE
% solver crashes, choose a smaller step size and resolve it.
test_sim = size(x,1);
% sam_ind
if test_sim < N
    % ode options
    options = odeset('MaxStep',100);

    % Solve dae
    [~,x] = ode15i(@(t,x,x_p) ...
                   bp_reg_mod(t,x,x_p,pars,tchange,varargin_input{:}), ...
                   tspan, x0, x_p0, options);
    test_sim = size(x,1);
    if test_sim < N
        % ode options
        options = odeset('MaxStep',10);

        % Solve dae
        [~,x] = ode15i(@(t,x,x_p) ...
                       bp_reg_mod(t,x,x_p,pars,tchange,varargin_input{:}), ...
                       tspan, x0, x_p0, options);
        test_sim = size(x,1);
        if test_sim < N
            % ode options
            options = odeset('MaxStep',1);

            % Solve dae
            [~,x] = ode15i(@(t,x,x_p) ...
                           bp_reg_mod(t,x,x_p,pars,tchange,varargin_input{:}), ...
                           tspan, x0, x_p0, options);
        end
    end
end

%% Retrieve Mean Arterial Pressure.

% Substract MAP by baseline for each sex and all scenarios.
% MAP = (points, samples, sex)
MAP(:,sam_ind,sex_ind) = x(:,42) - x(1,42);

fprintf('%s simulation = %s out of %s. \n',...
        sex{sex_ind},num2str(sam_ind),num2str(num_sample))

end % VP
end % sex

% Shut down parallel pool on cluster. 
delete(gcp)

%% Post processing

% Retrieve male and female.
MAP_m = reshape(MAP(:,:,1), [N,num_sample]); 
MAP_f = reshape(MAP(:,:,2), [N,num_sample]); 

% MAP mean and standard deviation
MAP_mean_m  = mean(MAP_m,2);
MAP_std_m   = std (MAP_m ,0,2); 
MAP_lower_m = MAP_mean_m - 2*MAP_std_m; 
MAP_upper_m = MAP_mean_m + 2*MAP_std_m;
% ---
MAP_mean_f  = mean(MAP_f,2);
MAP_std_f   = std (MAP_f ,0,2); 
MAP_lower_f = MAP_mean_f - 2*MAP_std_f; 
MAP_upper_f = MAP_mean_f + 2*MAP_std_f;

% Ang II infusion data
load_data_name_AngII = sprintf('%s_male_AngII_data_bs_rep.mat'  , species{spe_ind});
load(load_data_name_AngII, 'AngII_data_rep'); AngII_data_rep = AngII_data_rep';
MAPdata_mean_m = mean(AngII_data_rep,2);
MAPdata_std_m  = std (AngII_data_rep ,0,2);
load_data_name_AngII = sprintf('%s_female_AngII_data_bs_rep.mat', species{spe_ind});
load(load_data_name_AngII, 'AngII_data_rep'); AngII_data_rep = AngII_data_rep';
MAPdata_mean_f = mean(AngII_data_rep,2);
MAPdata_std_f  = std (AngII_data_rep ,0,2);
tdata = [0+1 , 1+1 , 2+1 , 3+1 , 4+1 , 5+1 , 6+1 ,...
         7+1 , 8+1 , 9+1 , 10+1, 11+1, 12+1, 13+1, 14+1];

%% Plot

% x-axis limits
xlower = t0; xupper = tend; 

% Convert minutes to days for longer simulations.
tspan = tspan/1440; tchange = tchange/1440; 
xlower = xlower/1440; xupper = xupper/1440; 

fig = figure('DefaultAxesFontSize',16*1.5);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7, 5]);
p1 = plot(tspan,MAP_mean_m,'-', 'Color',[0.203, 0.592, 0.835], 'LineWidth',4.5);
xlim([xlower, xupper]); ylim([-1, 60]);
ax = gca;
ax.XTick = (tchange+0*(1) : 2 : tchange+days*(1));
ax.XTickLabel = {'0','2','4','6','8','10','12','14'};
xlabel('Time (days)'); ylabel('\DeltaMAP (mmHg)');
hold on
p2 = plot(tspan,MAP_mean_f , '-', 'Color',[0.835, 0.203, 0.576], 'LineWidth',4.5);
p3 = plot(tspan,MAP_lower_m, ':', 'Color',[0.203, 0.592, 0.835], 'LineWidth',1.5);
p4 = plot(tspan,MAP_lower_f, ':', 'Color',[0.835, 0.203, 0.576], 'LineWidth',1.5);
p5 = plot(tspan,MAP_upper_m, ':', 'Color',[0.203, 0.592, 0.835], 'LineWidth',1.5);
p6 = plot(tspan,MAP_upper_f, ':', 'Color',[0.835, 0.203, 0.576], 'LineWidth',1.5);
p7 = errorbar(tdata,MAPdata_mean_m,2*MAPdata_std_m,'o', 'DisplayName',  'Male data', 'Color',[0.203, 0.592, 0.835], 'MarkerSize',9, 'LineWidth',1.5);
p8 = errorbar(tdata,MAPdata_mean_f,2*MAPdata_std_f,'o', 'DisplayName','Female data', 'Color',[0.835, 0.203, 0.576], 'MarkerSize',9, 'LineWidth',1.5);
[~, hobj, ~, ~] = legend([p1 p2 p7 p8],{'Male sim','Female sim','Male data','Female data'}, 'FontSize',8*1.5,'Location','Northwest');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',2.25);
title('A', 'FontWeight','normal')

%% Save figures.

save_data_name = sprintf('Pri_hyp_sim_%s.png', ...
                         sim_scenario{exact_sim_scen});
save_data_name = strcat('Figures/', save_data_name);
exportgraphics(fig, save_data_name)

end % man_figs
































