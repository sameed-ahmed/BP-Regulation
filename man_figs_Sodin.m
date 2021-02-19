% This simulates the blood pressure regulation model bp_reg_mod.m for 
% sodium loading across the entire virtual population. The result is
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

function man_figs_Sodin

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
% Sodin    - sodium loading
sim_scenario = {'Sodin'};
exact_sim_scen = 1;

% Species
spe_ind = 2;

% Number of variables
num_vars = 93;
% Number of samples
num_sample = 1000;

% Sodium loading initialization -------------------------------------------
% Number of iterations below/above baseline.
iteration = 11; % must be odd number for symmetry
% Fold decrease/increase.
lower = 1/4; upper = 4;
% Arbitrary value for time to input, greater than tchange + deltat.
% Time at which to change place holder.
t = 30; tchange = 0;

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
% MAP = (iteration, samples, sex)
MAP = zeros(2*iteration-1,num_sample,2);

change  = {'decrease', 'increase'};

species = {'human', 'rat'   };
sex     = {'male' , 'female'};

% Create parallel pool on cluster. 
parpool

for sex_ind = 1:2 % sex

varargin_input = {scenario{1},true};

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

% Baseline/range of sodium intake.
Phi_sodin_bl = pars_rep(17,1);
Phi_sodin_range = Phi_sodin_bl * iter_range;

%% Run simulation across whole VP.

parfor sam_ind = 1:num_sample % VP
% for sam_ind = 3:3 % test

% Initialize variables.
% MAP_temp = (iteration)
MAP_temp = zeros(2*iteration-1,1);

% Retreive parameters and steady state data.
pars = pars_rep(:,sam_ind); SSdataIG = SSdata_rep(:,sam_ind);

for cha_ind = 1:2 % change
for iter = 1:iteration % range

%% Parameters

% Vary sodium intake.
if     strcmp(change{cha_ind}, 'decrease')
    Phi_sodin = Phi_sodin_range(iteration-iter+1);
elseif strcmp(change{cha_ind}, 'increase')
    Phi_sodin = Phi_sodin_range(iteration+iter-1);
end
pars(17) = Phi_sodin;

%% Variables initial guess

% Initial guess for the variables.
% Find the steady state solution, so the derivative is 0.
x0 = SSdataIG; x_p0 = zeros(num_vars,1);

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

%% Retrieve Mean Arterial Pressure.

% MAP_temp = (iteration)
if     strcmp(change{cha_ind}, 'decrease')
    MAP_temp(iteration-iter+1) = SSdata(42);
elseif strcmp(change{cha_ind}, 'increase')
    MAP_temp(iteration+iter-1) = SSdata(42);
end

% To avoid the solver not converging, the initial guess for solution to the
% system is taken as the previous solution value. That is, IG_i = SOL_i-1.
% Update next initial guess as current solution.
SSdataIG = SSdata;

end % range
end % change

MAP(:,sam_ind,sex_ind) = MAP_temp;

fprintf('%s simulation = %s out of %s. \n',...
        sex{sex_ind},num2str(sam_ind),num2str(num_sample))

end % VP
end % sex

% Shut down parallel pool on cluster. 
delete(gcp)

%% Post processing

% Retrieve male and female.
% MAP = (iteration, samples, sex)
MAP_m = reshape(MAP(:,:,1), [2*iteration-1,num_sample]); 
MAP_f = reshape(MAP(:,:,2), [2*iteration-1,num_sample]); 

% MAP mean and standard deviation
MAP_mean_m  = mean(MAP_m,2);
% MAP_mean_m  = max(MAP_m,[],2); % test
MAP_std_m   = std (MAP_m ,0,2); 
MAP_lower_m = MAP_mean_m - 2*MAP_std_m; 
MAP_upper_m = MAP_mean_m + 2*MAP_std_m;
% ---
MAP_mean_f  = mean(MAP_f,2);
% MAP_mean_f  = max(MAP_f,[],2); % test
MAP_std_f   = std (MAP_f ,0,2); 
MAP_lower_f = MAP_mean_f - 2*MAP_std_f; 
MAP_upper_f = MAP_mean_f + 2*MAP_std_f;

% Sodium loading data
MAP_bl_m = MAP_mean_m(iteration);
sod_point_m = [1,4];
MAPpoint_mid_m   = [MAP_bl_m,(15+MAP_bl_m + 25+MAP_bl_m)/2];
MAPpoint_range_m = [5  ,(25+  5 - 15+  5)/2];
% ---
MAP_bl_f = MAP_mean_f(iteration);
sod_point_f = [1,4];
MAPpoint_mid_f   = [MAP_bl_f,( 5+MAP_bl_f + 10+MAP_bl_f)/2];
MAPpoint_range_f = [5  ,(10+  5 -  5+  5)/2];

%% Plot

% x-axis
xscale = iter_range;

fig = figure('DefaultAxesFontSize',16*1.5);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7, 5]);
p1 = plot(xscale,MAP_mean_m,'-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',4.5);
ylim([125, 185])
xlim([lower, upper+0.01])
ax = gca;
ax.YTick = (130 : 10 : 180);
ylabel('MAP (mmHg)')
xlabel({'Fold change in sodium excretion'})
hold on
p2 = plot(xscale,MAP_mean_f , '-', 'Color',[0.835, 0.203, 0.576], 'LineWidth',4.5);
p3 = plot(xscale,MAP_lower_m, ':', 'Color',[0.203, 0.592, 0.835], 'LineWidth',1.5);
p4 = plot(xscale,MAP_lower_f, ':', 'Color',[0.835, 0.203, 0.576], 'LineWidth',1.5);
p5 = plot(xscale,MAP_upper_m, ':', 'Color',[0.203, 0.592, 0.835], 'LineWidth',1.5);
p6 = plot(xscale,MAP_upper_f, ':', 'Color',[0.835, 0.203, 0.576], 'LineWidth',1.5);
p7 = errorbar(sod_point_m,MAPpoint_mid_m,MAPpoint_range_m,'o', 'DisplayName',  'Male data', 'Color',[0.203, 0.592, 0.835], 'MarkerSize',9, 'LineWidth',1.5);
p8 = errorbar(sod_point_f,MAPpoint_mid_f,MAPpoint_range_f,'o', 'DisplayName','Female data', 'Color',[0.835, 0.203, 0.576], 'MarkerSize',9, 'LineWidth',1.5);
[~, hobj, ~, ~] = legend([p1 p2 p7 p8],{'Male sim','Female sim','Male data','Female data'}, 'FontSize',8*1.5,'Location','Northwest');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',2.25);
hold off
title('B', 'FontWeight','normal')

%% Save figures.

save_data_name = sprintf('Pri_hyp_sim_%s.png', ...
                         sim_scenario{exact_sim_scen});
save_data_name = strcat('Figures/', save_data_name);
exportgraphics(fig, save_data_name)

end % man_figs
































