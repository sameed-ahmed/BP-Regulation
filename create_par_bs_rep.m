% This script fits the subset of parameters corresponding to the
% pathophysiology of hypertension.

% The data fitted are three-fold. First, the constraints are satisfied for
% certain key physiological variables to be within range of measurements in
% the SHR. The other two datasets to be fit are perturbation experiments.
% One is dynamic response of MAP to Ang II infusion. The other is steady
% state MAP response to sodium loading.

% This script calls upon solve_ss_hyp_fit.m.

% Warning: This script takes a really long time.

function create_par_bs_rep

clc

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

% Scenarios
% Normal - Normal conditions
scenario = {'Normal'};
% Index of scenario to fix.
fixed_ss = 1;

% Species
spe_ind = 2;

% sample_num = random('Discrete Uniform',1000)
% sample_num = 158

species = {'human', 'rat'   };
sex     = {'male' , 'female'};

% Create parallel pool on cluster. 
% To be used by global optimizer in solve_ss_hyp_fit.m
parpool

% Fit parameters for each sex.
for sex_ind = 1:2 % sex

%% Parameters

% Optional parameter input
varargin_input = {scenario{fixed_ss},true};

% Parameter input
pars0 = get_pars(species{spe_ind}, sex{sex_ind}, varargin_input{:});

%% Variables initial guess

% Load data for steady state initial guess. 
load_data_name_IG = sprintf('%s_%s_ss_data_scenario_%s.mat', ...
                            species{spe_ind},sex{sex_ind},scenario{fixed_ss});
load(load_data_name_IG, 'SSdata');

%% Load bootstrap replicate data created by create_data_bs_rep.m.

load_data_name_data = sprintf('%s_%s_AngII_data_bs_rep.mat', ...
                              species{spe_ind},sex{sex_ind});
load(load_data_name_data, 'AngII_data_rep');

num_sample = size(AngII_data_rep,1);

%% Find fitted parameters for given data set.

% Initialize parameter replicates corresponding to data replicates.
pars_rep = zeros(153,num_sample);
% Initialize residual error.
residual_pars = zeros(1,num_sample);

% Oddly, the solver converges to bad fits, despite using the global
% optimizer. Or it gets stuck at certain internal steps depending upon a
% bad starting value for the parameters. Or it converges to a point for
% which the model is undefined. Therefore, I have nested the solver within 
% a while loop that ensures that the solver reruns with a new starting
% point for the parameters when the time is too long, the model is
% undefined, or the residual error is too large.

% Initialize the iteration count. This is just as a sanity check to see how
% many whole optimizer iterations are needed to find a good set of
% parameters for a given virtual individual. It also provides a saftey to
% see that if the iteration count is too high, then the solver is likely
% stuck.

% Intialize exitflag. Need this because solver gets stuck for certain bad
% initial guesses. So a maximum time is set in the solver options. If
% reached, exitflag = 0. Then this script reruns the solver with a new
% random initial guess.

% Initialize residual error. Keeping track of this because despite using
% the global optimizer, the optimizer may still converge to a bad local
% minimum. Therefore the optimizer is instructed to rerun for a given
% virtual individual if the residual error is too large.

tic
% for j = sample_num:sample_num
for j = 1:num_sample 
    iter = 0;
    exitflag_pars = 0;
    residual_pars(j) = 0.40;
    while exitflag_pars <= 0 || residual_pars(j) >= 0.40
        iter = iter + 1;
        [pars_rep(:,j),residual_pars(j),exitflag_pars] = ...
            solve_ss_hyp_fit(sex_ind,varargin_input,pars0,SSdata,AngII_data_rep(j,:));
        fprintf('********** %s while loop iteration = %s ********** \n', ...
                sex{sex_ind},num2str(iter))
    end
    
    % Sanity check to see script's progress. Also a check for where to
    % troubleshoot in case the solver gets stuck.
    fprintf('********** %s iteration = %s out of %s ********** \n', ...
            sex{sex_ind},num2str(j),num2str(num_sample))
end
bs_rep_fit_time = toc

%% % Find steady state soltuion if single individual. 
% % Comment out if entire population.
% 
% % Number of variables.
% num_vars = 93;
% % Initial guess for the variables.
% % Find the steady state solution, so the derivative is 0.
% % Arbitrary value for time to input, greater than tchange + deltat.
% % Time at which to change place holder.
% x0 = SSdata; x_p0 = zeros(num_vars,1); t = 30; tchange = 0;
% clear SSdata
% 
% options_ss = optimset('Display','off', 'MaxFunEvals',2000);
% [SSdata, ~, ~, ~] = ...
%     fsolve(@(x) ...
%            bp_reg_mod(t,x,x_p0,pars_rep(:,sample_num),tchange,varargin_input{:}), ...
%            x0, options_ss);

%% Save data.

save_data_name = sprintf('%s_%s_pars_scenario_Pri_Hyp_bs_rep1000.mat', ...
                         species{spe_ind},sex{sex_ind});
save_data_name = strcat('Data/', save_data_name);
save(save_data_name, 'pars_rep', 'residual_pars', 'num_sample', 'bs_rep_fit_time')

end % sex

% Shut down parallel pool on cluster. 
delete(gcp)

end % create_par_bs_rep































