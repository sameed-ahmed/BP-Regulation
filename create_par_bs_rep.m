% 

function create_par_bs_rep

clc

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

% Scenarios
% Normal - Normal conditions
% m_RAS  - male RAS pars
% m_Reab - male fractional sodium and water reabsorption
scenario = {'Normal', 'm_RAS', 'm_Reab', 'm_RAS_m_Reab'};
% Index of scenario to fix.
fixed_ss = 1;

% Species
spe_ind = 2;

% sample_num = random('Discrete Uniform',1000)

species = {'human', 'rat'   };
sex     = {'male' , 'female'};

parpool
for sex_ind = 1:2 % sex

%% Parameters

varargin_input = {scenario{fixed_ss},true};

% Parameter input
pars0 = get_pars(species{spe_ind}, sex{sex_ind}, varargin_input{:});

%% Variables initial guess

% Load data for steady state initial guess. 
% Set name for data file to be loaded based upon sex.    
load_data_name_IG = sprintf('%s_%s_ss_data_scenario_%s.mat', ...
                         species{spe_ind},sex{sex_ind},scenario{fixed_ss});
load(load_data_name_IG, 'SSdata');

%% Load bootstrap replicate data created by create_data_bs_rep.m.

load_data_name_data = sprintf('%s_%s_AngII_data_bs_rep.mat', species{spe_ind},sex{sex_ind});
load(load_data_name_data, 'AngII_data_rep');

num_sample = size(AngII_data_rep,1);
% num_sample = 10;

%% Find fitted parameters for given data set.

% Initialize parameter replicates corresponding to data replicates.
pars_rep = zeros(153,num_sample);
% Initialize residual error.
residual_pars = zeros(1,num_sample);
% Intialize exitflag. Need this because solver gets stuck for certain bad
% initial guesses. So a maximum time is set in the solver options. If
% reached, exitflag = 0. Then this script reruns the solver with a new
% random initial guess.

tic
% for j = sample_num:sample_num
% for j = 1:10
for j = 801:1000
% parfor j = 1:10
% [SSdata, pars] = solve_ss_hyp_fit2(sex_ind,AngII_MAP_data);

    iter = 0;
    exitflag_pars = 0;
%     residual_pars = 0.35;
    while exitflag_pars <= 0 %|| residual_pars >= 0.35
        iter = iter + 1;
%         [pars_rep(:,j),exitflag_pars] = ...
%             solve_ss_hyp_fit2(sex_ind,AngII_data_rep(j,:));
        [pars_rep(:,j),residual_pars(j),exitflag_pars] = ...
            solve_ss_hyp_fit2(sex_ind,varargin_input,pars0,SSdata,AngII_data_rep(j,:));
        fprintf('********** %s while loop iteration = %s ********** \n', ...
                sex{sex_ind},num2str(iter))
    end
    
%     [pars_rep(:,j),exitflag_pars] = ...
%             solve_ss_hyp_fit2(sex_ind,AngII_data_rep(j,:));
%     exitflag_pars
        
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

% save_data_name = sprintf('%s_%s_ss_data_scenario_Pri_Hyp.mat', ...
%                          species{spe_ind},sex{sex_ind});
% save_data_name = strcat('Data/', save_data_name);
% save(save_data_name, 'SSdata')
% save_data_name = sprintf('%s_%s_pars_scenario_Pri_Hyp.mat', ...
%                          species{spe_ind},sex{sex_ind});
% save_data_name = strcat('Data/', save_data_name);
% save(save_data_name, 'pars_rep', 'sample_num', 'bs_rep_fit_time')

save_data_name = sprintf('%s_%s_pars_scenario_Pri_Hyp_bs_rep1000.mat', ...
                         species{spe_ind},sex{sex_ind});
save_data_name = strcat('Data/', save_data_name);
save(save_data_name, 'pars_rep', 'residual_pars', 'num_sample', 'bs_rep_fit_time')

end % sex
delete(gcp)
% delete(myCluster.Jobs)

end































