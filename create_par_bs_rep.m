function create_par_bs_rep

% clc

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

sample_num = random('Discrete Uniform',1000)

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

% if     strcmp(sex{sex_ind}, 'male'  )
%     AngII_MAP_data = [0   , 1.1 , 2.3 , 8.9 , 15.5, 18.3, 22.7, 22.6, ...
%                28.6, 31.2, 30.9, 32.8, 37.4, 41.4, 40.3];
% elseif strcmp(sex{sex_ind}, 'female')
%     AngII_MAP_data = [0   , 5.2 ,  5.3,  3.9,  3.6,  5.9,    8,   13, ...
%                15.7, 17.4, 19.8, 23.7, 25.8,  23.5,  24];
% end

load_data_name_data = sprintf('%s_%s_AngII_data_bs_rep.mat', species{spe_ind},sex{sex_ind});
load(load_data_name_data, 'AngII_data_rep');

num_sample = size(AngII_data_rep,1);
% num_sample = 10;

%% Find fitted parameters for given data set.

% Initialize parameter replicates corresponding to data replicates.
pars_rep = zeros(152,num_sample);
% Intialize exitflag. Need this because solver gets stuck for certain bad
% initial guesses. So a maximum time is set in the solver options. If
% reached, exitflag = 0. Then this script reruns the solver with a new
% random initial guess.

tic
for j = sample_num:sample_num
% for j = 1:10
% for j = 1:num_sample
% parfor j = 1:10
% [SSdata, pars] = solve_ss_hyp_fit2(sex_ind,AngII_MAP_data);

    iter = 0;
    exitflag_pars = 0;
    while exitflag_pars <= 0
        iter = iter + 1;
%         [pars_rep(:,j),exitflag_pars] = ...
%             solve_ss_hyp_fit2(sex_ind,AngII_data_rep(j,:));
        [pars_rep(:,j),exitflag_pars] = ...
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

% %% Save data.
% 
% save_data_name = sprintf('%s_%s_ss_data_scenario_Pri_Hyp.mat', ...
%                          species{spe_ind},sex{sex_ind});
% save_data_name = strcat('Data/', save_data_name);
% save(save_data_name, 'SSdata')
% save_data_name = sprintf('%s_%s_pars_scenario_Pri_Hyp.mat', ...
%                          species{spe_ind},sex{sex_ind});
% save_data_name = strcat('Data/', save_data_name);
% save(save_data_name, 'pars')
% 
save_data_name = sprintf('%s_%s_pars_scenario_Pri_Hyp_bs_rep.mat', ...
                         species{spe_ind},sex{sex_ind});
save_data_name = strcat('Data/', save_data_name);
save(save_data_name, 'pars_rep', 'num_sample', 'bs_rep_fit_time')

end % sex
delete(gcp)
% delete(myCluster.Jobs)




end































