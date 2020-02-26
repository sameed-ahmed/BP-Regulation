function create_par_bs_rep

clc

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

% Species
spe_ind = 2;

species = {'human', 'rat'   };
sex     = {'male' , 'female'};

parpool
for sex_ind = 1:1 % sex

%% Load bootstrap replicate data created by create_data_bs_rep.m.

% if     strcmp(sex{sex_ind}, 'male'  )
%     AngII_MAP_data = [0   , 1.1 , 2.3 , 8.9 , 15.5, 18.3, 22.7, 22.6, ...
%                28.6, 31.2, 30.9, 32.8, 37.4, 41.4, 40.3];
% elseif strcmp(sex{sex_ind}, 'female')
%     AngII_MAP_data = [0   , 5.2 ,  5.3,  3.9,  3.6,  5.9,    8,   13, ...
%                15.7, 17.4, 19.8, 23.7, 25.8,  23.5,  24];
% end

load_data_name = sprintf('%s_%s_AngII_data_bs_rep.mat', species{spe_ind},sex{sex_ind});
load(load_data_name, 'AngII_data_rep');

num_sample = size(AngII_data_rep,1);
% num_sample = 10;

%% Find fitted parameters for given data set.

pars_rep = zeros(152,num_sample);

tic
for j = 1:10
% [SSdata, pars] = solve_ss_hyp_fit2(sex_ind,AngII_MAP_data);
    pars_rep(:,j) = solve_ss_hyp_fit2(sex_ind,AngII_data_rep(j,:));
    
    % Sanity check to see script's progress. Also a check for where to
    % troubleshoot in case the solver gets stuck.
    fprintf('********** %s iteration = %s out of %s ********** \n', ...
            sex{sex_ind},num2str(j),num2str(num_sample))
end
bs_rep_fit_time = toc

%% Save data.

% save_data_name = sprintf('%s_%s_ss_data_scenario_Pri_Hyp.mat', ...
%                          species{spe_ind},sex{sex_ind});
% save_data_name = strcat('Data/', save_data_name);
% save(save_data_name, 'SSdata')
% save_data_name = sprintf('%s_%s_pars_scenario_Pri_Hyp.mat', ...
%                          species{spe_ind},sex{sex_ind});
% save_data_name = strcat('Data/', save_data_name);
% save(save_data_name, 'pars')

% save_data_name = sprintf('%s_%s_pars_scenario_Pri_Hyp_bs_rep.mat', ...
%                          species{spe_ind},sex{sex_ind});
% save_data_name = strcat('Data/', save_data_name);
% save(save_data_name, 'pars_rep', 'num_sample', 'bs_rep_fit_time')

end % sex
delete(gcp)
% delete(myCluster.Jobs)




end































