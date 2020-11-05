% This script creates the dataset on which to run the decision tree
% supervised learning classification algorithm that uses the
% pathophysiological profile as predictors for choosing the optimal drug
% class.
% The saved dataset can then be inputted into MATLAB's classification
% learner app.
% 
% Parameters and steady state data for the virtual population are 
% calculated by create_par_bs_rep.m and create_vp.m, respectively.

% Input : threshold percent decrease in MAP at which to compare all drugs
% Output: saves dataset that uses 
%   virtual individuals as instances
%   pathophysiological profile as predictors
%   optimal drug classes as classes

function create_opt_drugs_class_dataset

close all

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))
mypath = pwd;
mypath = strcat(mypath, '/Data/Large Files');
addpath(genpath(mypath))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Begin user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Species
spe_ind = 2;

% Mean arterial pressure threshold
MAP_th = -20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Drug scenarios
% Normal - Normal conditions
% ACEi   - Angiotensin converting enzyme inhibitor % 
% ARB1   - Angiotensin receptor 1 blocker % 
% CCB    - Calcium channel blocker % 
% TZD    - Thiazide diuretic % 
scenario = {'ACEi', 'ARB1', 'CCB', 'TZD'};
num_scen = length(scenario);

species = {'human', 'rat'   };
sex     = {'male' , 'female'};

% Pathophysiologically perturbed parameters
pars_ind = [13;14;4;21;18;3;15;41];
% par names + units
pars_names   = {'$AR$'  , '$VR$'  , ...
                '$AAR$' , '$PRC$' , ...
                '$ALD$' , '$RSNA$', ...
                '$ADH$' , '$MYO$' };

% Drug dose
num_dose = 100;

% Bootstrap replicate sample number
num_samples = 1000;

% Initialize matrix to store relatiave change and baseline values for all drugs
MAP_TARGET_M = zeros(         num_samples,num_scen);
% ---
MAP_TARGET_F = zeros(         num_samples,num_scen);

%% Load bootstrap replicate parameters.

% Parameters
load_data_name_pars = sprintf(...
      '%s_male_pars_scenario_Pri_Hyp_bs_rep1000.mat', species{spe_ind});
load(load_data_name_pars, 'pars_rep');
pars_hyp_m = pars_rep(pars_ind,:);
% ---
load_data_name_pars = sprintf(...
    '%s_female_pars_scenario_Pri_Hyp_bs_rep1000.mat', species{spe_ind});
load(load_data_name_pars, 'pars_rep');
pars_hyp_f = pars_rep(pars_ind,:);

% Run through each drug.
for scen_ind = 1:num_scen % drugs

%% Load bootstrap replicate variables before and after drug dose.

% Variables after drug dose and hypertensive baseline before drug dose
% X_m/f = (variable, sample, scenario)
load_data_name_vars = sprintf('%s_male_ss_data_scenario_Pri_Hyp_%s.mat'  , ...
                              species{spe_ind},scenario{scen_ind});
load(load_data_name_vars, ...
     'X_ss_m' , 'X_ss_mean_m' , 'X_ss_std_m' , ...
     'X_bl_m' , 'X_bl_mean_m' , 'X_bl_std_m' , ...
     'X_rel_m', 'X_rel_mean_m', 'X_rel_std_m');
% ---
load_data_name_vars = sprintf('%s_female_ss_data_scenario_Pri_Hyp_%s.mat', ...
                              species{spe_ind},scenario{scen_ind});
load(load_data_name_vars, ...
     'X_ss_f' , 'X_ss_mean_f' , 'X_ss_std_f' , ...
     'X_bl_f' , 'X_bl_mean_f' , 'X_bl_std_f' , ...
     'X_rel_f', 'X_rel_mean_f', 'X_rel_std_f');

%% Target MAP vs Par/Var 

MAP_rel_mean_m = X_rel_mean_m(42,:);
dose_target_ind_m = find(abs(MAP_rel_mean_m) <= abs(MAP_th), 1,'last');
MAP_rel_m = reshape(X_rel_m(42,:,:), [num_samples,num_dose]);
MAP_TARGET_M(:,scen_ind) = abs(MAP_rel_m(:,dose_target_ind_m));
% ---
MAP_rel_mean_f = X_rel_mean_f(42,:);
dose_target_ind_f = find(abs(MAP_rel_mean_f) <= abs(MAP_th), 1,'last');
MAP_rel_f = reshape(X_rel_f(42,:,:), [num_samples,num_dose]);
MAP_TARGET_F(:,scen_ind) = abs(MAP_rel_f(:,dose_target_ind_f));

end

%% Set up dataset for decision tree classification.

% Find best drug for each virtual patient. --------------------------------
drug_str = {'ARB', 'CCB', 'TZD'};
MAP_TARGET_M(:,1) = '';
[~,max_ind_m] = max(MAP_TARGET_M,[],2);
best_drug_m = drug_str(max_ind_m)';
% ---
MAP_TARGET_F(:,1) = '';
[~,max_ind_f] = max(MAP_TARGET_F,[],2);
best_drug_f = drug_str(max_ind_f)';

% Create dataset. ---------------------------------------------------------
% Pathophysiological pars/vars target
% R_a       var 38
% R_bv      par 02
% R_aa      var 73
% PRC       var 65
% C_al      var 55
% rsna      var 01
% C_adh     var 47
% Sigma_myo var 75
pars_reg_ind = [02];
vars_reg_ind = [38; 73; 65; 55; 01; 47; 75];
pars_names = [pars_names, 'Best Drug'];

pars_reg_m = round(pars_hyp_m(pars_reg_ind,:)',3,'significant');
vars_reg_m = round(X_bl_m    (vars_reg_ind,:)',3,'significant');
data_m = table(vars_reg_m(:,1), pars_reg_m(:,1), vars_reg_m(:,2), ...
               vars_reg_m(:,3), vars_reg_m(:,4), vars_reg_m(:,5), ...
               vars_reg_m(:,6), vars_reg_m(:,7), best_drug_m, ...
               'VariableNames',pars_names);
% ---
pars_reg_f = round(pars_hyp_f(pars_reg_ind,:)',3,'significant');
vars_reg_f = round(X_bl_f    (vars_reg_ind,:)',3,'significant');
data_f = table(vars_reg_f(:,1), pars_reg_f(:,1), vars_reg_f(:,2), ...
               vars_reg_f(:,3), vars_reg_f(:,4), vars_reg_f(:,5), ...
               vars_reg_f(:,6), vars_reg_f(:,7), best_drug_f, ...
               'VariableNames',pars_names);

%% Save data.

save_data_name = sprintf('rat_male_ss_data_reg_best_drug_MAP%s.xls', num2str(MAP_th));
save_data_name = strcat('Data/', save_data_name);
writetable(data_m, save_data_name)
% ---
save_data_name = sprintf('rat_female_ss_data_reg_best_drug_MAP%s.xls', num2str(MAP_th));
save_data_name = strcat('Data/', save_data_name);
writetable(data_f, save_data_name)

end






























