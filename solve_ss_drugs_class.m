% This simulates the blood pressure regulation model bp_reg.m for drug administration.
% 
% Steady state data is calculated by solve_ss_scenario.m.

function solve_ss_drugs_class

close all

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))
mypath = pwd;
mypath = strcat(mypath, '/Data/Large Files');
addpath(genpath(mypath))

%% Variable names for plotting.
var_names  = {'$rsna$'; '$\alpha_{map}$'; '$\alpha_{rap}$'; '$R_{r}$'; ...
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
              '$\Phi_{u}$'; '$\Phi_{win}$'; '$R_{ea}/R_{r}$'};
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
num_scen1 = length(scenario1);
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
% scenario2 = {'Normal', 'ACEi', 'ARB1', 'CCB', 'DIU', ...
%              'ARB2'  , 'DRI' , 'MRB' , 'RSS', 'AngII'};
scenario2 = {'ACEi', 'ARB1', 'CCB', 'DIU'};
fixed_ss2 = [2];
num_scen2 = length(scenario2);

% Species
spe_ind = 2;

% Bootstrap replicate sample number
% sample_num = random('Discrete Uniform',1000)
% sample_num = 208
% fixed_sample = 1;
num_samples = 1000;

% Drug dose
num_dose = 100;
drug_dose = linspace(0,0.99,num_dose);

% Mean arterial pressure threshold
% MAP_th = 120;
MAP_th = -20;
% Glomerular filtration rate threshold
GFR_th = 25;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pathophysiologically perturbed parameters
pars_ind = [13;14;4;21;18;3;15;41];
pars_hyp_num = length(pars_ind);
%% par names + units
pars_names   = {'$AR$'  , '$VR$'  , ...
                '$AAR$' , '$PRC$' , ...
                '$ALD$' , '$RSNA$', ...
                '$ADH$' , '$MYO$' };
pars_units   = {'$\frac{mmHg}{ml/min}$', '$\frac{mmHg}{ml/min}$', ...
                '$\frac{mmHg}{ml/min}$', '$-$'                  , ...
                '$-$'                  , '$-$'                  , ...
                '$-$'                  , '$-$'                  };
pars_names_des = {'Arterial resist.'           , 'Venous resist.'              , ...
                  'Afferent arteriolar resist.', 'Renin sec. rate'             , ...
                  'Aldosterone sec. rate'      , 'Renal sympathetic nerve act.', ...
                  'Antidiuretic hor. sec. rate', 'Myogenic effect strength'    };
%%

species = {'human', 'rat'   };
sex     = {'male' , 'female'};

% Number of variables
num_vars = 93;

% Initialize matrix to store relatiave change and baseline values for all drugs
MAP_TARGET_M = zeros(         num_samples,num_scen2);
% X_REL_MEAN_M = zeros(num_vars,num_dose   ,num_scen2);
% X_BL_M       = zeros(num_vars,num_samples,num_scen2);
% ---
MAP_TARGET_F = zeros(         num_samples,num_scen2);
% X_REL_MEAN_F = zeros(num_vars,num_dose   ,num_scen2);
% X_BL_F       = zeros(num_vars,num_samples,num_scen2);

%% Load bootstrap replicate parameters.

% Parameters
% load_data_name_pars = sprintf(...
%     '%s_male_pars_scenario_Pri_Hyp_bs_rep1000OLD.mat', species{spe_ind});
load_data_name_pars = sprintf(...
      '%s_male_pars_scenario_Pri_Hyp_bs_rep1000NEWNEW.mat', species{spe_ind});
load(load_data_name_pars, 'pars_rep');
num_pars = size(pars_rep,1);
pars_hyp_m = pars_rep(pars_ind,:);
% ---
% load_data_name_pars = sprintf(...
%     '%s_female_pars_scenario_Pri_Hyp_bs_rep1000OLD.mat', species{spe_ind});
load_data_name_pars = sprintf(...
    '%s_female_pars_scenario_Pri_Hyp_bs_rep1000NEWNEW.mat', species{spe_ind});
load(load_data_name_pars, 'pars_rep');
num_pars = size(pars_rep,1);
pars_hyp_f = pars_rep(pars_ind,:);

% Baseline parameters for relative change
varargin_input = {'Normal',true};
pars_bl_m = get_pars(species{spe_ind}, 'male'  , varargin_input{:}); 
pars_bl_m = pars_bl_m(pars_ind); 
pars_rel_m = pars_hyp_m(:,:) ./ pars_bl_m(:);
% ---
pars_bl_f = get_pars(species{spe_ind}, 'female', varargin_input{:});
pars_bl_f = pars_bl_f(pars_ind); 
pars_rel_f = pars_hyp_f(:,:) ./ pars_bl_f(:);

% Run through each drug.
for scen_ind = 1:num_scen2 % drugs

%% Load bootstrap replicate variables before and after drug dose.

% Variables after drug dose and hypertensiv baseline before drug dose
% X_m/f = (variable, sample, scenario)
% load_data_name_vars = sprintf('%s_male_ss_data_scenario_Pri_Hyp_%s.mat'  , ...
%                               species{spe_ind},scenario2{fixed_ss2});
load_data_name_vars = sprintf('%s_male_ss_data_scenario_Pri_Hyp_%sNEW.mat'  , ...
                              species{spe_ind},scenario2{scen_ind});
load(load_data_name_vars, ...
     'X_ss_m' , 'X_ss_mean_m' , 'X_ss_std_m' , ...
     'X_bl_m' , 'X_bl_mean_m' , 'X_bl_std_m' , ...
     'X_rel_m', 'X_rel_mean_m', 'X_rel_std_m');
% ---
% load_data_name_vars = sprintf('%s_female_ss_data_scenario_Pri_Hyp_%s.mat', ...
%                               species{spe_ind},scenario2{fixed_ss2});
load_data_name_vars = sprintf('%s_female_ss_data_scenario_Pri_Hyp_%sNEW.mat', ...
                              species{spe_ind},scenario2{scen_ind});
load(load_data_name_vars, ...
     'X_ss_f' , 'X_ss_mean_f' , 'X_ss_std_f' , ...
     'X_bl_f' , 'X_bl_mean_f' , 'X_bl_std_f' , ...
     'X_rel_f', 'X_rel_mean_f', 'X_rel_std_f');

% % Retrieve vars for each drug.
% % X_REL_MEAN = zeros(num_vars,num_dose   ,num_scen2);
% % X_BL       = zeros(num_vars,num_samples,num_scen2);
% X_REL_MEAN_M(:,:,scen_ind) = X_rel_mean_m;
% X_BL_M      (:,:,scen_ind) = X_bl_m;
% % ---
% X_REL_MEAN_F(:,:,scen_ind) = X_rel_mean_f;
% X_BL_F      (:,:,scen_ind) = X_bl_f;

%% Target MAP vs Par/Var 

% MAP target --------------------------------------------------------------
% MAP success threshold indices - % change
% X_REL_MEAN = zeros(num_vars,num_dose   ,num_scen2);
% X_BL       = zeros(num_vars,num_samples,num_scen2);

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

%% Set up matrix for decision tree classification.

% Find best drug for each virtual patient. --------------------------------
drug_str = {'ARB', 'CCB', 'TZD'};
MAP_TARGET_M(:,1) = '';
[~,max_ind_m] = max(MAP_TARGET_M,[],2);
best_drug_m = drug_str(max_ind_m)';
% ---
MAP_TARGET_F(:,1) = '';
[~,max_ind_f] = max(MAP_TARGET_F,[],2);
best_drug_f = drug_str(max_ind_f)';

% Create matrix of predictors. --------------------------------------------
% Pars/Vars target
% Interesting pars/vars
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
% feat_m = [vars_reg_m(:,1), pars_reg_m(:,1), vars_reg_m(:,2), ...
%           vars_reg_m(:,3), vars_reg_m(:,4), vars_reg_m(:,5), ...
%           vars_reg_m(:,6), vars_reg_m(:,7)];
data_m = table(vars_reg_m(:,1), pars_reg_m(:,1), vars_reg_m(:,2), ...
               vars_reg_m(:,3), vars_reg_m(:,4), vars_reg_m(:,5), ...
               vars_reg_m(:,6), vars_reg_m(:,7), best_drug_m, ...
               'VariableNames',pars_names);
% ---
pars_reg_f = round(pars_hyp_f(pars_reg_ind,:)',3,'significant');
vars_reg_f = round(X_bl_f    (vars_reg_ind,:)',3,'significant');
% feat_f = [vars_reg_f(:,1), pars_reg_f(:,1), vars_reg_f(:,2), ...
%           vars_reg_f(:,3), vars_reg_f(:,4), vars_reg_f(:,5), ...
%           vars_reg_f(:,6), vars_reg_f(:,7)];
data_f = table(vars_reg_f(:,1), pars_reg_f(:,1), vars_reg_f(:,2), ...
               vars_reg_f(:,3), vars_reg_f(:,4), vars_reg_f(:,5), ...
               vars_reg_f(:,6), vars_reg_f(:,7), best_drug_f, ...
               'VariableNames',pars_names);

%% Save figures and data.

% Save regression -------------------------------------------------------
save_data_name = sprintf('rat_male_ss_data_reg_best_drug_MAP%s.xls', num2str(MAP_th));
save_data_name = strcat('Data/', save_data_name);
writetable(data_m, save_data_name)
% ---
save_data_name = sprintf('rat_female_ss_data_reg_best_drug_MAP%s.xls', num2str(MAP_th));
save_data_name = strcat('Data/', save_data_name);
writetable(data_f, save_data_name)

end






























