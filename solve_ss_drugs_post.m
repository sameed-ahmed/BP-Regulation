% This simulates the blood pressure regulation model bp_reg.m for drug administration.
% 
% Steady state data is calculated by solve_ss_scenario.m.

function solve_ss_drugs_post

close all

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
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
num_scen = length(scenario1);
% Drug scenarios
% Normal - Normal conditions
% ACEi   - Angiotensin converting enzyme inhibitor % 95
% ARB1   - Angiotensin receptor 1 blocker % 94
% CCB    - Calcium channel blocker % 84? 70?
% ARB2   - Angiotensin receptor 2 blocker %
% DRI    - Direct renin inhibitor %
% MRB    - Aldosterone blocker (MR?) %
% RSS    - Renin secretion stimulator (thiazide?) % % NOT COMPLETE
% AngII  - Ang II infusion fmol/(ml min)
scenario2 = {'Normal', 'ACEi', 'ARB1', 'CCB', ...
             'ARB2'  , 'DRI' , 'MRB' , 'RSS', 'AngII'};
fixed_ss2 = [4];

% Species
spe_ind = 2;

% Bootstrap replicate sample number
% sample_num = random('Discrete Uniform',1000)
% sample_num = 208
% sample_num = 655
% fixed_sample = 655;
num_samples = 1000;
% num_samples = 5;
fixed_sample = 1;

% Drug dose
drug_dose = 0.84

% Mean arterial pressure threshold
MAP_th = -25
% Glomerular filtration rate threshold
GFR_th = 25;

% Pathophysiologically perturbed parameters
pars_ind     = [13;14;4;21;18;3];
pars_hyp_num = length(pars_ind);
pars_names   = {'K_{bar}$'            , 'R_{bv}$'             , ...
                'R_{aa-ss}$'          , 'N_{rs}$'             , ...
                'N_{als}^{eq}$'       , 'N_{rsna}$'           };
pars_units   = {'$\frac{mmHg}{ml/min}$', '$\frac{mmHg}{ml/min}$', ...
                '$\frac{mmHg}{ml/min}$', '$-$'                  , ...
                '$-$'                  , '$-$'                  };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

species = {'human', 'rat'   };
sex     = {'male' , 'female'};

% Number of variables
num_vars = 93;

%% Load bootstrap replicate parameters & variables before and after drug dose.

% Parameters
load_data_name_pars = sprintf(  '%s_male_pars_scenario_Pri_Hyp_bs_rep1000.mat', species{spe_ind});
load(load_data_name_pars, 'pars_rep');
num_pars = size(pars_rep,1);
pars_hyp_m = pars_rep(pars_ind,:);
load_data_name_pars = sprintf('%s_female_pars_scenario_Pri_Hyp_bs_rep1000.mat', species{spe_ind});
load(load_data_name_pars, 'pars_rep');
num_pars = size(pars_rep,1);
pars_hyp_f = pars_rep(pars_ind,:);

% Baseline parameters for relative change
varargin_input = {'Normal',true};
pars_bl_m = get_pars(species{spe_ind}, 'male'  , varargin_input{:}); 
pars_bl_f = get_pars(species{spe_ind}, 'female', varargin_input{:});
pars_bl_m = pars_bl_m(pars_ind); pars_bl_f = pars_bl_f(pars_ind); 

% Variables after drug dose and hypertensiv baseline before drug dose
% X_m/f = (variable, sample, scenario)
load_data_name_vars = sprintf('%s_male_ss_data_scenario_Pri_Hyp_%s%s%%.mat'  , ...
                         species{spe_ind},scenario2{fixed_ss2},num2str(drug_dose*100))
load(load_data_name_vars, 'X_ss_m', 'X_bl_m', 'X_rel_m');
load_data_name_vars = sprintf('%s_female_ss_data_scenario_Pri_Hyp_%s%s%%.mat', ...
                         species{spe_ind},scenario2{fixed_ss2},num2str(drug_dose*100))
load(load_data_name_vars, 'X_ss_f', 'X_bl_f', 'X_rel_f');

%% Post processing

% ---

%% Low

% MAP intervals
MAP_ind_low_m = find(142 < X_bl_m(42,:) & X_bl_m(42,:) <= 146);
X_bl_low_m = X_bl_m(:,MAP_ind_low_m); 
X_ss_low_m = X_ss_m(:,MAP_ind_low_m); 
X_rel_low_m = (X_ss_low_m(:,:) - X_bl_low_m(:,:)) ./ X_bl_low_m(:,:) * 100;
% ---
MAP_ind_low_f = find(132 < X_bl_f(42,:) & X_bl_f(42,:) <= 136);
X_bl_low_f = X_bl_f(:,MAP_ind_low_f); 
X_ss_low_f = X_ss_f(:,MAP_ind_low_f); 
X_rel_low_f = (X_ss_low_f(:,:) - X_bl_low_f(:,:)) ./ X_bl_low_f(:,:) * 100;

% MAP success threshold indices
MAP_rel_low_m = X_rel_low_m(42,:);
MAP_success_ind_low_m = find(MAP_rel_low_m - MAP_th <= 0);
MAP_failure_ind_low_m = 1:length(MAP_rel_low_m); MAP_failure_ind_low_m(MAP_success_ind_low_m) = '';
num_success_low_m = length(MAP_success_ind_low_m); num_failure_low_m = length(MAP_failure_ind_low_m); 
% ---
MAP_rel_low_f = X_rel_low_f(42,:);
MAP_success_ind_low_f = find(MAP_rel_low_f - MAP_th <= 0);
MAP_failure_ind_low_f = 1:length(MAP_rel_low_f); MAP_failure_ind_low_f(MAP_success_ind_low_f) = '';
num_success_low_f = length(MAP_success_ind_low_f); num_failure_low_f = length(MAP_failure_ind_low_f); 

% Success/failure hypertensive variables before drug dose
X_success_low_m = X_bl_low_m(:,MAP_success_ind_low_m);
X_failure_low_m = X_bl_low_m(:,MAP_failure_ind_low_m);
% ---
X_success_low_f = X_bl_low_f(:,MAP_success_ind_low_f);
X_failure_low_f = X_bl_low_f(:,MAP_failure_ind_low_f);

% Success/failure hypertensive variables before drug dose mean and standard deviation
X_success_mean_low_m = mean(X_success_low_m,  2);
X_failure_mean_low_m = mean(X_failure_low_m,  2);
X_success_std_low_m  = std (X_success_low_m,0,2);
X_failure_std_low_m  = std (X_failure_low_m,0,2);
% ---
X_success_mean_low_f = mean(X_success_low_f,  2);
X_failure_mean_low_f = mean(X_failure_low_f,  2);
X_success_std_low_f  = std (X_success_low_f,0,2);
X_failure_std_low_f  = std (X_failure_low_f,0,2);

% Relative change between variables before drug dose
X_sf_rel_low_m = (X_success_mean_low_m - X_failure_mean_low_m) ./ X_success_mean_low_m * 100;
% ---
X_sf_rel_low_f = (X_success_mean_low_f - X_failure_mean_low_f) ./ X_success_mean_low_f * 100;

% Round to 3 significant digits.
X_success_mean_low_m = round(X_success_mean_low_m,3,'significant');
X_failure_mean_low_m = round(X_failure_mean_low_m,3,'significant');
X_success_std_low_m  = round(X_success_std_low_m, 3,'significant');
X_failure_std_low_m  = round(X_failure_std_low_m, 3,'significant');
X_sf_rel_low_m       = round(X_sf_rel_low_m,      3,'significant');
X_success_low_m      = round(X_success_low_m,     3,'significant');
X_failure_low_m      = round(X_failure_low_m,     3,'significant');
% ---
X_success_mean_low_f = round(X_success_mean_low_f,3,'significant');
X_failure_mean_low_f = round(X_failure_mean_low_f,3,'significant');
X_success_std_low_f  = round(X_success_std_low_f, 3,'significant');
X_failure_std_low_f  = round(X_failure_std_low_f, 3,'significant');
X_sf_rel_low_f       = round(X_sf_rel_low_f,      3,'significant');
X_success_low_f      = round(X_success_low_f,     3,'significant');
X_failure_low_f      = round(X_failure_low_f,     3,'significant');

% Compute statistical significance of difference between success and
% failure variables before treatment
[h_low_m,p_low_m] = ttest2(X_success_low_m',X_failure_low_m');
sig_low_m = [p_low_m',h_low_m'];
% ---
[h_low_f,p_low_f] = ttest2(X_success_low_f',X_failure_low_f');
sig_low_f = [p_low_f',h_low_f'];

%% Medium

% MAP intervals
MAP_ind_med_m = find(146 < X_bl_m(42,:) & X_bl_m(42,:) <= 150);
X_bl_med_m = X_bl_m(:,MAP_ind_med_m); 
X_ss_med_m = X_ss_m(:,MAP_ind_med_m); 
X_rel_med_m = (X_ss_med_m(:,:) - X_bl_med_m(:,:)) ./ X_bl_med_m(:,:) * 100;
% ---
MAP_ind_med_f = find(136 < X_bl_f(42,:) & X_bl_f(42,:) <= 140);
X_bl_med_f = X_bl_f(:,MAP_ind_med_f); 
X_ss_med_f = X_ss_f(:,MAP_ind_med_f); 
X_rel_med_f = (X_ss_med_f(:,:) - X_bl_med_f(:,:)) ./ X_bl_med_f(:,:) * 100;

% MAP success threshold indices
MAP_rel_med_m = X_rel_med_m(42,:);
MAP_success_ind_med_m = find(MAP_rel_med_m - MAP_th <= 0);
MAP_failure_ind_med_m = 1:length(MAP_rel_med_m); MAP_failure_ind_med_m(MAP_success_ind_med_m) = '';
num_success_med_m = length(MAP_success_ind_med_m); num_failure_med_m = length(MAP_failure_ind_med_m); 
% ---
MAP_rel_med_f = X_rel_med_f(42,:);
MAP_success_ind_med_f = find(MAP_rel_med_f - MAP_th <= 0);
MAP_failure_ind_med_f = 1:length(MAP_rel_med_f); MAP_failure_ind_med_f(MAP_success_ind_med_f) = '';
num_success_med_f = length(MAP_success_ind_med_f); num_failure_med_f = length(MAP_failure_ind_med_f); 

% Success/failure hypertensive variables before drug dose
X_success_med_m = X_bl_med_m(:,MAP_success_ind_med_m);
X_failure_med_m = X_bl_med_m(:,MAP_failure_ind_med_m);
% ---
X_success_med_f = X_bl_med_f(:,MAP_success_ind_med_f);
X_failure_med_f = X_bl_med_f(:,MAP_failure_ind_med_f);

% Success/failure hypertensive variables before drug dose mean and standard deviation
X_success_mean_med_m = mean(X_success_med_m,  2);
X_failure_mean_med_m = mean(X_failure_med_m,  2);
X_success_std_med_m  = std (X_success_med_m,0,2);
X_failure_std_med_m  = std (X_failure_med_m,0,2);
% ---
X_success_mean_med_f = mean(X_success_med_f,  2);
X_failure_mean_med_f = mean(X_failure_med_f,  2);
X_success_std_med_f  = std (X_success_med_f,0,2);
X_failure_std_med_f  = std (X_failure_med_f,0,2);

% Relative change between variables before drug dose
X_sf_rel_med_m = (X_success_mean_med_m - X_failure_mean_med_m) ./ X_success_mean_med_m * 100;
% ---
X_sf_rel_med_f = (X_success_mean_med_f - X_failure_mean_med_f) ./ X_success_mean_med_f * 100;

% Round to 3 significant digits.
X_success_mean_med_m = round(X_success_mean_med_m,3,'significant');
X_failure_mean_med_m = round(X_failure_mean_med_m,3,'significant');
X_success_std_med_m  = round(X_success_std_med_m, 3,'significant');
X_failure_std_med_m  = round(X_failure_std_med_m, 3,'significant');
X_sf_rel_med_m       = round(X_sf_rel_med_m,      3,'significant');
X_success_med_m      = round(X_success_med_m,     3,'significant');
X_failure_med_m      = round(X_failure_med_m,     3,'significant');
% ---
X_success_mean_med_f = round(X_success_mean_med_f,3,'significant');
X_failure_mean_med_f = round(X_failure_mean_med_f,3,'significant');
X_success_std_med_f  = round(X_success_std_med_f, 3,'significant');
X_failure_std_med_f  = round(X_failure_std_med_f, 3,'significant');
X_sf_rel_med_f       = round(X_sf_rel_med_f,      3,'significant');
X_success_med_f      = round(X_success_med_f,     3,'significant');
X_failure_med_f      = round(X_failure_med_f,     3,'significant');

% Compute statistical significance of difference between success and
% failure variables before treatment
[h_med_m,p_med_m] = ttest2(X_success_med_m',X_failure_med_m');
sig_med_m = [p_med_m',h_med_m'];
% ---
[h_med_f,p_med_f] = ttest2(X_success_med_f',X_failure_med_f');
sig_med_f = [p_med_f',h_med_f'];

%% High

% MAP intervals
MAP_ind_hii_m = find(150 < X_bl_m(42,:) & X_bl_m(42,:) <= 154);
X_bl_hii_m = X_bl_m(:,MAP_ind_hii_m); 
X_ss_hii_m = X_ss_m(:,MAP_ind_hii_m); 
X_rel_hii_m = (X_ss_hii_m(:,:) - X_bl_hii_m(:,:)) ./ X_bl_hii_m(:,:) * 100;
% ---
MAP_ind_hii_f = find(140 < X_bl_f(42,:) & X_bl_f(42,:) <= 144);
X_bl_hii_f = X_bl_f(:,MAP_ind_hii_f); 
X_ss_hii_f = X_ss_f(:,MAP_ind_hii_f); 
X_rel_hii_f = (X_ss_hii_f(:,:) - X_bl_hii_f(:,:)) ./ X_bl_hii_f(:,:) * 100;

% MAP success threshold indices
MAP_rel_hii_m = X_rel_hii_m(42,:);
MAP_success_ind_hii_m = find(MAP_rel_hii_m - MAP_th <= 0);
MAP_failure_ind_hii_m = 1:length(MAP_rel_hii_m); MAP_failure_ind_hii_m(MAP_success_ind_hii_m) = '';
num_success_hii_m = length(MAP_success_ind_hii_m); num_failure_hii_m = length(MAP_failure_ind_hii_m); 
% ---
MAP_rel_hii_f = X_rel_hii_f(42,:);
MAP_success_ind_hii_f = find(MAP_rel_hii_f - MAP_th <= 0);
MAP_failure_ind_hii_f = 1:length(MAP_rel_hii_f); MAP_failure_ind_hii_f(MAP_success_ind_hii_f) = '';
num_success_hii_f = length(MAP_success_ind_hii_f); num_failure_hii_f = length(MAP_failure_ind_hii_f); 

% Success/failure hypertensive variables before drug dose
X_success_hii_m = X_bl_hii_m(:,MAP_success_ind_hii_m);
X_failure_hii_m = X_bl_hii_m(:,MAP_failure_ind_hii_m);
% ---
X_success_hii_f = X_bl_hii_f(:,MAP_success_ind_hii_f);
X_failure_hii_f = X_bl_hii_f(:,MAP_failure_ind_hii_f);

% Success/failure hypertensive variables before drug dose mean and standard deviation
X_success_mean_hii_m = mean(X_success_hii_m,  2);
X_failure_mean_hii_m = mean(X_failure_hii_m,  2);
X_success_std_hii_m  = std (X_success_hii_m,0,2);
X_failure_std_hii_m  = std (X_failure_hii_m,0,2);
% ---
X_success_mean_hii_f = mean(X_success_hii_f,  2);
X_failure_mean_hii_f = mean(X_failure_hii_f,  2);
X_success_std_hii_f  = std (X_success_hii_f,0,2);
X_failure_std_hii_f  = std (X_failure_hii_f,0,2);

% Relative change between variables before drug dose
X_sf_rel_hii_m = (X_success_mean_hii_m - X_failure_mean_hii_m) ./ X_success_mean_hii_m * 100;
% ---
X_sf_rel_hii_f = (X_success_mean_hii_f - X_failure_mean_hii_f) ./ X_success_mean_hii_f * 100;

% Round to 3 significant digits.
X_success_mean_hii_m = round(X_success_mean_hii_m,3,'significant');
X_failure_mean_hii_m = round(X_failure_mean_hii_m,3,'significant');
X_success_std_hii_m  = round(X_success_std_hii_m, 3,'significant');
X_failure_std_hii_m  = round(X_failure_std_hii_m, 3,'significant');
X_sf_rel_hii_m       = round(X_sf_rel_hii_m,      3,'significant');
X_success_hii_m      = round(X_success_hii_m,     3,'significant');
X_failure_hii_m      = round(X_failure_hii_m,     3,'significant');
% ---
X_success_mean_hii_f = round(X_success_mean_hii_f,3,'significant');
X_failure_mean_hii_f = round(X_failure_mean_hii_f,3,'significant');
X_success_std_hii_f  = round(X_success_std_hii_f, 3,'significant');
X_failure_std_hii_f  = round(X_failure_std_hii_f, 3,'significant');
X_sf_rel_hii_f       = round(X_sf_rel_hii_f,      3,'significant');
X_success_hii_f      = round(X_success_hii_f,     3,'significant');
X_failure_hii_f      = round(X_failure_hii_f,     3,'significant');

% Compute statistical significance of difference between success and
% failure variables before treatment
[h_hii_m,p_hii_m] = ttest2(X_success_hii_m',X_failure_hii_m');
sig_hii_m = [p_hii_m',h_hii_m'];
% ---
[h_hii_f,p_hii_f] = ttest2(X_success_hii_f',X_failure_hii_f');
sig_hii_f = [p_hii_f',h_hii_f'];

%% Overall

% MAP success threshold indices
MAP_rel_m = X_rel_m(42,:);
MAP_success_ind_m = find(MAP_rel_m - MAP_th <= 0);
MAP_failure_ind_m = 1:num_samples; MAP_failure_ind_m(MAP_success_ind_m) = '';
num_success_m = length(MAP_success_ind_m); num_failure_m = length(MAP_failure_ind_m); 
% ---
MAP_rel_f = X_rel_f(42,:);
MAP_success_ind_f = find(MAP_rel_f - MAP_th <= 0);
MAP_failure_ind_f = 1:num_samples; MAP_failure_ind_f(MAP_success_ind_f) = '';
num_success_f = length(MAP_success_ind_f); num_failure_f = length(MAP_failure_ind_f); 

% Success/failure hypertensive variables before drug dose.
X_success_m = X_bl_m(:,MAP_success_ind_m,fixed_ss1);
X_failure_m = X_bl_m(:,MAP_failure_ind_m,fixed_ss1);
% ---
X_success_f = X_bl_f(:,MAP_success_ind_f,fixed_ss1);
X_failure_f = X_bl_f(:,MAP_failure_ind_f,fixed_ss1);

% Success/failure hypertensive variables before drug dose mean and standard deviation.
X_success_mean_m = mean(X_success_m,  2); 
X_failure_mean_m = mean(X_failure_m,  2); 
X_success_std_m  = std (X_success_m,0,2); 
X_failure_std_m  = std (X_failure_m,0,2); 
% ---
X_success_mean_f = mean(X_success_f,  2);
X_failure_mean_f = mean(X_failure_f,  2);
X_success_std_f  = std (X_success_f,0,2);
X_failure_std_f  = std (X_failure_f,0,2);

% Relative change between variables before drug dose
X_sf_rel_m = (X_success_mean_m - X_failure_mean_m) ./ X_success_mean_m * 100;
% ---
X_sf_rel_f = (X_success_mean_f - X_failure_mean_f) ./ X_success_mean_f * 100;

% Round to 3 significant digits.
X_success_mean_m = round(X_success_mean_m,3,'significant');
X_failure_mean_m = round(X_failure_mean_m,3,'significant');
X_success_std_m  = round(X_success_std_m, 3,'significant');
X_failure_std_m  = round(X_failure_std_m, 3,'significant');
X_sf_rel_m       = round(X_sf_rel_m,      3,'significant');
X_success_m      = round(X_success_m,     3,'significant');
X_failure_m      = round(X_failure_m,     3,'significant');
% ---
X_success_mean_f = round(X_success_mean_f,3,'significant');
X_failure_mean_f = round(X_failure_mean_f,3,'significant');
X_success_std_f  = round(X_success_std_f, 3,'significant');
X_failure_std_f  = round(X_failure_std_f, 3,'significant');
X_sf_rel_f       = round(X_sf_rel_f,      3,'significant');
X_success_f      = round(X_success_f,     3,'significant');
X_failure_f      = round(X_failure_f,     3,'significant');

% Compute statistical significance of difference between success and
% failure variables before treatment.
[h_m,p_m] = ttest2(X_success_m',X_failure_m');
sig_m = [p_m',h_m']; 
% ---
[h_f,p_f] = ttest2(X_success_f',X_failure_f');
sig_f = [p_f',h_f'];

% Compute relative change in parameters.
pars_rel_m = (pars_hyp_m(:,:) - pars_bl_m(:)) ./ pars_bl_m(:) * 100;
% ---
pars_rel_f = (pars_hyp_f(:,:) - pars_bl_f(:)) ./ pars_bl_f(:) * 100;

% Success/failure relative change in parameters from normotensive to hypertensive.
pars_success_m = pars_rel_m(:,MAP_success_ind_m);
pars_failure_m = pars_rel_m(:,MAP_failure_ind_m);
% ---
pars_success_f = pars_rel_f(:,MAP_success_ind_f);
pars_failure_f = pars_rel_f(:,MAP_failure_ind_f);

%% Plot success and failure.

% Plot hypertensive perturbed parameters success and failure. -------------

fp1 = figure('DefaultAxesFontSize',14);
sp1 = gobjects(pars_hyp_num);
for i = 1:pars_hyp_num
    sp1(i) = subplot(3,2,i);
    h1 = histogram(sp1(i),pars_success_m(i,:));
    hold(sp1(i), 'on')
    h2 = histogram(sp1(i),pars_failure_m(i,:));
    hold(sp1(i), 'off')
    
    h1.Normalization = 'probability'; h2.Normalization = 'probability';  
    h1.BinWidth = 10; h2.BinWidth = 10; 
    h1.FaceColor = [0.070, 0.886, 0.333]; h2.FaceColor = [0.917, 0.121, 0.121];

    xlabel_name = strcat('$ \% \Delta ', {' '}, pars_names(i));
    xlabel(sp1(i), xlabel_name, 'Interpreter','latex', 'FontSize',16)    
end
hist_title = sprintf('Male Parameters %s %s%%',scenario2{fixed_ss2},num2str(drug_dose*100));
sgtitle(hist_title, 'FontSize',14)

fp2 = figure('DefaultAxesFontSize',14);
sp2 = gobjects(pars_hyp_num);
for i = 1:pars_hyp_num
    sp2(i) = subplot(3,2,i);
    h1 = histogram(sp2(i),pars_success_f(i,:));
    hold(sp2(i), 'on')
    h2 = histogram(sp2(i),pars_failure_f(i,:));
    hold(sp2(i), 'off')
    
    h1.Normalization = 'probability'; h2.Normalization = 'probability';  
    h1.BinWidth = 10; h2.BinWidth = 10; 
    h1.FaceColor = [0.070, 0.886, 0.333]; h2.FaceColor = [0.917, 0.121, 0.121];

    xlabel_name = strcat('$ \% \Delta ', {' '}, pars_names(i));
    xlabel(sp2(i), xlabel_name, 'Interpreter','latex', 'FontSize',16)    
end
hist_title = sprintf('Female Parameters %s %s%%',scenario2{fixed_ss2},num2str(drug_dose*100));
sgtitle(hist_title, 'FontSize',14)

% Plot all variables success and failure. ---------------------------------

fv1 = gobjects(7,1);
sv1 = gobjects(7,15);
% Loop through each set of subplots.
for i = 1:7
    fv1(i) = figure('pos',[750 500 650 450]);
%     This is to avoid the empty plots in the last subplot set.
    if i == 7
        last_plot = mod(num_vars, 15);
    else
        last_plot = 15;
    end
%     Loop through each subplot within a set of subplots.
    for j = 1:last_plot
        sv1(i,j) = subplot(3,5,j);
%         s(i,j).Position = s(i,j).Position + [0 0 0.01 0];
        
        h1 = histogram(sv1(i,j),X_success_m((i-1)*15 + j,:),3);
        hold on
        h2 = histogram(sv1(i,j),X_failure_m((i-1)*15 + j,:),3);
        hold off
        
        h1.Normalization = 'probability'; h2.Normalization = 'probability'; 
%         h1.BinWidth = 1.0; h2.BinWidth = 1.0; 
        h1.FaceColor = [0.070, 0.886, 0.333]; h2.FaceColor = [0.917, 0.121, 0.121];

        xlabel_name = strcat(var_names((i-1)*15 + j));
        xlabel(sv1(i,j), xlabel_name, 'Interpreter','latex', 'FontSize',16)

%         legend('Male', 'Female')
    end
    hist_title = sprintf('Male Variables %s %s%%',scenario2{fixed_ss2},num2str(drug_dose*100));
    sgtitle(hist_title, 'FontSize',14)
end

fv2 = gobjects(7,1);
sv2 = gobjects(7,15);
% Loop through each set of subplots.
for i = 1:7
    fv2(i) = figure('pos',[750 500 650 450]);
%     This is to avoid the empty plots in the last subplot set.
    if i == 7
        last_plot = mod(num_vars, 15);
    else
        last_plot = 15;
    end
%     Loop through each subplot within a set of subplots.
    for j = 1:last_plot
        sv2(i,j) = subplot(3,5,j);
%         s(i,j).Position = s(i,j).Position + [0 0 0.01 0];
        
        h1 = histogram(sv2(i,j),X_success_f((i-1)*15 + j,:),3);
        hold on
        h2 = histogram(sv2(i,j),X_failure_f((i-1)*15 + j,:),3);
        hold off
        
        h1.Normalization = 'probability'; h2.Normalization = 'probability'; 
%         h1.BinWidth = 1.0; h2.BinWidth = 1.0; 
        h1.FaceColor = [0.070, 0.886, 0.333]; h2.FaceColor = [0.917, 0.121, 0.121];

        xlabel_name = strcat(var_names((i-1)*15 + j));
        xlabel(sv2(i,j), xlabel_name, 'Interpreter','latex', 'FontSize',16)

%         legend('Male', 'Female')
    end
    hist_title = sprintf('Female Variables %s %s%%',scenario2{fixed_ss2},num2str(drug_dose*100));
    sgtitle(hist_title, 'FontSize',14)
end

% Plot mean arterial pressure success and failure. ------------------------

g1 = figure('DefaultAxesFontSize',14);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 4]);
t1 = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');

edges_m = categorical(        {'142-146','146-150','150-154'});
edges_m = reordercats(edges_m,{'142-146','146-150','150-154'});
% ---
edges_f = categorical(        {'132-136','136-140','140-144'});
edges_f = reordercats(edges_f,{'132-136','136-140','140-144'});

succ_bar_m = [num_success_low_m; num_success_med_m; num_success_hii_m];
fail_bar_m = [num_failure_low_m; num_failure_med_m; num_failure_hii_m];
bar_m = [succ_bar_m, fail_bar_m];
% ---
succ_bar_f = [num_success_low_f; num_success_med_f; num_success_hii_f];
fail_bar_f = [num_failure_low_f; num_failure_med_f; num_failure_hii_f];
bar_f = [succ_bar_f, fail_bar_f];

nexttile
b1 = bar(edges_m,bar_m, 1, 'FaceColor','flat');
b1(1).CData = [0, 1, 0];
b1(2).CData = [1, 0, 0];
xlabel('MAP (mmHg)');
ylim([0,500]);
legend('Success','Failure')
title('A')

nexttile
b1 = bar(edges_f,bar_f, 1, 'FaceColor','flat');
b1(1).CData = [0, 1, 0];
b1(2).CData = [1, 0, 0];
xlabel('MAP (mmHg)');
title('B')

%% % Plot some interesting variables
% 
% R_bl_m = reshape(X_bl_m(74,:,:) ./ X_bl_m(4,:,:), [num_samples,num_scen]);
% R_bl_f = reshape(X_bl_f(74,:,:) ./ X_bl_f(4,:,:), [num_samples,num_scen]);
% size(R_bl_m(:,fixed_ss1))
% 
% FRNA_bl_m = reshape((X_bl_m(11,:,:) - X_bl_m(27,:,:)) ./ X_bl_m(11,:,:), [num_samples,num_scen]) * 100;
% FRNA_bl_f = reshape((X_bl_f(11,:,:) - X_bl_f(27,:,:)) ./ X_bl_f(11,:,:), [num_samples,num_scen]) * 100;
% 
% FRW_bl_m = reshape((X_bl_m( 7,:,:) - X_bl_m(92,:,:)) ./ X_bl_m( 7,:,:), [num_samples,num_scen]) * 100;
% FRW_bl_f = reshape((X_bl_f( 7,:,:) - X_bl_f(92,:,:)) ./ X_bl_f( 7,:,:), [num_samples,num_scen]) * 100;
% 
% g1 = figure('DefaultAxesFontSize',14);
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 12, 4]);
% s_mech(1) = subplot(1,3,1); 
% s_mech(2) = subplot(1,3,2); 
% s_mech(3) = subplot(1,3,3); 
% 
% h1 = histogram(s_mech(1),R_bl_m(:,fixed_ss1),10);
% hold(s_mech(1), 'on')
% h2 = histogram(s_mech(1),R_bl_f(:,fixed_ss1),10);
% hold(s_mech(1), 'off')
% h1.Normalization = 'probability'; h2.Normalization = 'probability';  
% h1.BinWidth = 0.01; h2.BinWidth = 0.01; 
% h1.FaceColor = [0.203, 0.592, 0.835]; h2.FaceColor = [0.835, 0.203, 0.576];
% legend(s_mech(1), [h1, h2],{'Male','Female'}, 'FontSize',10,'Location','Northeast');
% xlabel(s_mech(1), 'R_{EA}/R_R');
% title(s_mech(1), 'A')
% 
% h1 = histogram(s_mech(2),FRNA_bl_m(:,fixed_ss1),10);
% hold(s_mech(2), 'on')
% h2 = histogram(s_mech(2),FRNA_bl_f(:,fixed_ss1),10);
% hold(s_mech(2), 'off')
% h1.Normalization = 'probability'; h2.Normalization = 'probability';  
% h1.BinWidth = 0.01; h2.BinWidth = 0.01; 
% h1.FaceColor = [0.203, 0.592, 0.835]; h2.FaceColor = [0.835, 0.203, 0.576];
% xlabel(s_mech(2), 'FR_{Na^+}');
% title(s_mech(2), 'B')
% 
% h1 = histogram(s_mech(3),FRW_bl_m(:,fixed_ss1),10);
% hold(s_mech(3), 'on')
% h2 = histogram(s_mech(3),FRW_bl_f(:,fixed_ss1),10);
% hold(s_mech(3), 'off')
% h1.Normalization = 'probability'; h2.Normalization = 'probability';  
% h1.BinWidth = 0.02; h2.BinWidth = 0.02; 
% h1.FaceColor = [0.203, 0.592, 0.835]; h2.FaceColor = [0.835, 0.203, 0.576];
% xlabel(s_mech(3), 'FR_{U}');
% title(s_mech(3), 'C')

%% % Plot variables that explain mechanims. ---------------------------------
% 
% R_ss_m = reshape(X_ss_m(74,:,:) ./ X_ss_m(4,:,:), [num_samples,num_scen]);
% R_ss_f = reshape(X_ss_f(74,:,:) ./ X_ss_f(4,:,:), [num_samples,num_scen]);
% R_bl_m = reshape(X_bl_m(74,:,:) ./ X_bl_m(4,:,:), [num_samples,num_scen]);
% R_bl_f = reshape(X_bl_f(74,:,:) ./ X_bl_f(4,:,:), [num_samples,num_scen]);
% 
% R_rel_m = (R_ss_m(:,fixed_ss1) - R_bl_m(:,fixed_ss1)) ...
%           ./ R_bl_m(:,fixed_ss1) * 100;
% R_rel_f = (R_ss_f(:,fixed_ss1) - R_bl_f(:,fixed_ss1)) ...
%           ./ R_bl_f(:,fixed_ss1) * 100;
% 
% FRNA_ss_m = reshape((X_ss_m(11,:,:) - X_ss_m(27,:,:)) ./ X_ss_m(11,:,:), [num_samples,num_scen]) * 100;
% FRNA_ss_f = reshape((X_ss_f(11,:,:) - X_ss_f(27,:,:)) ./ X_ss_f(11,:,:), [num_samples,num_scen]) * 100;
% FRNA_bl_m = reshape((X_bl_m(11,:,:) - X_bl_m(27,:,:)) ./ X_bl_m(11,:,:), [num_samples,num_scen]) * 100;
% FRNA_bl_f = reshape((X_bl_f(11,:,:) - X_bl_f(27,:,:)) ./ X_bl_f(11,:,:), [num_samples,num_scen]) * 100;
% 
% FRNA_rel_m = (FRNA_ss_m(:,fixed_ss1) - FRNA_bl_m(:,fixed_ss1)) ...
%              ./ FRNA_bl_m(:,fixed_ss1) * 100;
% FRNA_rel_f = (FRNA_ss_f(:,fixed_ss1) - FRNA_bl_f(:,fixed_ss1)) ...
%              ./ FRNA_bl_f(:,fixed_ss1) * 100;
% 
% mean(FRNA_rel_m)
% mean(FRNA_rel_f)
% 
% FRW_ss_m = reshape((X_ss_m( 7,:,:) - X_ss_m(92,:,:)) ./ X_ss_m( 7,:,:), [num_samples,num_scen]) * 100;
% FRW_ss_f = reshape((X_ss_f( 7,:,:) - X_ss_f(92,:,:)) ./ X_ss_f( 7,:,:), [num_samples,num_scen]) * 100;
% FRW_bl_m = reshape((X_bl_m( 7,:,:) - X_bl_m(92,:,:)) ./ X_bl_m( 7,:,:), [num_samples,num_scen]) * 100;
% FRW_bl_f = reshape((X_bl_f( 7,:,:) - X_bl_f(92,:,:)) ./ X_bl_f( 7,:,:), [num_samples,num_scen]) * 100;
% 
% FRW_rel_m = (FRW_ss_m(:,fixed_ss1) - FRW_bl_m(:,fixed_ss1)) ...
%             ./ FRW_bl_m(:,fixed_ss1) * 100;
% FRW_rel_f = (FRW_ss_f(:,fixed_ss1) - FRW_bl_f(:,fixed_ss1)) ...
%             ./ FRW_bl_f(:,fixed_ss1) * 100;
% 
% mean(FRW_rel_m)
% mean(FRW_rel_f)
% 
% g1 = figure('DefaultAxesFontSize',14);
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 12, 4]);
% s_mech(1) = subplot(1,3,1); 
% s_mech(2) = subplot(1,3,2); 
% s_mech(3) = subplot(1,3,3); 
% 
% h1 = histogram(s_mech(1),R_rel_m(:),10);
% hold(s_mech(1), 'on')
% h2 = histogram(s_mech(1),R_rel_f(:),10);
% hold(s_mech(1), 'off')
% h1.Normalization = 'probability'; h2.Normalization = 'probability';  
% h1.BinWidth = 5.0; h2.BinWidth = 5.0; 
% h1.FaceColor = [0.203, 0.592, 0.835]; h2.FaceColor = [0.835, 0.203, 0.576];
% legend(s_mech(1), [h1, h2],{'Male','Female'}, 'FontSize',10,'Location','Northeast');
% xlabel(s_mech(1), 'R_{EA}/R_R');
% title(s_mech(1), 'A')
% 
% h1 = histogram(s_mech(2),FRNA_rel_m(:),10);
% hold(s_mech(2), 'on')
% h2 = histogram(s_mech(2),FRNA_rel_f(:),10);
% hold(s_mech(2), 'off')
% h1.Normalization = 'probability'; h2.Normalization = 'probability';  
% h1.BinWidth = 0.05; h2.BinWidth = 0.05; 
% h1.FaceColor = [0.203, 0.592, 0.835]; h2.FaceColor = [0.835, 0.203, 0.576];
% xlabel(s_mech(2), 'FR_{Na^+}');
% title(s_mech(2), 'B')
% 
% h1 = histogram(s_mech(3),FRW_rel_m(:),10);
% hold(s_mech(3), 'on')
% h2 = histogram(s_mech(3),FRW_rel_f(:),10);
% hold(s_mech(3), 'off')
% h1.Normalization = 'probability'; h2.Normalization = 'probability';  
% h1.BinWidth = 0.05; h2.BinWidth = 0.05; 
% h1.FaceColor = [0.203, 0.592, 0.835]; h2.FaceColor = [0.835, 0.203, 0.576];
% xlabel(s_mech(3), 'FR_{W}');
% title(s_mech(3), 'C')

%% % Plot Mean Arterial Pressure distribution. ------------------------------
% 
% % Actual, change, and % change in MAP.
% % X_m/f = (variable, sample, scenario)
% MAP_ac_m = zeros(num_samples,num_scen); MAP_ac_f = zeros(num_samples,num_scen);
% MAP_ch_m = zeros(num_samples,num_scen); MAP_ch_f = zeros(num_samples,num_scen);
% MAP_pc_m = zeros(num_samples,num_scen); MAP_pc_f = zeros(num_samples,num_scen);
% for i = 1:num_scen
%     MAP_ac_m(:,i) = (X_ss_m(42,:,i)                 )                        ;
%     MAP_ac_f(:,i) = (X_ss_f(42,:,i)                 )                        ;
%     
%     MAP_ch_m(:,i) = (X_ss_m(42,:,i) - X_bl_m(42,:,i))                        ;
%     MAP_ch_f(:,i) = (X_ss_f(42,:,i) - X_bl_f(42,:,i))                        ;
%     
%     MAP_pc_m(:,i) = (X_ss_m(42,:,i) - X_bl_m(42,:,i)) ./ X_bl_m(42,1,i) * 100;
%     MAP_pc_f(:,i) = (X_ss_f(42,:,i) - X_bl_f(42,:,i)) ./ X_bl_f(42,1,i) * 100;
% end
% 
% g2 = figure('DefaultAxesFontSize',14);
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 12, 4]);
% s_map1(1) = subplot(1,3,1); 
% s_map1(2) = subplot(1,3,2); 
% s_map1(3) = subplot(1,3,3);
% 
% h1 = histogram(s_map1(1),MAP_ac_m(:,fixed_ss1),10);
% hold(s_map1(1), 'on')
% h2 = histogram(s_map1(1),MAP_ac_f(:,fixed_ss1),10);
% hold(s_map1(1), 'off')
% h1.Normalization = 'probability'; h2.Normalization = 'probability';  
% h1.BinWidth = 5.0; h2.BinWidth = 5.0; 
% h1.FaceColor = [0.203, 0.592, 0.835]; h2.FaceColor = [0.835, 0.203, 0.576];
% legend(s_map1(1), [h1, h2],{'Male','Female'}, 'FontSize',10,'Location','Northwest');
% xlabel(s_map1(1), 'MAP (mmHg)');
% title(s_map1(1), 'A')
% 
% h1 = histogram(s_map1(2),MAP_ch_m(:,fixed_ss1),10);
% hold(s_map1(2), 'on')
% h2 = histogram(s_map1(2),MAP_ch_f(:,fixed_ss1),10);
% hold(s_map1(2), 'off')
% h1.Normalization = 'probability'; h2.Normalization = 'probability';  
% h1.BinWidth = 5.0; h2.BinWidth = 5.0; 
% h1.FaceColor = [0.203, 0.592, 0.835]; h2.FaceColor = [0.835, 0.203, 0.576];
% xlabel(s_map1(2), '\DeltaMAP (mmHg)');
% title(s_map1(2), 'B')
% 
% h1 = histogram(s_map1(3),MAP_pc_m(:,fixed_ss1),10);
% hold(s_map1(3), 'on')
% h2 = histogram(s_map1(3),MAP_pc_f(:,fixed_ss1),10);
% hold(s_map1(3), 'off')
% h1.Normalization = 'probability'; h2.Normalization = 'probability';  
% h1.BinWidth = 5.0; h2.BinWidth = 5.0; 
% h1.FaceColor = [0.203, 0.592, 0.835]; h2.FaceColor = [0.835, 0.203, 0.576];
% xlabel(s_map1(3), '% \DeltaMAP');
% title(s_map1(3), 'C')

%% Save figures and data.

% save_data_name = sprintf('success_failure_distribution_%s%s%%.fig', ...
%                          scenario2{fixed_ss2},num2str(drug_dose*100));
% save_data_name = strcat('Figures/', save_data_name);
% savefig([fp1;fp2;fv1;fv2;g1], save_data_name)

% save_data_name = sprintf('%s_male_succfail_scenario_Pri_Hyp_%s%s%%.mat'  , ...
%                          species{spe_ind},scenario2{fixed_ss2},num2str(drug_dose*100));
% save_data_name = strcat('Data/', save_data_name);
% save(save_data_name, 'MAP_th', ...
%                      'sig_m', 'num_success_m', 'num_failure_m', ...
%                      'X_success_mean_m', 'X_success_std_m', ...
%                      'X_failure_mean_m', 'X_failure_std_m', ... 
%                      'X_sf_rel_m', ... % ---
%                      'sig_low_m', 'num_success_low_m', 'num_failure_low_m', ...
%                      'X_success_mean_low_m', 'X_success_std_low_m', ...
%                      'X_failure_mean_low_m', 'X_failure_std_low_m', ... 
%                      'X_sf_rel_low_m', ... % ---
%                      'sig_med_m', 'num_success_med_m', 'num_failure_med_m', ...
%                      'X_success_mean_med_m', 'X_success_std_med_m', ...
%                      'X_failure_mean_med_m', 'X_failure_std_med_m', ... 
%                      'X_sf_rel_med_m', ... % ---
%                      'sig_hii_m', 'num_success_hii_m', 'num_failure_hii_m', ...
%                      'X_success_mean_hii_m', 'X_success_std_hii_m', ...
%                      'X_failure_mean_hii_m', 'X_failure_std_hii_m', ...
%                      'X_sf_rel_hii_m')
% % ---
% save_data_name = sprintf('%s_female_succfail_scenario_Pri_Hyp_%s%s%%.mat', ...
%                          species{spe_ind},scenario2{fixed_ss2},num2str(drug_dose*100));
% save_data_name = strcat('Data/', save_data_name);
% save(save_data_name, 'MAP_th', ...
%                      'sig_f', 'num_success_f', 'num_failure_f', ...
%                      'X_success_mean_f', 'X_success_std_f', ...
%                      'X_failure_mean_f', 'X_failure_std_f', ... 
%                      'X_sf_rel_f', ... % ---
%                      'sig_low_f', 'num_success_low_f', 'num_failure_low_f', ...
%                      'X_success_mean_low_f', 'X_success_std_low_f', ...
%                      'X_failure_mean_low_f', 'X_failure_std_low_f', ... 
%                      'X_sf_rel_low_f', ... % ---
%                      'sig_med_f', 'num_success_med_f', 'num_failure_med_f', ...
%                      'X_success_mean_med_f', 'X_success_std_med_f', ...
%                      'X_failure_mean_med_f', 'X_failure_std_med_f', ... 
%                      'X_sf_rel_med_f', ... % ---
%                      'sig_hii_f', 'num_success_hii_f', 'num_failure_hii_f', ...
%                      'X_success_mean_hii_f', 'X_success_std_hii_f', ...
%                      'X_failure_mean_hii_f', 'X_failure_std_hii_f', ...
%                      'X_sf_rel_hii_f')

end






























