% This script administers each drug for the entire inhibition range across
% the entire virtual population. It then saves the resulting steady state
% values of the model variables after drug administration.

% Warning: This script takes a long time.

% Input
% fixed_ss: index of which drug to administer
% Output
% saves model variable values before drug dose, after drug dose, and 
% relative change, along with mean and standard deviation

function solve_ss_drugs_dose_res

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

% Drug scenarios
% Normal - Normal conditions
% ACEi   - Angiotensin converting enzyme inhibitor % 95
% ARB1   - Angiotensin receptor 1 blocker % 94
% CCB    - Calcium channel blocker % 84
% TZD    - Thiazide diuretic % 0.5 1?
scenario = {'Normal', 'ACEi', 'ARB1', 'CCB', 'TZD'};
fixed_ss = [5];

% Species
spe_ind = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

species = {'human', 'rat'   };
sex     = {'male' , 'female'};

% Number of variables
num_vars = 93;

% Drug dose
% Number of intervals for dose
num_dose = 100;
% Inibition level for primary effect
drug_dose = linspace(0,0.99,num_dose);
% TZD effect on vasodilation
drug_dose_vaso = zeros(1,num_dose);
% TZD effect on renin secretion
a = 11/9; b = 1/9;
% drug_dose_rsec = drug_dose + 0.5 % TZD
% drug_dose_rsec = 2*drug_dose
drug_dose_rsec = a * drug_dose ./ (b + drug_dose);

% Create parallel pool on cluster. 
% To be used by parallel for loop running through bootstrap samples
parpool

for sex_ind = 1:2        % sex

%% Load bootstrap replicate parameters & variables created by create_par_bs_rep.m.

% Parameters
load_data_name_pars = sprintf('%s_%s_pars_scenario_Pri_Hyp_bs_rep1000.mat', ...
                              species{spe_ind},sex{sex_ind});
load(load_data_name_pars, 'pars_rep');
num_pars   = size(pars_rep,1);
% Bootstrap replicate sample number
num_sample = size(pars_rep,2);

% Initialize variables.
% X = (variable, sample, iteration, sex)
X_ss = zeros(num_vars,num_sample,num_dose,2);
X_bl = zeros(num_vars,num_sample,         2);

% Variables
load_data_name_vars = sprintf('%s_%s_ss_data_scenario_Pri_Hyp_bs_rep1000.mat', ...
                              species{spe_ind},sex{sex_ind});
load(load_data_name_vars, 'SSdata_rep');
num_vars = size(SSdata_rep,1);
% Store baseline value to compute relative change.
X_bl(:,:,sex_ind) = SSdata_rep(:,1:num_sample);

parfor sam_iter = 1:num_sample % bootstrap samples
for    dose_iter = 1:num_dose  % dose range

%% Drugs

% Optional parameters.
varargin_input = {'Normal',true};

for i = 1:length(fixed_ss)
    if     strcmp(scenario{fixed_ss(i)}, 'ACEi' )
            varargin_input = [varargin_input, 'ACEi' ,drug_dose(dose_iter)]; % 
    elseif strcmp(scenario{fixed_ss(i)}, 'ARB1' )
            varargin_input = [varargin_input, 'ARB1' ,drug_dose(dose_iter)]; % 
    elseif strcmp(scenario{fixed_ss(i)}, 'CCB'  )
            varargin_input = [varargin_input, 'CCB'  ,[drug_dose(dose_iter),2/3]]; % 
    elseif strcmp(scenario{fixed_ss(i)}, 'TZD'  )
        if     strcmp(sex{sex_ind}, 'male')
            varargin_input = [varargin_input, 'TZD'  ,[drug_dose(dose_iter)/1.0,drug_dose_vaso(dose_iter),drug_dose_rsec(dose_iter)]]; % 
        elseif strcmp(sex{sex_ind}, 'female')
            varargin_input = [varargin_input, 'TZD'  ,[drug_dose(dose_iter)/1.0,drug_dose_vaso(dose_iter),drug_dose_rsec(dose_iter)]]; % 
        end
    end
end

%% Solve system steady state

% Initialization

% Initial guess for the variables.
% Find the steady state solution, so the derivative is 0.
% Arbitrary value for time to input, sufficiently greater than tchange so
% that drug dose reaches steady state value (tanh).
x0_ss = SSdata_rep(:,sam_iter); x_p0_ss = zeros(num_vars,1); t_ss = 2000;

% Time at which to change place holder.
tchange_ss = 0;

% Solver options
options_ss = optimset('Display','off');
% Solve system
[SSdata, residual, ...
 exitflag, output] = fsolve(@(x) ...
                            bp_reg_mod(t_ss,x,x_p0_ss,...
                                       pars_rep(:,sam_iter),tchange_ss,...
                                       varargin_input{:}), ...
                            x0_ss, options_ss);

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

% Store solution.
% X = (variable, sample, iteration, sex)
X_ss(:,sam_iter,dose_iter,sex_ind) = SSdata;

end % range

% Sanity check to see script's progress. Also a check for where to
% troubleshoot in case the solver does not converge.
fprintf('%s sample = %s out of %s \n', ...
        sex{sex_ind},num2str(sam_iter),num2str(num_sample))

end % samples
end % sex

% Shut down parallel pool on cluster. 
delete(gcp)

%% Post processing

% Retrieve male and female data. Find relative change, mean, standard
% deviation, etc. in order to save data for loading in post 
% processing script.

% Retrieve male and female.
% Delete other scenarios for now.
X_ss_m = reshape(X_ss(:,:,:,1), [num_vars,num_sample,num_dose]); 
X_ss_f = reshape(X_ss(:,:,:,2), [num_vars,num_sample,num_dose]); 
X_bl_m = reshape(X_bl(:,:,  1), [num_vars,num_sample         ]); 
X_bl_f = reshape(X_bl(:,:,  2), [num_vars,num_sample         ]); 

% Compute relative change.
X_rel_m = zeros(num_vars,num_sample,num_dose);
X_rel_f = zeros(num_vars,num_sample,num_dose);
for i = 1:num_dose
    X_rel_m(:,:,i) = (X_ss_m(:,:,i) - X_bl_m) ./ X_bl_m * 100;
    X_rel_f(:,:,i) = (X_ss_f(:,:,i) - X_bl_f) ./ X_bl_f * 100;
end

% Compute mean and standard deviation.
X_bl_mean_m  = mean(X_bl_m ,  2); 
X_bl_mean_f  = mean(X_bl_f ,  2);
X_bl_std_m   = std (X_bl_m ,0,2); 
X_bl_std_f   = std (X_bl_f ,0,2);
% ---
X_ss_mean_m  = reshape(mean(X_ss_m ,  2), [num_vars,num_dose]); 
X_ss_mean_f  = reshape(mean(X_ss_f ,  2), [num_vars,num_dose]); 
X_ss_std_m   = reshape(std (X_ss_m ,0,2), [num_vars,num_dose]); 
X_ss_std_f   = reshape(std (X_ss_f ,0,2), [num_vars,num_dose]); 
% ---
X_rel_mean_m = reshape(mean(X_rel_m,  2), [num_vars,num_dose]); 
X_rel_mean_f = reshape(mean(X_rel_f,  2), [num_vars,num_dose]); 
X_rel_std_m  = reshape(std (X_rel_m,0,2), [num_vars,num_dose]); 
X_rel_std_f  = reshape(std (X_rel_f,0,2), [num_vars,num_dose]); 

%% Save figures and data. 

save_data_name = sprintf('%s_male_ss_data_scenario_Pri_Hyp_%s.mat'  , ...
                         species{spe_ind},scenario{fixed_ss});
save_data_name = strcat('Data/Large Files/', save_data_name);
save(save_data_name, 'X_bl_m' , 'X_bl_mean_m' , 'X_bl_std_m' , ...
                     'X_ss_m' , 'X_ss_mean_m' , 'X_ss_std_m' , ...
                     'X_rel_m', 'X_rel_mean_m', 'X_rel_std_m')

save_data_name = sprintf('%s_female_ss_data_scenario_Pri_Hyp_%s.mat', ...
                         species{spe_ind},scenario{fixed_ss});
save_data_name = strcat('Data/Large Files/', save_data_name);
save(save_data_name, 'X_bl_f' , 'X_bl_mean_f' , 'X_bl_std_f' , ...
                     'X_ss_f' , 'X_ss_mean_f' , 'X_ss_std_f' , ...
                     'X_rel_f', 'X_rel_mean_f', 'X_rel_std_f')

end






























