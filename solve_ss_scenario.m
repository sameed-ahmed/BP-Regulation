% This script calculates the steady state values of the blood pressure 
% regulation model bp_reg_mod.m for some baseline and alternative scenarios.

% Input:  scenario
% Output: saves data files for steady values of model variables

function solve_ss_scenario

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Begin user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scenarios
% Normal - Normal conditions
% m_RSNA - male RSNA
% m_AT2R - male AT2R
% m_RAS  - male RAS pars
% m_Reab - male fractional sodium/water reabsorption

% m_RAS_m_Reab  - male RAS pars & fractional sodium/water reabsorption
% m_RSNA_m_Reab - male RSNA     & fractional sodium/water reabsorption

% AngII - Ang II infusion
% ACEi  - Angiotensin convernting enzyme inhibitor
% ARB   - Angiotensin receptor blocker
scenario = {'Normal', 'm_RSNA', 'm_AT2R', 'm_RAS', 'm_Reab', ...
            'm_RAS_m_Reab', 'm_RSNA_m_Reab', ...
            'AngII', 'ACEi', 'ARB', ...
            'Pri_Hyp'};
num_scen = length(scenario);
fixed_scen = 11;

% Species
spe_ind = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of variables.
num_vars = 93;

species = {'human', 'rat'   };
sex     = {'male' , 'female'};

for sce_ind = fixed_scen:fixed_scen % scenario
for sex_ind = 1:2        % sex

%% Parameters

varargin_input = {scenario{sce_ind},true};

% Parameter input
pars = get_pars(species{spe_ind}, sex{sex_ind}, varargin_input{:});

%% Drugs

% drugs = [Ang II inf rate fmol/(ml min), ACEi target level, ARB target level]
if     strcmp(scenario{sce_ind}, 'AngII')
    if     strcmp(sex{sex_ind}, 'male'  )
        varargin_input = {'AngII',2022}; % Sampson 2008
    elseif strcmp(sex{sex_ind}, 'female')
        varargin_input = {'AngII',2060}; % Sampson 2008
    end
elseif strcmp(scenario{sce_ind}, 'ACEi' )
        varargin_input = {'ACEi' ,0.78 }; % Leete 2018
elseif strcmp(scenario{sce_ind}, 'ARB'  )
        varargin_input = {'ARB'  ,0.67 }; % Leete 2018
end

%% Variables initial guess

% Load data for steady state initial guess. 
% Set name for data file to be loaded based upon sex.

load_data_name = sprintf('%s_%s_ss_data_scenario_%s.mat', ...
                         species{spe_ind},sex{sex_ind},scenario{sce_ind});
load(load_data_name, 'SSdata');
SSdataIG     = SSdata;
clear SSdata

% load_data_name = sprintf('%s_%s_ss_data_IG.mat', ...
%                          species{spe_ind},sex{sex_ind});
% load(load_data_name, 'SSdataIG');

% Order
% x  = [rsna; alpha_map; alpha_rap; R_r; beta_rsna; Phi_rb; Phi_gfilt; ...
%       P_f; P_gh; Sigma_tgf; Phi_filsod; Phi_ptsodreab; eta_ptsodreab; ...
%       gamma_filsod; gamma_at; gamma_rsna; Phi_mdsod; Phi_dtsodreab; ...
%       eta_dtsodreab; psi_al; Phi_dtsod; Phi_cdsodreab; eta_cdsodreab; ...
%       lambda_dt; lambda_anp; lambda_al; Phi_usod; Phi_sodin; V_ecf; ...
%       V_b; P_mf; Phi_vr; Phi_co; P_ra; vas; vas_f; vas_d; R_a; R_ba; ...
%       R_vr; R_tp; P_ma; epsilon_aum; a_auto; a_chemo; a_baro; C_adh; ...
%       N_adh; N_adhs; delta_ra; ...
%       M_sod; C_sod; nu_mdsod; nu_rsna; C_al; N_al; N_als; xi_ksod; ...
%       xi_map; xi_at; hatC_anp; AGT; nu_AT1; R_sec; PRC; PRA; AngI; ...
%       AngII; AT1R; AT2R; Ang17; AngIV; R_aa; R_ea; Sigma_myo; ...
%       Psi_AT1RAA; Psi_AT1REA; Psi_AT2RAA; Psi_AT2REA; ...
%       Phi_ptwreab; eta_ptwreab; mu_ptsodreab; ...
%       Phi_mdu; Phi_dtwreab; eta_dtwreab; mu_dtsodreab; Phi_dtu; ...
%       Phi_cdwreab; eta_cdwreab; mu_cdsodreab; mu_adh; Phi_u; Phi_win];

% Initial guess for the variables.
% Find the steady state solution, so the derivative is 0.
% Arbitrary value for time to input, greater than tchange + deltat.
x0 = SSdataIG; x_p0 = zeros(num_vars,1); t = 30;

% Time at which to change place holder.
tchange = 0;

%% Find steady state solution

options = optimset();
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

%% Save values.

save_data_name = sprintf('%s_%s_ss_data_scenario_%s.mat', ...
                         species{spe_ind},sex{sex_ind},scenario{sce_ind});
save_data_name = strcat('Data/', save_data_name);
save(save_data_name, 'SSdata', 'residual', 'exitflag', 'output')

end % sex
end % scenario

end






























