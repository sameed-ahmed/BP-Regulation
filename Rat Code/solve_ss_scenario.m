% This script calculates the steady state values of the blood pressure 
% regulation model bp_reg_solve_baseline.m for some alternative scenarios.
% 
% Steady state data for the intial guess is calculated by solve_ss_baseline.m.

function solve_ss_scenario

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Begin user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Experiment
experiment = 'steady state';

% Scenarios
% AngII  - Ang II infusion
% ACEi   - Angiotensin convernting enzyme inhibitor
% ARB    - Angiotensin receptor blocker
% AT2R-  - Block AT2R through decay
scenario = {'Normal', 'm_RSNA', 'm_AT2R', 'm_RAS', 'm_Reab', ...
            'm_RAS_&_m_Reab', 'm_RSNA_&_m_Reab', ...
            'AngII', 'ACEi', 'ARB', 'AT2R-'};
num_scen = length(scenario);
% Index of scenario to fix.
fixed_ss = 1;

% Species
sp = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of variables.
num_vars = 93;

species = {'human', 'rat'   };
gender  = {'male' , 'female'};

for ss = 1:num_scen % scenario
for gg = 1:2        % gender

%% Parameters

% Parameter input
pars = get_pars(species{sp}, gender{gg}, scenario{ss});

%% Drugs

% drugs = [Ang II inf rate fmol/(ml min), ACEi target level, ARB target level, AT2R decay rate]
if     strcmp(scenario{ss}, 'AngII')
    if     strcmp(gender{gg}, 'male')
        drugs = [2022, 0, 0, 0]; % Sampson 2008
    elseif strcmp(gender{gg}, 'female')
        drugs = [2060, 0, 0, 0]; % Sampson 2008
    end
elseif strcmp(scenario{ss}, 'ACEi')
        drugs = [0, 0.78, 0, 0]; % Leete 2018
elseif strcmp(scenario{ss}, 'ARB')
        drugs = [0, 0, 0.67, 0]; % Leete 2018
elseif strcmp(scenario{ss}, 'AT2R-')
        drugs = [0, 0, 0, 10];
else
        drugs = [0, 0, 0, 0];
end

%% Variables initial guess

% Load data for steady state initial guess. 
% Set name for data file to be loaded based upon gender.    
load_data_name = sprintf('%s_%s_ss_data_scenario_Normal.mat', species{sp}, gender{gg});
% load_data_name = sprintf('NEW%s_ss_data_IG.mat', gender{gg});
load(load_data_name, 'SSdata');
SSdataIG     = SSdata;
SSdata_input = SSdata;
clear SSdata

% load_data_name = sprintf('NEW%s_%s_ss_data_IG.mat', species{sp},gender{gg});
% load(load_data_name, 'SSdataIG');

% Ang II infusion steady state value is sensitive to such a huge change.
% Thus it is done incrementally and iteratively.
if strcmp(scenario{ss}, 'AngII')
    load_data_name = sprintf('%s_%s_ss_data_scenario_AngII.mat', species{sp},gender{gg});
    load(load_data_name, 'SSdata');
    SSdataIG = SSdata;
    clear SSdata;
end

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
% Arbitrary value for time to input.
x0 = SSdataIG; x_p0 = zeros(num_vars,1); t = 0;

% Time at which to change and renal perfusion pressure perturbation 
% place holders.
RPP_per = 0; tchange = 0;

%% Find steady state solution

options = optimset();
[SSdata, residual, ...
 exitflag, output] = fsolve(@(x) bp_reg_mod(t,x,x_p0,pars,SSdata_input ,...
                                            tchange,drugs,RPP_per      ,...
                                            scenario{ss},experiment)   ,...
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

% Save values.
save_data_name = sprintf('%s_%s_ss_data_scenario_%s.mat', ...
                         species{sp}, gender{gg},scenario{ss});
save_data_name = strcat('Data/', save_data_name);
save(save_data_name, 'SSdata', 'residual', 'exitflag', 'output')

end % gender
end % scenario

end






























