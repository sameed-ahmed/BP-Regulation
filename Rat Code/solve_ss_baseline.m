% This script calculates the steady state values of the blood pressure 
% regulation model bp_reg_solve_baseline.m by using fsolve.
% It is adopted with modifications from Karaaslan 2005 and Leete 2018.
% 
% Steady state data for the intial guess is inputted by solver_initial_guess_data.m.

function solve_ss_baseline

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Begin user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scenarios
% Normal - Normal conditions
% m_RAS  - male RAS pars
% m_Reab - male fractional sodium and water reabsorption
scenario = {'Normal', 'm_RSNA', 'm_AT2R', 'm_RAS', 'm_Reab', ...
            'm_RAS_&_m_Reab', 'm_RSNA_&_m_Reab'};
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

for gg = 1:2 % gender

%% Parameters

% Parameter input
pars = get_pars(species{sp}, gender{gg}, scenario{fixed_ss});

%% Drugs

drugs = [0, 0, 0, 0];

%% Variables initial guess

% Load data for steady state initial guess. 
% Set name for data file to be loaded based upon gender.    
load_data_name = sprintf('%s_%s_ss_data_IG.mat', species{sp}, gender{gg});
load(load_data_name, 'SSdataIG');

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

%% Find steady state solution

options = optimset();
[SSdata, residual, ...
 exitflag, output] = fsolve(@(x) bp_reg_solve_baseline(t,x,x_p0,pars,drugs), ...
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

% Save values.
save_data_name = sprintf('NEW%s_%s_ss_data_scenario_%s.mat', ...
                         species{sp}, gender{gg},scenario{fixed_ss});
save_data_name = strcat('Data/', save_data_name);
save(save_data_name, 'SSdata', 'residual', 'exitflag', 'output')

end % gender

end






























