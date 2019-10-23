% This script calculates the steady state values of the blood pressure 
% regulation model bp_reg_solve_baseline.m by using fsolve.
% It is adopted with modifications from Karaaslan 2005 and Leete 2018.
% 
% Steady state data for the intial guess is inputted by solver_initial_guess_data.m.

function solve_ss_hyp

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Begin user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters to perturb
% K_bar            - par 14; - 0%, +200%
% R_bv             - par 15; - 0%, +200%
% % C_gcf            - par 6 ; -20%
% R_aass           - par 2 ; - 0%, +100%
% N_rs             - par 22; - 0%, +100%
% N_als_eq         - par 19; - 0%, +100%
% N_rsna           - par 1 ; - 0%, +100%
% Indices
par_ind = [14;15;2;22;19;1];
par_num = length(par_ind);
% Range for parameters
par_range_lower = [0  ;0  ;0  ;0  ;0  ;0  ]/100;
par_range_upper = [200;200;100;100;100;100]/100;
% Transport parameters to perturb
% eta_ptsodreab_eq - par 7 ; 0.4, 0.9
% eta_dtsodreab_eq - par 8 ; 0.4, 0.9
% eta_cdsodreab_eq - par 9 ; 
% Range for transport parameters
trans_par_range_lower = [0.4; 0.4];
trans_par_range_upper = [0.9; 0.9];

% Variables to check
% P_ma      - var 42; +40, +50
% Phi_co    - var 33; -5%, +5%
% Phi_rb    - var 6 ; -5%, +5%
% Phi_gfilt - var 7 ; -5%, +5%
% Phi_u     - var 63; -5%, +5%
% Phi_usod  - var 27; -5%, +5%
% C_sod     - var 65; -2%, +2%
% Indices
var_ind = [42;33;6;7;63;27;65];
var_num = length(var_ind);
% Range for variables
var_range_lower_change = [40*100;5;5;5;5;5;2]/100;
var_range_upper_change = [50*100;5;5;5;5;5;2]/100;

% Scenarios
% Normal - Normal conditions
% m_RAS  - male RAS pars
% m_Reab - male fractional sodium and water reabsorption
scenario = {'Normal', 'm_RAS', 'm_Reab', 'm_RAS_&_m_Reab'};
% Index of scenario to fix.
fixed_ss = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of variables; number of parameters; 
num_vars = 92; num_pars = 48;

gender = {'male', 'female'};

for gg = 1:2 % gender

%% Parameters

% Parameter input
pars = get_pars(gender{gg}, scenario{fixed_ss});

% Set interval bounds for parameters.
lower = pars(par_ind) - par_range_lower .* pars(par_ind);
upper = pars(par_ind) + par_range_upper .* pars(par_ind);

%% Drugs

% drugs = [Ang II inf rate fmol/(ml min), ACEi target level, ARB target level, AT2R decay rate]
if     strcmp(scenario{fixed_ss}, 'Normal')
    drugs = [0, 0, 0, 0];
elseif strcmp(scenario{fixed_ss}, 'AngII')
    if     strcmp(gender{gg}, 'male')
        drugs = [2022, 0, 0, 0]; % Sampson 2008
    elseif strcmp(gender{gg}, 'female')
        drugs = [2060, 0, 0, 0]; % Sampson 2008
    end
elseif strcmp(scenario{fixed_ss}, 'ACEi')
%     drugs = [0, 1, 0, 0 ]; % Hall 1980
    drugs = [0, 0.78, 0, 0]; % Leete 2018
elseif strcmp(scenario{fixed_ss}, 'ARB')
    drugs = [0, 0, 0.67, 0]; % Leete 2018
elseif strcmp(scenario{fixed_ss}, 'AT2R-')
    drugs = [0, 0, 0, 10];
end

%% Variables initial guess

% Load data for steady state initial guess. 
% Set name for data file to be loaded based upon gender.    
load_data_name = sprintf('%s_ss_data_scenario_Normal.mat', gender{gg});
load(load_data_name, 'SSdata');
% Retrieve and replace parameters in fixed variable equations.
fixed_ind = [2, 10, 14, 24, 44, 49, 66, 71, 88];
fixed_var_pars = SSdata(fixed_ind);
phicophico = SSdata(33); cadhcadh = SSdata(47);
fixed_var_pars = [fixed_var_pars; cadhcadh; phicophico];
SSdata(fixed_ind) = 1;
% Store SS data for initial guess.
SSdataIG = SSdata;

% Set acceptable range for certain variables.
% MAP needs special treatment because it is an additive range, 
% not a relative range.
var_bl_value    = SSdata(var_ind(2:end));
var_range_lower = var_bl_value - var_range_lower_change(2:end) .* var_bl_value;
var_range_upper = var_bl_value + var_range_upper_change(2:end) .* var_bl_value;
MAP_lower       = SSdata(var_ind(1)) + var_range_lower_change(1);
MAP_upper       = SSdata(var_ind(1)) + var_range_upper_change(1);
var_range_lower = [MAP_lower; var_range_lower];
var_range_upper = [MAP_upper; var_range_upper];

% Find filtered sodium load and urine sodium for CD calculation.
Phi_filsod = SSdata(11); Phi_usod = SSdata(27); 

clear SSdata

% Order
% x  = [rsna; alpha_map; alpha_rap; R_r; beta_rsna; Phi_rb; Phi_gfilt; ...
%       P_f; P_gh; Sigma_tgf; Phi_filsod; Phi_ptsodreab; eta_ptsodreab; ...
%       gamma_filsod; gamma_at; gamma_rsna; Phi_mdsod; Phi_dtsodreab; ...
%       eta_dtsodreab; psi_al; Phi_dtsod; Phi_cdsodreab; eta_cdsodreab; ...
%       lambda_dt; lambda_anp; lambda_al; Phi_usod; Phi_win; V_ecf; V_b; ...
%       P_mf; Phi_vr; Phi_co; P_ra; vas; vas_f; vas_d; R_a; R_ba; R_vr; ...
%       R_tp; P_ma; epsilon_aum; a_auto; a_chemo; a_baro; C_adh; N_adh; ...
%       N_adhs; delta_ra; Phi_ptwreab; eta_ptwreab; mu_ptsodreab; ...
%       Phi_mdu; Phi_dtwreab; eta_dtwreab; mu_dtsodreab; Phi_dtu; ...
%       Phi_cdwreab; eta_cdwreab; mu_cdsodreab; mu_adh; Phi_u; M_sod; ...
%       C_sod; nu_mdsod; nu_rsna; C_al; N_al; N_als; xi_ksod; xi_map; ...
%       xi_at; hatC_anp; AGT; nu_AT1; R_sec; PRC; PRA; AngI; AngII; ...
%       AT1R; AT2R; Ang17; AngIV; R_aa; R_ea; Sigma_myo; Psi_AT1RAA; ...
%       Psi_AT1REA; Psi_AT2RAA; Psi_AT2REA];

% Initial guess for the variables.
% Find the steady state solution, so the derivative is 0.
% Arbitrary value for time to input.
x0 = SSdataIG; x_p0 = zeros(num_vars,1); t = 0;

% Failure tracker for each gender; output check boolean;
iter_fail = 0; check_vars = 1;

while check_vars % iteration

    %% Uniformly randomly sample perturbed parameters from range.

    % Sample between 0 and 1.
    ran_vec = random('unif',0,1,par_num,1);
    % Sample within interval.
    ran_vec = lower + ran_vec .* (upper - lower);
    % Replace input parameters with newly sampled parameter.
    pars(par_ind) = ran_vec;
    
    % Sample between 0 and 1.
    trans_ran_vec = random('unif',0,1,2,1);
    % Sample within interval.
    trans_ran_vec = trans_par_range_lower + ...
                    trans_ran_vec .* (trans_par_range_upper - trans_par_range_lower);
    % Calculate CD sodium reabsorption based upon constraint.
    eta_cd = 1 - Phi_usod / (Phi_filsod * (1 - trans_ran_vec(1)) * (1 - trans_ran_vec(2)));
    trans_ran_vec = [trans_ran_vec; eta_cd];
    % Replace input parameters with newly sampled parameter.
    pars([7;8;9]) = trans_ran_vec;

    %% Find steady state solution

%     options = optimset();
    options = optimset('Display','off');
    [SSdata, residual, ...
     exitflag, output] = fsolve(@(x) bp_reg_solve_scenario(t,x,x_p0,pars ,...
                                                           fixed_var_pars,...
                                                           drugs)        ,...
                                x0, options);

    % Check for solver convergence.
    if exitflag == 0
        continue
    end

    % Check for imaginary solution.
    if not (isreal(SSdata))
        continue
    end

    % Set any values that are within machine precision of 0 equal to 0.
    for i = 1:length(SSdata)
        if abs(SSdata(i)) < eps*100
            SSdata(i) = 0;
        end
    end

    % Check if variables of interest are within range.
    if       var_range_lower <= SSdata(var_ind)
        if   SSdata(var_ind) <= var_range_upper
            % Declare success.
            check_vars = 0;
        else
            % Update failure iteration.
            iter_fail = iter_fail + 1;
        end
    else
        % Update failure iteration.
        iter_fail = iter_fail + 1;
    end
    
    % Sanity check to see script progress
    fprintf('%s %s iteration = %s \n', scenario{fixed_ss},gender{gg},num2str(iter_fail))
        
end % iteration

% Save values.
save_data_name = sprintf('%s_ss_hyp_data_&_pars_scenario_%s.mat', gender{gg},scenario{fixed_ss});
save_data_name = strcat('Data/', save_data_name);
save(save_data_name, 'SSdata', 'pars', 'iter_fail')

end % gender

end






























