% This script calculates the steady state values of the blood pressure 
% regulation model bp_reg_solve_baseline.m by using fsolve.
% It is adopted with modifications from Karaaslan 2005 and Leete 2018.
% 
% Steady state data for the intial guess is inputted by solver_initial_guess_data.m.

function solve_ss_vir_pop

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Begin user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters to perturb
% K_bar            - par 14; , 
% R_bv             - par 15; , 
% C_gcf            - par 6 ; , 
% R_aass           - par 2 ; , 
% R_eass           - par 3 ; , 
% eta_ptsodreab_eq - par 7 ; , 
% eta_dtsodreab_eq - par 8 ; , 
% eta_cdsodreab_eq - par 9 ; , 
% eta_ptwreab_eq   - par 10; , 
% eta_dtwreab_eq   - par 11; , 
% eta_cdwreab_eq   - par 12; , 
% N_rs             - par 22; , 
% X_PRCPRA         - par 23; , 
% k_AGT            - par 32; , 
% c_ACE            - par 33; , 
% c_Chym           - par 34; , 
% c_NEP            - par 35; , 
% c_ACE2           - par 36; , 
% c_IIIV           - par 37; , 
% c_AT1R           - par 38; , 
% c_AT2R           - par 39; , 
% N_als_eq         - par 19; , 
% N_rsna           - par 1 ; , 
% N_adhs_eq        - par 16; , 
% C_K              - par 20; , 
% Indices
par_ind = [14;15;6;2;3;7;8;9;10;11;12;22;23;32;33;34;35;36;37;38;39;19;1;16;20];
par_num = length(par_ind);
% Range for percent change
par_range_lower = [10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10]/100;
par_range_upper = [10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10]/100;

% Variables to check
% P_ma      - var 42; , 
% Phi_co    - var 33; , 
% R_tp      - var 41; , 
% Phi_rb    - var 6 ; , 
% Phi_gfilt - var 7 ; , 
% R_r       - var 4 ; , 
% PRA       - var 79; , 
% C_al      - var 68; , 
% V_ecf     - var 29; , 
% C_adh     - var 47; , 
% C_sod     - var 65; , 
% Indices
var_ind = [42;33;41;6;7;4;79;68;29;47;65];
var_num = length(var_ind);
% Acceptable percent change
var_per_change = 10/100;

% Size of virtual population
N_VP = 10;

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

% Range for acceptable variable value
load(  'male_ss_data_scenario_Normal.mat', 'SSdata');
var_bl_value_m    = SSdata(var_ind);
var_range_lower_m = var_bl_value_m - var_per_change .* var_bl_value_m;
var_range_upper_m = var_bl_value_m + var_per_change .* var_bl_value_m;
clear SSdata
load('female_ss_data_scenario_Normal.mat', 'SSdata');
var_bl_value_f    = SSdata(var_ind);
var_range_lower_f = var_bl_value_f - var_per_change .* var_bl_value_f;
var_range_upper_f = var_bl_value_f + var_per_change .* var_bl_value_f;
clear SSdata

% Number of variables; number of parameters; 
num_vars = 92; num_pars = 45;

gender = {'male', 'female'};

for gg = 1:1:2 % gender

%% VP iteration

% Set VP iteration for each gender.
iter_succ = 0; iter_fail = 0;

% Initialize parameter sets for VP.
pars_VP = zeros(num_pars,N_VP);

if     strcmp(gender{gg}, 'male')
    var_range_lower = var_range_lower_m;
    var_range_upper = var_range_upper_m;
elseif strcmp(gender{gg}, 'female')
    var_range_lower = var_range_lower_f;
    var_range_upper = var_range_upper_f;
end

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
SSdataIG = SSdata;
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

while iter_succ < N_VP % iteration

    %% Uniformly randomly sample perturbed parameters from range.
    
    % Sample between 0 and 1.
    ran_vec = random('unif',0,1,par_num,1);
    % Sample within interval.
    ran_vec = lower + ran_vec .* (upper - lower);
    % Replace input parameters with newly sampled parameter.
    pars(par_ind) = ran_vec;

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
    if       var_range_lower < SSdata(var_ind)
        if   SSdata(var_ind) < var_range_upper
            % Update success iteration.
            iter_succ = iter_succ + 1;
            % Sanity check to see script's progress. 
            fprintf('%s %s iteration = %s out of %s \n', ...
                    scenario{fixed_ss},gender{gg},num2str(iter_succ),num2str(N_VP))
            % Store successful set of parameters.
            pars_VP(:,iter_succ) = pars;
        else
            % Update failure iteration.
            iter_fail = iter_fail + 1;
        end
    else
        % Update failure iteration.
        iter_fail = iter_fail + 1;
    end
    
end % iteration

% Keep track of total number of iterations needed to create VP.
iter_tot = iter_succ + iter_fail;

% Save values.
save_data_name = sprintf('%s_vir_pop_pars_scenario_%s.mat', gender{gg},scenario{fixed_ss});
save_data_name = strcat('Data/', save_data_name);
save(save_data_name, 'pars_VP', 'iter_tot')

end % gender

end






























