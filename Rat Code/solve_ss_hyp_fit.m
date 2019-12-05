% This script calculates the steady state values of the blood pressure 
% regulation model bp_reg_solve_baseline.m by using fsolve.
% It is adopted with modifications from Karaaslan 2005 and Leete 2018.
% 
% Steady state data for the intial guess is inputted by solver_initial_guess_data.m.

function solve_ss_hyp_fit

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Begin user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters to perturb
% K_bar            - par 13; - 0%, +200%
% R_bv             - par 14; - 0%, +200%
% % C_gcf            - par 8 ; -20%
% R_aass           - par 4 ; - 0%, +100%
% N_rs             - par 21; - 0%, +100%
% N_als_eq         - par 18; - 0%, +100%
% N_rsna           - par 3 ; - 0%, +100%
% Indices
par_ind = [13;14;4;21;18;3];
par_num = length(par_ind);
% Range for parameters
par_range_lower = [0  ;0  ;0  ;0  ;0  ;0  ]/100;
par_range_upper = [200;200;100;100;100;100]/100;

% Variables to check
% P_ma      - var 42; +30,40, +40,50
% Phi_co    - var 33; -5%   , +5%
% Phi_rb    - var 6 ; -5%   , +5%
% Phi_gfilt - var 7 ; -5%   , +5%
% Phi_u     - var 92; -5%   , +5%
% Phi_usod  - var 27; -5%   , +5%
% C_sod     - var 52; -2%   , +2%
% Indices
var_ind = [42;33;6;7;92;27;52];
var_num = length(var_ind);
% Range for variables for each sex.
var_range_lower_change_m = [40*100;5;5;5;5;5;2]/100;
var_range_upper_change_m = [50*100;5;5;5;5;5;2]/100;
var_range_lower_change_f = [30*100;5;5;5;5;5;2]/100;
var_range_upper_change_f = [40*100;5;5;5;5;5;2]/100;

% Scenarios
% Normal - Normal conditions
% m_RAS  - male RAS pars
% m_Reab - male fractional sodium and water reabsorption
scenario = {'Normal', 'm_RAS', 'm_Reab', 'm_RAS_m_Reab'};
% Index of scenario to fix.
fixed_ss = 1;

% Species
spe_ind = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of variables; number of parameters; 
num_vars = 93; num_pars = 46; % + SF + fixed_var_pars + SSdata

species = {'human', 'rat'   };
sex     = {'male' , 'female'};

for sex_ind = 1:2 % sex

%% Parameters

varargin_input = {scenario{fixed_ss},true};

% Parameter input
pars0 = get_pars(species{spe_ind}, sex{sex_ind}, varargin_input{:});

% Set interval bounds for parameters.
lower = pars0(par_ind) - par_range_lower .* pars0(par_ind);
upper = pars0(par_ind) + par_range_upper .* pars0(par_ind);

%% Drugs

% drugs = [Ang II inf rate fmol/(ml min), ACEi target level, ARB target level]
if     strcmp(scenario{fixed_ss}, 'AngII')
    if     strcmp(sex{sex_ind}, 'male'  )
        varargin_input = {'AngII',2022}; % Sampson 2008
    elseif strcmp(sex{sex_ind}, 'female')
        varargin_input = {'AngII',2060}; % Sampson 2008
    end
elseif strcmp(scenario{fixed_ss}, 'ACEi' )
        varargin_input = {'ACEi' ,0.78 }; % Leete 2018
elseif strcmp(scenario{fixed_ss}, 'ARB'  )
        varargin_input = {'ARB'  ,0.67 }; % Leete 2018
end

%% Variables initial guess

% Load data for steady state initial guess. 
% Set name for data file to be loaded based upon sex.    
load_data_name = sprintf('%s_%s_ss_data_scenario_%s.mat', ...
                         species{spe_ind},sex{sex_ind},scenario{fixed_ss});
load(load_data_name, 'SSdata');
SSdataIG = SSdata;

% Set acceptable range for certain variables.
% MAP needs special treatment because it is an additive range, 
% not a relative range.
var_bl_value    = SSdata(var_ind(2:end));
if     strcmp(sex{sex_ind},'male')
    var_range_lower_change = var_range_lower_change_m;
    var_range_upper_change = var_range_upper_change_m;
elseif strcmp(sex{sex_ind},'female')
    var_range_lower_change = var_range_lower_change_f;
    var_range_upper_change = var_range_upper_change_f;
end
var_range_lower = var_bl_value - var_range_lower_change(2:end) .* var_bl_value;
var_range_upper = var_bl_value + var_range_upper_change(2:end) .* var_bl_value;
MAP_lower       = SSdata(var_ind(1)) + var_range_lower_change(1);
MAP_upper       = SSdata(var_ind(1)) + var_range_upper_change(1);
var_range_lower = [MAP_lower; var_range_lower];
var_range_upper = [MAP_upper; var_range_upper];

clear SSdata

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

%% Uniformly randomly sample perturbed parameters from range.

% Sample between 0 and 1.
ran_vec = random('unif',0,1,par_num,1);
% Sample within interval.
ran_vec = lower + ran_vec .* (upper - lower);
% Replace input parameters with newly sampled parameter.
pars0(par_ind) = ran_vec;
pars0_est = ran_vec;

%% Optimize.

% pars0
% spe_par = pars0(1);
% sex_par = pars0(2);

% Place holders for fmincon.
A = []; b = []; Aeq = []; beq = []; nonlcon = [];
% Lower and upper bounds for parameters in fmincon.
lb = lower;
ub = upper;
% Edit options for optimizer.
options1 = optimoptions('fmincon', 'Display','iter');
[pars_est_min, residual_pars, exitflag_pars, output_pars] = ...
    fmincon(@(pars_est) ...
            cost_fun(t,x0,x_p0,pars0,pars_est,par_ind,tchange,varargin_input, ...
                     var_ind,var_range_lower,var_range_upper), ...
            pars0_est,A,b,Aeq,beq,lb,ub,nonlcon,options1); %#ok<ASGLU>

% Place estimated pars in proper location.
pars = pars0;
pars(par_ind) = pars_est_min;

% pars_min
% spe_par = pars_min(1);
% sex_par = pars_min(2);

% Solve system with found pars.
options2 = optimset();
[SSdata, residual_ss, exitflag_ss, output_ss  ] = ...
    fsolve(@(x) ...
           bp_reg_mod(t,x,x_p0,pars,tchange,varargin_input{:}), ...
           x0, options2); %#ok<ASGLU>

%% Save values.

% Steady state data
save_data_name = sprintf('%s_%s_ss_data_scenario_Pri_Hyp.mat', ...
                         species{spe_ind},sex{sex_ind});
save_data_name = strcat('Data/', save_data_name);
save(save_data_name, 'SSdata', 'residual_ss', 'exitflag_ss', 'output_ss')
% Parameters
save_data_name = sprintf('%s_%s_pars_scenario_Pri_Hyp.mat', ...
                         species{spe_ind},sex{sex_ind});
save_data_name = strcat('Data/', save_data_name);
save(save_data_name, 'pars', 'residual_pars', 'exitflag_pars', 'output_pars')

end % sex

end

% -------------------------------------------------------------------------
% Cost function
% -------------------------------------------------------------------------

function tot_err = cost_fun(t,x0,x_p0,pars,pars_est,par_ind,tchange,varargin_input, ...
                            var_ind,var_range_lower,var_range_upper)

%% Find steady state solution

% Place estimated pars in proper location.
pars(par_ind) = pars_est;

%     options = optimset();
options = optimset('Display','off');
[SSdata  , ~, ...
 exitflag, ~] = fsolve(@(x) ...
                       bp_reg_mod(t,x,x_p0,pars,tchange,varargin_input{:}), ...
                       x0, options);

% Check for solver convergence.
if exitflag == 0
    tot_err = 1;
    return
end

% Check for imaginary solution.
if not (isreal(SSdata))
    tot_err = 1;
    return
end

% Set any values that are within machine precision of 0 equal to 0.
for i = 1:length(SSdata)
    if abs(SSdata(i)) < eps*100
        SSdata(i) = 0;
    end
end

% Compute error within range of specified variables.
num_vars = length(var_ind);
err = zeros(num_vars,1);
for i = 1:num_vars
    err(i) = max( ( (SSdata(var_ind(i)) - (var_range_lower(i) + var_range_upper(i))/2)^2 - ...
                    (                     (var_range_upper(i) - var_range_lower(i))/2)^2 ) ...
                  / SSdata(var_ind(i))^2, 0 );
end

% Error
tot_err = sqrt(sum(err)) / num_vars;

end






























