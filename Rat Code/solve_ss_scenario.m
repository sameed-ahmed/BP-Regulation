% This script calculates the steady state solution to the system using
% fsolve. Some previous values are used as an initial guess. These are
% taken from Jessica, which are taken in part from some paper (Karaaslan
% 2005?).

% function SSdata = solve_ss
function solve_ss_scenario

% Scenarios
% AngII - Ang II infusion
% ACEi  - Angiotensin convernting enzyme inhibitor
% ARB   - Angiotensin receptor blocker
% AT2R- - Block AT2R through decay
% RHyp  - Renal hypertension due to increased afferent arteriolar resistance
scenario = {'AngII', 'ACEi', 'ARB', 'AT2R-', 'RHyp'};
ss = 1;

num_vars = 92;

gender = {'male', 'female'};

X        = zeros(num_vars,2);
RESIDUAL = zeros(num_vars,2);
EXITFLAG = zeros(1 ,2);
OUTPUT   = cell (1 ,2);

for gg = 1:2 % gender

%% Parameters

if     strcmp(gender{gg}, 'male')
    gen = 1;
elseif strcmp(gender{gg}, 'female')
    gen = 0;
end

% Scaling factor
% Rat flow = Human flow x SF
if     strcmp(gender{gg}, 'male')
    SF = 4.5*10^(-3)*10^(3);
elseif strcmp(gender{gg}, 'female')
    SF = 4.5*10^(-3)*10^(3);
end
% Rat resistance = Human resistance x SF
% Note: This includes conversion from l to ml.
if     strcmp(gender{gg}, 'male')
    SF_R = 0.343;
elseif strcmp(gender{gg}, 'female')
    SF_R = 0.537;
end
% Rat volume = Human volume x SF
% Note: This includes conversion from l to ml.
if     strcmp(gender{gg}, 'male')
    SF_V = 3;
elseif strcmp(gender{gg}, 'female')
    SF_V = 2.4;
end

N_rsna    = 1.00;
% R_aass    = 31.67 / SF;   % mmHg min / ml
% R_eass    = 51.66 / SF;   % mmHg min / ml
if     strcmp(gender{gg}, 'male')
R_aass    = 10.87;   % mmHg min / ml
R_eass    = 17.74;   % mmHg min / ml
elseif strcmp(gender{gg}, 'female')
R_aass    = 17.02;   % mmHg min / ml
R_eass    = 27.76;   % mmHg min / ml
end
if     strcmp(scenario{ss}, 'RHyp') %|| strcmp(scenario{ss}, 'ACEi') || strcmp(scenario{ss}, 'ARB')
    R_aass = R_aass*3.5;
end
P_B       = 18;           % mmHg
P_go      = 28;           % mmHg
% C_gcf     = 0.00781 * SF;
if     strcmp(gender{gg}, 'male')
C_gcf     = 0.068;
elseif strcmp(gender{gg}, 'female')
C_gcf     = 0.047;
end
if     strcmp(gender{gg}, 'male')
    eta_ptsodreab_eq = 0.93; 
    eta_dtsodreab_eq = 0.77; 
    eta_cdsodreab_eq = 0.15;
elseif strcmp(gender{gg}, 'female')
    eta_ptsodreab_eq = 0.90;
    eta_dtsodreab_eq = 0.77; 
    eta_cdsodreab_eq = 0.15;
end
if     strcmp(gender{gg}, 'male')
    eta_ptwreab_eq = 0.86; 
    eta_dtwreab_eq = 0.60; 
    eta_cdwreab_eq = 0.78;
elseif strcmp(gender{gg}, 'female')
    eta_ptwreab_eq = 0.80;
    eta_dtwreab_eq = 0.60; 
    eta_cdwreab_eq = 0.78;
end
% K_vd      = 0.00001;
K_vd      = 0.01;
% K_bar     = 16.6 / SF;    % mmHg min / ml
K_bar     = 16.6 * SF_R;    % mmHg min / ml
% R_bv      = 3.4 / SF;     % mmHg min / ml
R_bv      = 3.4 * SF_R;     % mmHg min / ml
T_adh     = 6;            % min
% Phi_sodin = 1.2278;       % microEq / min
Phi_sodin = 2.3875;       % microEq / min
C_K       = 5;            % microEq / ml 
T_al      = 30;           % min LISTED AS 30 IN TABLE %listed as 60 in text will only change dN_al
N_rs      = 1;            % ng / ml / min

% RAS
h_renin   = 12;      % min
h_AGT     = 10*60;   % min
h_AngI    = 0.5;     % min
h_AngII   = 0.66;    % min
h_Ang17   = 30;      % min
h_AngIV   = 0.5;     % min
h_AT1R    = 12;      % min
h_AT2R    = 12;      % min

% Male and female different parameters for RAS
if     strcmp(gender{gg}, 'male')
    X_PRCPRA = 135.59/17.312;
    k_AGT    = 801.02;
    c_ACE    = 0.096833;
    c_Chym   = 0.010833;
    c_NEP    = 0.012667;
    c_ACE2   = 0.0026667;
    c_IIIV   = 0.29800;
    c_AT1R   = 0.19700;
    c_AT2R   = 0.065667;
    AT1R_eq  = 20.4807902818665;
    AT2R_eq  = 6.82696474842298;
elseif strcmp(gender{gg}, 'female')
    X_PRCPRA = 114.22/17.312;
    k_AGT    = 779.63;
    c_ACE    = 0.11600;
    c_Chym   = 0.012833;
    c_NEP    = 0.0076667;
    c_ACE2   = 0.00043333;
    c_IIIV   = 0.29800;
    c_AT1R   = 0.19700;
    c_AT2R   = 0.065667;
    AT1R_eq  = 20.4538920068419;
    AT2R_eq  = 6.81799861123497;
end

pars = [N_rsna; R_aass; R_eass; P_B; P_go; C_gcf; eta_ptsodreab_eq; ...
        eta_dtsodreab_eq; eta_cdsodreab_eq; eta_ptwreab_eq; ...
        eta_dtwreab_eq; eta_cdwreab_eq; K_vd; K_bar; R_bv; T_adh; ...
        Phi_sodin; C_K; T_al; N_rs; X_PRCPRA; h_renin; h_AGT; h_AngI; ...
        h_AngII; h_Ang17; h_AngIV; h_AT1R; h_AT2R; k_AGT; c_ACE; ...
        c_Chym; c_NEP; c_ACE2; c_IIIV; c_AT1R; c_AT2R; AT1R_eq; ...
        AT2R_eq; gen; SF; SF_R; SF_V];

%% Drugs

% drugs = [Ang II inf rate fmol/(ml min), ACEi target level, ARB target level, AT2R decay rate]
if     strcmp(scenario{ss}, 'Normal') || strcmp(scenario{ss}, 'RHyp')
    drugs = [0, 0, 0, 0];
elseif strcmp(scenario{ss}, 'AngII')
    if     strcmp(gender{gg}, 'male')
        drugs = [2022, 0, 0, 0]; % Sampson 2008
    elseif strcmp(gender{gg}, 'female')
        drugs = [2060, 0, 0, 0]; % Sampson 2008
    end
elseif strcmp(scenario{ss}, 'ACEi')
%     drugs = [0, 1, 0, 0 ]; % Hall 1980
    drugs = [0, 0.78, 0, 0]; % Leete 2018
elseif strcmp(scenario{ss}, 'ARB')
    drugs = [0, 0, 0.67, 0]; % Leete 2018
elseif strcmp(scenario{ss}, 'AT2R-')
    drugs = [0, 0, 0, 10];
end

%% Variables initial guess

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

% Load data for steady state initial value. 
if     strcmp(gender{gg}, 'male')
    load(  'male_ss_data_scenario_Normal.mat', 'SSdata');
elseif strcmp(gender{gg}, 'female')
    load('female_ss_data_scenario_Normal.mat', 'SSdata');
end

% Retrieve and replace parameters in fixed variable equations.
fixed_ind = [2, 10, 14, 24, 44, 49, 66, 71, 88];
fixed_var_pars = SSdata(fixed_ind);
phicophico = SSdata(33);
cadhcadh   = SSdata(47);
fixed_var_pars = [fixed_var_pars; cadhcadh; phicophico];
SSdata(fixed_ind) = 1;
SSdataIG = SSdata;
clear SSdata

% Ang II infusion steady state value is sensitive to such a huge change.
% Thus it is done incrementally and iteratively.
if     strcmp(scenario{ss}, 'AngII')
    if     strcmp(gender{gg}, 'male')
        load(  'male_ss_data_scenario_AngII.mat', 'SSdata');
        SSdataIG = SSdata;
        clear SSdata;
    elseif strcmp(gender{gg}, 'female')
        load('female_ss_data_scenario_AngII.mat', 'SSdata');
        SSdataIG = SSdata;
        clear SSdata;
    end
end

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

%% Find steady state solution

options = optimset(); %options = optimset('MaxFunEvals',num_vars*100+10000);
[SSdata, residual, ...
 exitflag, output] = fsolve(@(x) bp_reg_solve_scenario(t,x,x_p0,pars, ...
                                                       fixed_var_pars, ...
                                                       drugs), ...
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

% save_data_name = sprintf('%s_ss_data.mat', gender{gg});
% save_data_name = strcat('Data/', save_data_name);
% save(save_data_name, 'SSdata', 'residual', 'exitflag', 'output')

% save_data_name = sprintf('%s_ss_data_male_sodreab.mat', gender{gg});
% save_data_name = strcat('Data/', save_data_name);
% save(save_data_name, 'SSdata', 'residual', 'exitflag', 'output')

% save_data_name = sprintf('%s_ss_data_female_sodreab.mat', gender{gg});
% save_data_name = strcat('Data/', save_data_name);
% save(save_data_name, 'SSdata', 'residual', 'exitflag', 'output')

% save_data_name = sprintf('%s_ss_data_male_raas.mat', gender{gg});
% save_data_name = strcat('Data/', save_data_name);
% save(save_data_name, 'SSdata', 'residual', 'exitflag', 'output')

% save_data_name = sprintf('%s_ss_data_new_sigmamyo.mat', gender{gg});
% save_data_name = strcat('Data/', save_data_name);
% save(save_data_name, 'SSdata', 'residual', 'exitflag', 'output')

% save_data_name = sprintf('%s_ss_data_new_Nadhs.mat', gender{gg});
% save_data_name = strcat('Data/', save_data_name);
% save(save_data_name, 'SSdata', 'residual', 'exitflag', 'output')

% save_data_name = sprintf('%s_ss_data_new_Phitwreab.mat', gender{gg});
% save_data_name = strcat('Data/', save_data_name);
% save(save_data_name, 'SSdata', 'residual', 'exitflag', 'output')

% save_data_name = sprintf('%s_ss_data_new_Psi.mat', gender{gg});
% save_data_name = strcat('Data/', save_data_name);
% save(save_data_name, 'SSdata', 'residual', 'exitflag', 'output')

%

save_data_name = sprintf('%s_ss_data_scenario_%s.mat', gender{gg},scenario{ss});
save_data_name = strcat('Data/', save_data_name);
save(save_data_name, 'SSdata', 'residual', 'exitflag', 'output')

% X(:,g) = x;
% RESIDUAL(:,g) = residual;
% EXITFLAG(g) = exitflag;
% OUTPUT{g} = output;

end % gender

end






























