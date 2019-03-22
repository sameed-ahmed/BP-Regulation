% This script calculates the steady state solution to the system using
% fsolve. Some previous values are used as an initial guess. These are
% taken from Jessica, which are taken in part from some paper (Karaaslan
% 2005?).

function solve_steady_state_numerical_with_fixed

gender     = {'male', 'female'};
SS_data_IG = zeros(82,2);
X          = zeros(82,2);
RESIDUAL       = zeros(82,2);
EXITFLAG   = zeros(1 ,2);
OUTPUT     = cell (1 ,2);

for g = 1:2

%% Parameters

if     strcmp(gender{g}, 'male')
    gen = 1;
elseif strcmp(gender{g}, 'female')
    gen = 0;
end

N_rsna      = 1;
R_aass      = 31.67;   % mmHg min / l
R_eass      = 51.66;   % mmHg min / l
P_B         = 18;      % mmHg
P_go        = 28;      % mmHg
C_gcf       = 0.00781;
eta_etapt   = 0.8; 
eta_epsdt   = 0.5; 
eta_etacd   = 0.93; 
K_vd        = 0.00001;
K_bar       = 16.6;    % mmHg min / l
R_bv        = 3.4;     % mmHg min / l
T_adh       = 6;       % min
Phi_sodin   = 0.126;   % mEq / min
C_K         = 5;       % mEq / l
T_al        = 30;      % min LISTED AS 30 IN TABLE %listed as 60 in text will only change dN_al
N_rs        = 1;       % ng / ml / min
X_PRCPRA    = 61/60.0; % fmol / min / pg
h_renin     = 12;      % min
h_AGT       = 10*60;   % min
h_AngI      = 0.5;     % min
h_AngII     = 0.66;    % min
h_Ang17     = 30;      % min
h_AngIV     = 0.5;     % min
h_AT1R      = 12;      % min
h_AT2R      = 12;      % min
c_GPautoreg = 5;
P_ghnom     = 62;      % mmHg

% Male and female different parameters
if     strcmp(gender{g}, 'male')
    k_AGT   = 577.04;
    c_ACE   = 0.88492;
    c_Chym  = 0.09315;
    c_NEP   = 0.038189;
    c_ACE2  = 0.0078009;
    c_IIIV  = 0.25056;
    c_AT1R  = 0.17008;
    c_AT2R  = 0.065667;
    AT1R_eq = 13.99;
    AT2R_eq = 5.0854;
elseif strcmp(gender{g}, 'female')
    k_AGT   = 610.39;
    c_ACE   = 1.4079;
    c_Chym  = 0.1482;
    c_NEP   = 0.060759;
    c_ACE2  = 0.0037603;
    c_IIIV  = 0.038644;
    c_AT1R  = 0.027089;
    c_AT2R  = 0.038699;
    AT1R_eq = 3.78;
    AT2R_eq = 5.0854;
end

pars = [N_rsna; R_aass; R_eass; P_B; P_go; C_gcf; eta_etapt; eta_epsdt; ...
        eta_etacd; K_vd; K_bar; R_bv; T_adh; Phi_sodin; C_K; T_al; ...
        N_rs; X_PRCPRA; h_renin; h_AGT; h_AngI; h_AngII; h_Ang17; ...
        h_AngIV; h_AT1R; h_AT2R; c_GPautoreg; P_ghnom; k_AGT; c_ACE; ...
        c_Chym; c_NEP; c_ACE2; c_IIIV; c_AT1R; c_AT2R; AT1R_eq; ...
        AT2R_eq; gen];

%% Variables initial guess

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

% Load data for steady state initial value. 
% Need to first run transform_data.m on Jessica's data files.
% if     strcmp(gender{g}, 'male')
%     SS_data_IG(:,g) = csvread(  'male_ss_data_IG.txt');
% elseif strcmp(gender{g}, 'female')
%     SS_data_IG(:,g) = csvread('female_ss_data_IG.txt');
% end
% if     strcmp(gender{g}, 'male')
%     load(  'male_ss_data.mat', 'SSdata');
% elseif strcmp(gender{g}, 'female')
%     load('female_ss_data.mat', 'SSdata');
% end
if     strcmp(gender{g}, 'male')
    load(  'male_sim_ss_high_data_Pma=206.mat'  , 'male_ss_data_sim'  );
elseif strcmp(gender{g}, 'female')
    load('female_sim_ss_high_data_Pma=204.9.mat', 'female_ss_data_sim');
end

% Order
% x  = [rsna; alpha_map; alpha_rap; R_r; beta_rsna; Phi_rb; Phi_gfilt; ...
%       P_f; P_gh; Sigma_tgf; Phi_filsod; Phi_ptsodreab; eta_ptsodreab; ...
%       gamma_filsod; gamma_at; gamma_rsna; Phi_mdsod; Phi_dtsodreab; ...
%       eta_dtsodreab; psi_al; Phi_dtsod; Phi_cdsodreab; eta_cdsodreab; ...
%       lambda_dt; lambda_anp; Phi_usod; Phi_win; V_ecf; V_b; P_mf; ...
%       Phi_vr; Phi_co; P_ra; vas; vas_f; vas_d; R_a; R_ba; R_vr; R_tp; ...
%       P_ma; epsilon_aum; a_auto; a_chemo; a_baro; C_adh; N_adh; ...
%       N_adhs; delta_ra; Phi_twreab; mu_al; mu_adh; Phi_u; M_sod; ...
%       C_sod; nu_mdsod; nu_rsna; C_al; N_al; N_als; xi_ksod; xi_map; ...
%       xi_at; hatC_anp; AGT; nu_AT1; R_sec; PRC; PRA; AngI; AngII; ...
%       AT1R; AT2R; Ang17; AngIV; R_aa; R_ea; Sigma_myo; Psi_AT1RAA; ...
%       Psi_AT1REA; Psi_AT2RAA; Psi_AT2REA];

% Initial guess for the variables.
% Find the steady state solution, so the derivative is 0.
% Arbitrary value for time to input.
% x0 = SS_data_IG(:,g); x_p0 = zeros(82,1); t = 0;
% x0 = SSdata;         x_p0 = zeros(82,1); t = 0;
if     strcmp(gender{g}, 'male')
    x0 = male_ss_data_sim;   x_p0 = zeros(82,1); t = 0;
elseif strcmp(gender{g}, 'female')
    x0 = female_ss_data_sim; x_p0 = zeros(82,1); t = 0;
end

% Variable to be fixed
if     strcmp(gender{g}, 'male')
    P_ma = 206.0;
%     P_ma = 101.448543237097;
%     P_ma = 115.0;
    
%     P_ma = 99.21;
%     P_ma = 88.6;
elseif strcmp(gender{g}, 'female')
    P_ma = 204.9;
%     P_ma = 100.194934694651;
%     P_ma = 113.9;
    
%     P_ma = 97.97;
%     P_ma = 87.4;
end
% Variable position
position = 41;
x0(position) = ''; x_p0(position) = '';

%% Find steady state solution

options = optimset(); %options = optimset('MaxFunEvals',8100+20000);
[SSdata, residual, ...
 exitflag, output] = fsolve(@(x) blood_press_reg_solve_with_fixed(t,x,x_p0,pars,P_ma), ...
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
    
if     P_ma < 100
    save_data_name = sprintf('%s_ss_low_data_with_fixed.mat', gender{g});
    save_data_name = strcat('Data/', save_data_name);
    save(save_data_name, 'SSdata', 'residual', 'exitflag', 'output')
elseif P_ma > 105
    save_data_name = sprintf('%s_ss_high_data_with_fixed.mat', gender{g});
    save_data_name = strcat('Data/', save_data_name);
    save(save_data_name, 'SSdata', 'residual', 'exitflag', 'output')
else
    save_data_name = sprintf('%s_ss_norm_data_with_fixed.mat', gender{g});
    save_data_name = strcat('Data/', save_data_name);
    save(save_data_name, 'SSdata', 'residual', 'exitflag', 'output')
end

% X(:,g) = x;
% RESIDUAL(:,g) = residual;
% EXITFLAG(g) = exitflag;
% OUTPUT{g} = output;

end

end






























