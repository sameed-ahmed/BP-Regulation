function [dfdy_s,dfdy_p_s] = jac_spar

% dfdy_s   - Jacobian sparsity pattern wrt y
% dfdy_p_s - Jacobian sparsity pattern wrt y'

gender   = {'male', 'female'};
dfdy_s   = cell(1,2);
dfdy_p_s = cell(1,2);

num_pars = 39;
num_vars = 82;

pars_num = random('unif', 2,3, num_pars,1);
vars_num = random('unif', 2,3, num_vars,1);

for gg = 1:2 % gender

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

% Load analytic expression
if     strcmp(gender{gg}, 'male')
    load(  'male_jac_anal.mat', 'dfdy', 'dfdy_p');
elseif strcmp(gender{gg}, 'female')
    load('female_jac_anal.mat', 'dfdy', 'dfdy_p');
end

pars = sym('pars', [num_pars,1]); 
vars = sym('vars', [num_vars,1]); 

%% Set parameters by name.

% pars = [N_rsna; R_aass; R_eass; P_b; P_go; C_gcf; eta_etapt; eta_epsdt; ...
%         eta_etacd; K_vd; K_bar; R_bv; T_adh; Phi_sodin; C_k; T_al; ...
%         N_rs; X_PRC_PRA; h_renin; h_AGT; h_AngI; h_AngII; h_Ang17; ...
%         h_AngIV; h_AT1R; h_AT2R; c_gpautoreg; P_ghnom; k_AGT; c_ACE; ...
%         c_Chym; c_NEP; c_ACE2; c_IIIV; c_AT1R; c_AT2R; AT1R_eq; ...
%         AT2R_eq; gen];

pars(1 ) = 'N_rsna';      pars(2 ) = 'R_aazz';    pars(3 ) = 'R_eass';      
pars(4 ) = 'P_b';         pars(5 ) = 'P_go';      pars(6 ) = 'C_gcf';
pars(7 ) = 'eta_eta_pta'; pars(8 ) = 'eta_epsdt'; pars(9 ) = 'eta_etacd'; 
pars(10) = 'K_vd';        pars(11) = 'K_bar';     pars(12) = 'R_bv'; 
pars(13) = 'T_adh';       pars(14) = 'Phi_sodin'; pars(15) = 'C_k';         
pars(16) = 'T_al';        pars(17) = 'N_rs';      pars(18) = 'X_PRC_PRA'; 
pars(19) = 'h_renin';     pars(20) = 'h_AGT';     pars(21) = 'h_AngI';
pars(22) = 'h_AngII';     pars(23) = 'h_Ang17';   pars(24) = 'h_AngIV'; 
pars(25) = 'h_AT1R';      pars(26) = 'h_AT2R';    pars(27) = 'c_gpautoreg'; 
pars(28) = 'P_ghnom';     pars(29) = 'k_AGT';     pars(30) = 'c_ACE'; 
pars(31) = 'c_Chym';      pars(32) = 'c_NEP';     pars(33) = 'c_ACE2';
pars(34) = 'c_IIIV';      pars(35) = 'c_AT1R';    pars(36) = 'c_AT2R';
pars(37) = 'AT1R_eq';     pars(38) = 'AT2R_eq';   pars(39) = 'SF';

%% Set variables by name.

% vars  = [rsna; alpha_map; alpha_rap; R_r; beta_rsna; Phi_rb; Phi_gfilt; ...
%          P_f; P_gh; Sigma_tgf; Phi_filsod; Phi_ptsodreab; eta_ptsodreab; ...
%          gamma_filsod; gamma_at; gamma_rsna; Phi_mdsod; Phi_dtsodreab; ...
%          eta_dtsodreab; psi_al; Phi_dtsod; Phi_cdsodreab; eta_cdsodreab; ...
%          lambda_dt; lambda_anp; Phi_usod; Phi_win; V_ecf; V_b; P_mf; ...
%          Phi_vr; Phi_co; P_ra; vas; vas_f; vas_d; R_a; R_ba; R_vr; R_tp; ...
%          P_ma; epsilon_aum; a_auto; a_chemo; a_baro; C_adh; N_adh; ...
%          N_adhs; delta_ra; Phi_twreab; mu_al; mu_adh; Phi_u; M_sod; ...
%          C_sod; nu_mdsod; nu_rsna; C_al; N_al; N_als; xi_ksod; xi_map; ...
%          xi_at; hatC_anp; AGT; nu_AT1; R_sec; PRC; PRA; AngI; AngII; ...
%          AT1R; AT2R; Ang17; AngIV; R_aa; R_ea; Sigma_myo; Psi_AT1RAA; ...
%          Psi_AT1REA; Psi_AT2RAA; Psi_AT2REA];

vars(1 ) = 'rsna';          vars(2 ) = 'alpha_map';     vars(3 ) = 'alpha_rap'; 
vars(4 ) = 'R_r';           vars(5 ) = 'beta_rsna';     vars(6 ) = 'Phi_rb'; 
vars(7 ) = 'Phi_gfilt';     vars(8 ) = 'P_f';           vars(9 ) = 'P_gh'; 
vars(10) = 'Sigma_tgf';     vars(11) = 'Phi_filsod';    vars(12) = 'Phi_ptsodreab'; 
vars(13) = 'eta_ptsodreab'; vars(14) = 'gamma_filsod';  vars(15) = 'gamma_at'; 
vars(16) = 'gamma_rsna';    vars(17) = 'Phi_mdsod';     vars(18) = 'Phi_dtsodreab'; 
vars(19) = 'eta_dtsodreab'; vars(20) = 'psi_al';        vars(21) = 'Phi_dtsod'; 
vars(22) = 'Phi_cdsodreab'; vars(23) = 'eta_cdsodreab'; vars(24) = 'lambda_dt'; 
vars(25) = 'lambda_anp';    vars(26) = 'Phi_usod';      vars(27) = 'Phi_win'; 
vars(28) = 'V_ecf';         vars(29) = 'V_b';           vars(30) = 'P_mf'; 
vars(31) = 'Phi_vr';        vars(32) = 'Phi_co';        vars(33) = 'P_ra'; 
vars(34) = 'vas';           vars(35) = 'vas_f';         vars(36) = 'vas_d'; 
vars(37) = 'R_a';           vars(38) = 'R_ba';          vars(39) = 'R_vr'; 
vars(40) = 'R_tp';          vars(41) = 'P_ma';          vars(42) = 'epsilon_aum'; 
vars(43) = 'a_auto';        vars(44) = 'a_chemo';       vars(45) = 'a_baro'; 
vars(46) = 'C_adh';         vars(47) = 'N_adh';         vars(48) = 'N_adhs'; 
vars(49) = 'delta_ra';      vars(50) = 'Phi_twreab';    vars(51) = 'mu_al'; 
vars(52) = 'mu_adh';        vars(53) = 'Phi_u';         vars(54) = 'M_sod'; 
vars(55) = 'C_sod';         vars(56) = 'nu_mdsod';      vars(57) = 'nu_rsna'; 
vars(58) = 'C_al';          vars(59) = 'N_al';          vars(60) = 'N_als'; 
vars(61) = 'xi_ksod';       vars(62) = 'xi_map';        vars(63) = 'xi_at'; 
vars(64) = 'hatC_anp';      vars(65) = 'AGT';           vars(66) = 'nu_AT1'; 
vars(67) = 'R_sec';         vars(68) = 'PRC';           vars(69) = 'PRA'; 
vars(70) = 'AngI';          vars(71) = 'AngII';         vars(72) = 'AT1R'; 
vars(73) = 'AT2R';          vars(74) = 'Ang17';         vars(75) = 'AngIV'; 
vars(76) = 'R_aa';          vars(77) = 'R_ea';          vars(78) = 'Sigma_myo'; 
vars(79) = 'Psi_AT1RAA';    vars(80) = 'Psi_AT1REA';    vars(81) = 'Psi_AT2RAA'; 
vars(82) = 'Psi_AT2REA'; 

%% Plug values into analytic expression.

% Jacobian
dfdy   = double(subs(dfdy, [pars; vars], [pars_num; vars_num]));
dfdy_p = double(dfdy_p);

% Jacobian sparsity pattern
dfdy_s  {gg} = double(dfdy   ~= 0);
dfdy_p_s{gg} = double(dfdy_p ~= 0);

end % gender

end


























