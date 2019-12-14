% function [dfdy, dfdy_p] = jac_anal_find
function jac_anal_find

num_vars = 82;
num_pars = 39;

eqns   = sym('eqns', [num_vars,1]);
vars   = sym('vars', [num_vars,1]); 
vars_p = sym('vars', [num_vars,1]); 
pars   = sym('pars', [num_pars,1]); 

gender = {'male', 'female'};

for gg = 1:2 % gender

%% Retreive parameters by name.

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

N_rsna      = pars(1 );
R_aass      = pars(2 );
R_eass      = pars(3 );
P_B         = pars(4 );
P_go        = pars(5 );
C_gcf       = pars(6 );
eta_etapt   = pars(7 );
eta_epsdt   = pars(8 );
eta_etacd   = pars(9 );
K_vd        = pars(10);
K_bar       = pars(11);
R_bv        = pars(12);
T_adh       = pars(13);
Phi_sodin   = pars(14);
C_K         = pars(15);
T_al        = pars(16);
N_rs        = pars(17);
X_PRCPRA    = pars(18);
h_renin     = pars(19);
h_AGT       = pars(20);
h_AngI      = pars(21);
h_AngII     = pars(22);
h_Ang17     = pars(23);
h_AngIV     = pars(24);
h_AT1R      = pars(25);
h_AT2R      = pars(26);
c_GPautoreg = pars(27);
P_ghnom     = pars(28);
k_AGT       = pars(29);
c_ACE       = pars(30);
c_Chym      = pars(31);
c_NEP       = pars(32);
c_ACE2      = pars(33);
c_IIIV      = pars(34);
c_AT1R      = pars(35);
c_AT2R      = pars(36);
AT1R_eq     = pars(37);
AT2R_eq     = pars(38);
SF          = pars(39);

%% Retrieve variables by name.

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

rsna          = vars(1 );
alpha_map     = vars(2 );
alpha_rap     = vars(3 );
R_r           = vars(4 );
beta_rsna     = vars(5 );
Phi_rb        = vars(6 );
Phi_gfilt     = vars(7 );
P_f           = vars(8 );
P_gh          = vars(9 );
Sigma_tgf     = vars(10);
Phi_filsod    = vars(11);
Phi_ptsodreab = vars(12);
eta_ptsodreab = vars(13);
gamma_filsod  = vars(14);
gamma_at      = vars(15);
gamma_rsna    = vars(16);
Phi_mdsod     = vars(17);
Phi_dtsodreab = vars(18);
eta_dtsodreab = vars(19);
psi_al        = vars(20);
Phi_dtsod     = vars(21);
Phi_cdsodreab = vars(22);
eta_cdsodreab = vars(23);
lambda_dt     = vars(24);
lambda_anp    = vars(25);
Phi_usod      = vars(26);
Phi_win       = vars(27);
V_ecf         = vars(28);
V_b           = vars(29);
P_mf          = vars(30);
Phi_vr        = vars(31);
Phi_co        = vars(32);
P_ra          = vars(33);
vas           = vars(34);
vas_f         = vars(35);
vas_d         = vars(36);
R_a           = vars(37);
R_ba          = vars(38);
R_vr          = vars(39);
R_tp          = vars(40);
P_ma          = vars(41);
epsilon_aum   = vars(42);
a_auto        = vars(43);
a_chemo       = vars(44);
a_baro        = vars(45);
C_adh         = vars(46);
N_adh         = vars(47);
N_adhs        = vars(48);
delta_ra      = vars(49);
Phi_twreab    = vars(50);
mu_al         = vars(51);
mu_adh        = vars(52);
Phi_u         = vars(53);
M_sod         = vars(54);
C_sod         = vars(55);
nu_mdsod      = vars(56);
nu_rsna       = vars(57);
C_al          = vars(58);
N_al          = vars(59);
N_als         = vars(60);
xi_ksod       = vars(61);
xi_map        = vars(62);
xi_at         = vars(63);
hatC_anp      = vars(64);
AGT           = vars(65);
nu_AT1        = vars(66);
R_sec         = vars(67);
PRC           = vars(68);
PRA           = vars(69);
AngI          = vars(70);
AngII         = vars(71);
AT1R          = vars(72);
AT2R          = vars(73);
Ang17         = vars(74);
AngIV         = vars(75);
R_aa          = vars(76);
R_ea          = vars(77);
Sigma_myo     = vars(78);
Psi_AT1RAA    = vars(79);
Psi_AT1REA    = vars(80);
Psi_AT2RAA    = vars(81);
Psi_AT2REA    = vars(82);

vars_p(1 ) = 'rsna_p';          vars_p(2 ) = 'alpha_map_p';     vars_p(3 ) = 'alpha_rap_p'; 
vars_p(4 ) = 'R_r_p';           vars_p(5 ) = 'beta_rsna_p';     vars_p(6 ) = 'Phi_rb_p'; 
vars_p(7 ) = 'Phi_gfilt_p';     vars_p(8 ) = 'P_f_p';           vars_p(9 ) = 'P_gh_p'; 
vars_p(10) = 'Sigma_tgf_p';     vars_p(11) = 'Phi_filsod_p';    vars_p(12) = 'Phi_ptsodreab_p'; 
vars_p(13) = 'eta_ptsodreab_p'; vars_p(14) = 'gamma_filsod_p';  vars_p(15) = 'gamma_at_p'; 
vars_p(16) = 'gamma_rsna_p';    vars_p(17) = 'Phi_mdsod_p';     vars_p(18) = 'Phi_dtsodreab_p'; 
vars_p(19) = 'eta_dtsodreab_p'; vars_p(20) = 'psi_al_p';        vars_p(21) = 'Phi_dtsod_p'; 
vars_p(22) = 'Phi_cdsodreab_p'; vars_p(23) = 'eta_cdsodreab_p'; vars_p(24) = 'lambda_dt_p'; 
vars_p(25) = 'lambda_anp_p';    vars_p(26) = 'Phi_usod_p';      vars_p(27) = 'Phi_win_p'; 
vars_p(28) = 'V_ecf_p';         vars_p(29) = 'V_b_p';           vars_p(30) = 'P_mf_p'; 
vars_p(31) = 'Phi_vr_p';        vars_p(32) = 'Phi_co_p';        vars_p(33) = 'P_ra_p'; 
vars_p(34) = 'vas_p';           vars_p(35) = 'vas_f_p';         vars_p(36) = 'vas_d_p'; 
vars_p(37) = 'R_a_p';           vars_p(38) = 'R_ba_p';          vars_p(39) = 'R_vr_p'; 
vars_p(40) = 'R_tp_p';          vars_p(41) = 'P_ma_p';          vars_p(42) = 'epsilon_aum_p'; 
vars_p(43) = 'a_auto_p';        vars_p(44) = 'a_chemo_p';       vars_p(45) = 'a_baro_p'; 
vars_p(46) = 'C_adh_p';         vars_p(47) = 'N_adh_p';         vars_p(48) = 'N_adhs_p'; 
vars_p(49) = 'delta_ra_p';      vars_p(50) = 'Phi_twreab_p';    vars_p(51) = 'mu_al_p'; 
vars_p(52) = 'mu_adh_p';        vars_p(53) = 'Phi_u_p';         vars_p(54) = 'M_sod_p'; 
vars_p(55) = 'C_sod_p';         vars_p(56) = 'nu_mdsod_p';      vars_p(57) = 'nu_rsna_p'; 
vars_p(58) = 'C_al_p';          vars_p(59) = 'N_al_p';          vars_p(60) = 'N_als_p'; 
vars_p(61) = 'xi_ksod_p';       vars_p(62) = 'xi_map_p';        vars_p(63) = 'xi_at_p'; 
vars_p(64) = 'hatC_anp_p';      vars_p(65) = 'AGT_p';           vars_p(66) = 'nu_AT1_p'; 
vars_p(67) = 'R_sec_p';         vars_p(68) = 'PRC_p';           vars_p(69) = 'PRA_p'; 
vars_p(70) = 'AngI_p';          vars_p(71) = 'AngII_p';         vars_p(72) = 'AT1R_p'; 
vars_p(73) = 'AT2R_p';          vars_p(74) = 'Ang17_p';         vars_p(75) = 'AngIV_p'; 
vars_p(76) = 'R_aa_p';          vars_p(77) = 'R_ea_p';          vars_p(78) = 'Sigma_myo_p'; 
vars_p(79) = 'Psi_AT1RAA_p';    vars_p(80) = 'Psi_AT1REA_p';    vars_p(81) = 'Psi_AT2RAA_p'; 
vars_p(82) = 'Psi_AT2REA_p'; 

rsna_p          = vars_p(1 );
alpha_map_p     = vars_p(2 );
alpha_rap_p     = vars_p(3 );
R_r_p           = vars_p(4 );
beta_rsna_p     = vars_p(5 );
Phi_rb_p        = vars_p(6 );
Phi_gfilt_p     = vars_p(7 );
P_f_p           = vars_p(8 );
P_gh_p          = vars_p(9 );
Sigma_tgf_p     = vars_p(10);
Phi_filsod_p    = vars_p(11);
Phi_ptsodreab_p = vars_p(12);
eta_ptsodreab_p = vars_p(13);
gamma_filsod_p  = vars_p(14);
gamma_at_p      = vars_p(15);
gamma_rsna_p    = vars_p(16);
Phi_mdsod_p     = vars_p(17);
Phi_dtsodreab_p = vars_p(18);
eta_dtsodreab_p = vars_p(19);
psi_al_p        = vars_p(20);
Phi_dtsod_p     = vars_p(21);
Phi_cdsodreab_p = vars_p(22);
eta_cdsodreab_p = vars_p(23);
lambda_dt_p     = vars_p(24);
lambda_anp_p    = vars_p(25);
Phi_usod_p      = vars_p(26);
Phi_win_p       = vars_p(27);
V_ecf_p         = vars_p(28);
V_b_p           = vars_p(29);
P_mf_p          = vars_p(30);
Phi_vr_p        = vars_p(31);
Phi_co_p        = vars_p(32);
P_ra_p          = vars_p(33);
vas_p           = vars_p(34);
vas_f_p         = vars_p(35);
vas_d_p         = vars_p(36);
R_a_p           = vars_p(37);
R_ba_p          = vars_p(38);
R_vr_p          = vars_p(39);
R_tp_p          = vars_p(40);
P_ma_p          = vars_p(41);
epsilon_aum_p   = vars_p(42);
a_auto_p        = vars_p(43);
a_chemo_p       = vars_p(44);
a_baro_p        = vars_p(45);
C_adh_p         = vars_p(46);
N_adh_p         = vars_p(47);
N_adhs_p        = vars_p(48);
delta_ra_p      = vars_p(49);
Phi_twreab_p    = vars_p(50);
mu_al_p         = vars_p(51);
mu_adh_p        = vars_p(52);
Phi_u_p         = vars_p(53);
M_sod_p         = vars_p(54);
C_sod_p         = vars_p(55);
nu_mdsod_p      = vars_p(56);
nu_rsna_p       = vars_p(57);
C_al_p          = vars_p(58);
N_al_p          = vars_p(59);
N_als_p         = vars_p(60);
xi_ksod_p       = vars_p(61);
xi_map_p        = vars_p(62);
xi_at_p         = vars_p(63);
hatC_anp_p      = vars_p(64);
AGT_p           = vars_p(65);
nu_AT1_p        = vars_p(66);
R_sec_p         = vars_p(67);
PRC_p           = vars_p(68);
PRA_p           = vars_p(69);
AngI_p          = vars_p(70);
AngII_p         = vars_p(71);
AT1R_p          = vars_p(72);
AT2R_p          = vars_p(73);
Ang17_p         = vars_p(74);
AngIV_p         = vars_p(75);
R_aa_p          = vars_p(76);
R_ea_p          = vars_p(77);
Sigma_myo_p     = vars_p(78);
Psi_AT1RAA_p    = vars_p(79);
Psi_AT1REA_p    = vars_p(80);
Psi_AT2RAA_p    = vars_p(81);
Psi_AT2REA_p    = vars_p(82);

%% Equations

% rsna
% rsna0 = N_rsna * alpha_map * alpha_rap;
if     strcmp(gender{gg},'male')
    eqns(1 ) = rsna - N_rsna * alpha_map * alpha_rap;
elseif strcmp(gender{gg},'female')
    eqns(1 ) = rsna - (N_rsna * alpha_map * alpha_rap) ^ (1 / (N_rsna * alpha_map * alpha_rap));
end
% alpha_map
eqns(2 ) = alpha_map - ( 0.5 + 1 / (1 + exp((P_ma - 100) / 15)) );
% f(2 ) = alpha_map - ( 0.5 + 1.1 / (1 + exp((P_ma - 100) / 15)) );
% alpha_rap
eqns(3 ) = alpha_rap - ( 1 - 0.008 * P_ra );
% R_r
eqns(4 ) = R_r - ( R_aa + R_ea );
% beta_rsna
% f(5 ) = beta_rsna - ( 2 / (1 + exp(-3.16 * (rsna - 1))) );
eqns(5 ) = beta_rsna - ( 1.5 * (rsna - 1) + 1 );
% Phi_rb
eqns(6 ) = Phi_rb - ( P_ma / R_r );
% Phi_gfilt
eqns(7 ) = Phi_gfilt - ( P_f * C_gcf );
% P_f
eqns(8 ) = P_f - ( P_gh - (P_B + P_go) );
% P_gh
eqns(9 ) = P_gh - ( P_ma - Phi_rb * R_aa );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RAT
% Sigma_tgf
% f(10) = Sigma_tgf - ( 0.3408 + 3.449 / (3.88 + exp((Phi_mdsod - 3.859) / (-0.9617))) );
eqns(10) = Sigma_tgf - ( 0.3408 + 3.449 / (3.88 + exp((Phi_mdsod - 3.859 * SF) / (-0.9617 * SF) )) );
% Phi_filsod
eqns(11) = Phi_filsod - ( Phi_gfilt * C_sod );
% Phi_ptsodreab
eqns(12) = Phi_ptsodreab - ( Phi_filsod * eta_ptsodreab );
% eta_ptsodreab
eqns(13) = eta_ptsodreab - ( eta_etapt * gamma_filsod * gamma_at * gamma_rsna );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RAT
% gamma_filsod
% f(14) = gamma_filsod - ( 0.85 + 0.3 / (1 + exp((Phi_filsod - 18)/138)) );
eqns(14) = gamma_filsod - ( 0.85 + 0.3 / (1 + exp((Phi_filsod - 18 * SF)/(138 * SF) )) );
% gamma_at
eqns(15) = gamma_at - ( 0.95 + 0.12 / (1 + exp(2.6 - 1.8 * 1.301/20 * (AT1R*20/AT1R_eq))) );
% f(15) = gamma_at - ( 0.95 + 0.12 / (1 + exp(2.6 - 1.8 * log10(C_at))) );
% gamma_rsna
eqns(16) = gamma_rsna - ( 0.65+0.07 + 0.8*0.7 / (1 + exp((1 - rsna) / 2.18)) );
% f(16) = gamma_rsna - ( 0.5 + 0.7 / (1 + exp((1 - rsna) / 2.18)) );
% Phi_mdsod
eqns(17) = Phi_mdsod - ( Phi_filsod - Phi_ptsodreab );
% Phi_dtsodreab
eqns(18) = Phi_dtsodreab - ( Phi_mdsod * eta_dtsodreab );
% eta_dtsodreab
eqns(19) = eta_dtsodreab - ( eta_epsdt * psi_al );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RAT
% psi_al
% f(20) = psi_al - ( 0.17 + 0.94 / (1 + exp((0.48 - 1.2 * log10(C_al)) / 0.88)) );
eqns(20) = psi_al - ( 0.17 + 0.94 / (1 + exp((0.48 - 0.87 * log10(C_al)) / 0.88)) );
% Phi_dtsod
eqns(21) = Phi_dtsod - ( Phi_mdsod - Phi_dtsodreab );
% Phi_cdsodreab
eqns(22) = Phi_cdsodreab - ( Phi_dtsod * eta_cdsodreab );
% eta_cdsodreab
eqns(23) = eta_cdsodreab - ( eta_etacd * lambda_dt * lambda_anp );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RAT
% lambda_dt
% f(24) = lambda_dt - ( 0.82 + 0.39 / (1 + exp((Phi_dtsod - 1.7625) / 0.375)) );
eqns(24) = lambda_dt - ( 0.82 + 0.39 / (1 + exp((Phi_dtsod - 1.7625 * SF) / (0.375 * SF) )) );
% f(24) = lambda_dt - ( 0.82 + 0.39 / (1 + exp((Phi_dtsod - 1.6) / 2)) );
% lambda_anp
eqns(25) = lambda_anp - ( -0.1 * hatC_anp + 1.1 );
% f(25) = lambda_anp - ( -0.1 * hatC_anp + 1.1199 );
% Phi_usod
eqns(26) = Phi_usod - ( Phi_dtsod - Phi_cdsodreab );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RAT
% Phi_win
% f(27) = Phi_win - ( max( 0, 0.008 / (1 + 86.1*1.822*(C_adh^(3*-1.607))) - 0.005 ) );
% % f(27) = Phi_win - ( 0.008 / (1 + 1.822*(C_adh^(-1.607))) - 0.0053 );
% eqns(27) = Phi_win - ( max( 0, 0.008 * SF / (1 + 86.1*1.822*(C_adh^(3*-1.607))) - 0.005 * SF ) );
eqns(27) = Phi_win - ( 0.008 * SF / (1 + 86.1*1.822*(C_adh^(3*-1.607))) - 0.005 * SF );
% V_ecf
eqns(28) = V_ecf_p - ( Phi_win - Phi_u );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RAT
% V_b
% f(29) = V_b - ( 4.5479392962 + 2.431217 / (1 + exp(-(V_ecf - 18.11278) * 0.47437)) );
% % f(29) = V_b - ( 4.560227 + 2.431217 / (1 + exp(-(V_ecf - 18.11278) * 0.47437)) );
eqns(29) = V_b - ( 4.5479392962 * SF + 2.431217 * SF / (1 + exp(-(V_ecf - 18.11278 * SF) * (0.47437 / SF) )) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RAT
% P_mf
% f(30) = P_mf - ( (7.436 * V_b - 30.18) * epsilon_aum );
eqns(30) = P_mf - ( ( (7.436 / SF) * V_b - 30.18) * epsilon_aum );
% Phi_vr
eqns(31) = Phi_vr - ( (P_mf - P_ra) / R_vr );
% Phi_co
eqns(32) = Phi_co - ( Phi_vr );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RAT
% P_ra
% % f(33) = P_ra - ( 0.2787 * exp(Phi_co * 0.2281) -  );
% if     strcmp(gender,'male')
%     f(33) = P_ra - ( 0.2787 * exp(Phi_co * 0.2281) - 0.8268 );
% %     f(33) = P_ra - ( max( 0, 0.2787 * exp(Phi_co * 0.2281) - 0.8268 ) );
% elseif strcmp(gender,'female')
%     f(33) = P_ra - ( 0.2787 * exp(Phi_co * 0.2281) - 0.8245 );
% %     f(33) = P_ra - ( max( 0, 0.2787 * exp(Phi_co * 0.2281) - 0.8245 ) );
% end
if     strcmp(gender{gg},'male')
    eqns(33) = P_ra - ( 0.2787 * exp(Phi_co * 0.2281 / SF) - 0.8268 );
%     f(33) = P_ra - ( max( 0, 0.2787 * exp(Phi_co * 0.2281) - 0.8268 ) );
elseif strcmp(gender{gg},'female')
    eqns(33) = P_ra - ( 0.2787 * exp(Phi_co * 0.2281 / SF) - 0.8245 );
%     f(33) = P_ra - ( max( 0, 0.2787 * exp(Phi_co * 0.2281) - 0.8245 ) );
end
% vas
eqns(34) = vas_p - ( vas_f - vas_d );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RAT
% vas_f
% f(35) = vas_f - ( (11.312 * exp(-Phi_co * 0.4799)) / 100000 );
eqns(35) = vas_f - ( (11.312 * exp(-Phi_co * 0.4799 / SF)) / 100000 );
% vas_d
eqns(36) = vas_d - ( vas * K_vd );
% R_a
eqns(37) = R_a - ( R_ba * epsilon_aum );
% R_ba
eqns(38) = R_ba - ( K_bar / vas );
% R_vr
eqns(39) = R_vr - ( (8 * R_bv + R_a) / 31 );
% R_tp
eqns(40) = R_tp - ( R_a + R_bv );
% P_ma
eqns(41) = P_ma - ( Phi_co * R_tp );
% epsilon_aum
eqns(42) = epsilon_aum - ( a_chemo + a_baro ); 
% a_auto
eqns(43) = a_auto - ( 3.0042 * exp(-P_ma * 0.011) );
% f(43) = a_auto - ( 3.079 * exp(-P_ma * 0.011) );
% a_chemo
eqns(44) = a_chemo - ( 1/4 * a_auto );
% a_baro
eqns(45) = a_baro_p - ( 3/4 * (a_auto_p - 0.0000667 * (a_baro - 1)) );
% C_adh
eqns(46) = C_adh - ( 4 * N_adh );
% N_adh
eqns(47) = N_adh_p - ( 1/T_adh * (N_adhs - N_adh) );
% N_adhs
% eqns(48) = N_adhs - ( (C_sod - 141 + max( 1, epsilon_aum ) - 1 - delta_ra) / 3 );
f(48) = N_adhs - ( (C_sod - 141 + epsilon_aum - 1 - delta_ra) / 3 );
% delta_ra
eqns(49) = delta_ra_p - ( 0.2 * P_ra_p - 0.0007 * delta_ra );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RAT
% Phi_twreab
% f(50) = Phi_twreab - ( 0.025 - 0.001 / (mu_al * mu_adh) + 0.8 * Phi_gfilt );
eqns(50) = Phi_twreab - ( 0.025 * SF - 0.001 * SF / (mu_al * mu_adh) + 0.8 * Phi_gfilt );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RAT
% mu_al
% f(51) = mu_al - ( 0.17 + 0.94 / (1 + exp((0.48 - 1.2 * log10(C_al)) / 0.88)) );
eqns(51) = mu_al - ( 0.17 + 0.94 / (1 + exp((0.48 - 0.87 * log10(C_al)) / 0.88)) );
% mu_adh
eqns(52) = mu_adh - ( 0.3313 + 0.8 / (1 + exp(0.6 - 3.7 * log10(C_adh))) ); 
% f(52) = mu_adh - ( 0.37 + 0.8 / (1 + exp(0.6 - 3.7 * log10(C_adh))) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RAT
% Phi_u
% f(53) = Phi_u - ( max( 0.0003, Phi_gfilt - Phi_twreab ) );
% eqns(53) = Phi_u - ( max( 0.0003 * SF, Phi_gfilt - Phi_twreab ) );
f(53) = Phi_u - ( Phi_gfilt - Phi_twreab );
% M_sod
eqns(54) = M_sod_p - ( Phi_sodin - Phi_usod );
% C_sod
% eqns(55) = C_sod - ( max( 141, M_sod / V_ecf ) );
f(55) = C_sod - ( M_sod / V_ecf );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RAT
% nu_mdsod
% if     strcmp(gender,'male')
%     f(56) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.731) / 0.6056)) );
% elseif strcmp(gender,'female')
%     f(56) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.637) / 0.6056)) );
% end
% % f(56) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.667) / 0.6056)) );
if     strcmp(gender{gg},'male')
    eqns(56) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.731 * SF) / (0.6056 * SF) )) );
elseif strcmp(gender{gg},'female')
    eqns(56) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.637 * SF) / (0.6056 * SF) )) );
end
% nu_rsna
eqns(57) = nu_rsna - ( 1.822 - 2.056 / (1.358 + exp(rsna - 0.8667)) );
% f(57) = nu_rsna - ( 1.89 - 2.056 / (1.358 + exp(rsna - 0.8667)) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RAT
% C_al
% if     strcmp(gender,'male')
%     f(58) = C_al - ( max( 1, N_al * 85      ) );
% elseif strcmp(gender,'female')
%     f(58) = C_al - ( max( 1, N_al * 69.1775 ) );
% end
if     strcmp(gender{gg},'male')
%     eqns(58) = C_al - ( max( 1, N_al * 395.3 ) );
    eqns(58) = C_al - ( N_al * 395.3 );
elseif strcmp(gender{gg},'female')
%     eqns(58) = C_al - ( max( 1, N_al * 379.4 ) );
    eqns(58) = C_al - ( N_al * 379.4 );
end
% N_al
eqns(59) = N_al_p - ( 1/T_al * (N_als - N_al) );
% N_als
eqns(60) = N_als - ( xi_ksod * xi_map * xi_at );
% xi_ksod
% eqns(61) = xi_ksod - ( max( 0, (C_K / C_sod) / (C_K/144/(6+1)) - 6 ) ); 
eqns(61) = xi_ksod - ( (C_K / C_sod) / (C_K/144/(6+1)) - 6 ); 
% f(61) = xi_ksod - ( (C_K / C_sod) / 0.003525 - 9 );
% xi_map
% if P_ma <= 100
%     eqns(62) = xi_map - ( 70.1054 * exp(-0.0425 * P_ma) );
% %     f(62) = xi_map - ( 69.03 * exp(-0.0425 * P_ma) );
% else
%     eqns(62) = xi_map - ( 1 );
% end
eqns(62) = xi_map - ( 70.1054 * exp(-0.0425 * P_ma) );
% xi_at
eqns(63) = xi_at - ( 0.47 + 2.4 / (1 + exp((2.82 - 1.5 * 1.301/20 * (AT1R*20/AT1R_eq)) / 0.8)) );
% f(62) = xi_at - ( 0.4 + 2.4 / (1 + exp((2.82 - 1.5 * log10(C_at)) / 0.8)) );
% hatC_anp
eqns(64) = hatC_anp - ( 7.4052 - 6.554 / (1 + exp(P_ra - 3.762)) ); 
% f(63) = hatC_anp - ( 7.427 - 6.554 / (1 + exp(P_ra - 3.762)) );
% AGT
eqns(65) = AGT_p - ( k_AGT - PRA - log(2)/h_AGT * AGT );
% nu_AT1
eqns(66) = nu_AT1 - ( 10^(0.0102 - 0.95 * log10(AT1R / AT1R_eq)) );
% R_sec
eqns(67) = R_sec - ( N_rs * nu_mdsod * nu_rsna * nu_AT1 );
% PRC
eqns(68) = PRC_p - ( R_sec - log(2)/h_renin * PRC );
% PRA
eqns(69) = PRA - ( PRC * X_PRCPRA );
% AngI
eqns(70) = AngI_p - ( PRA - (c_ACE + c_Chym + c_NEP) * AngI - log(2)/h_AngI * AngI );
% AngII
eqns(71) = AngII_p - ( (c_ACE + c_Chym) * AngI - (c_ACE2 + c_IIIV + c_AT1R + c_AT2R) * AngII - log(2)/h_AngII * AngII );
% AT1R
eqns(72) = AT1R_p - ( c_AT1R * AngII - log(2)/h_AT1R * AT1R );
% AT2R
eqns(73) = AT2R_p - ( c_AT2R * AngII - log(2)/h_AT2R * AT2R );
% Ang17
eqns(74) = Ang17_p - ( c_NEP * AngI + c_ACE2 * AngII - log(2)/h_Ang17 * Ang17 );
% AngIV
eqns(75) = AngIV_p - ( c_IIIV * AngII - log(2)/h_AngIV * AngIV );
% R_aa
eqns(76) = R_aa - ( R_aass * beta_rsna * Sigma_tgf * Sigma_myo * Psi_AT1RAA * Psi_AT2RAA);
% R_ea
eqns(77) = R_ea - ( R_eass * Psi_AT1REA * Psi_AT2REA );
% Sigma_myo
% f(78) = Sigma_myo - ( 0.92 + 5 / ( 1 + exp(-2 * (P_gh - 64)) ) );
eqns(78) = Sigma_myo - ( c_GPautoreg * (P_gh / P_ghnom - 1) + 1 );
% Psi_AT1RAA
eqns(79) = Psi_AT1RAA - ( 0.8   + 0.1902*0.055 * (AT1R*20 / AT1R_eq) - 0.185 / (AT1R*20 / AT1R_eq) );
% Psi_AT1REA
eqns(80) = Psi_AT1REA - ( 0.925 + 0.0835*0.05  * (AT1R*20 / AT1R_eq) - 0.17  / (AT1R*20 / AT1R_eq) );
% Psi_AT2RAA
if     strcmp(gender{gg},'male')
    eqns(81) = Psi_AT2RAA - ( 1 );
elseif strcmp(gender{gg},'female')
    eqns(81) = Psi_AT2RAA - ( 0.025 * (AT2R_eq - AT2R) + 1 );
end
% Psi_AT2REA
if     strcmp(gender{gg},'male')
    eqns(82) = Psi_AT2REA - ( 1 );
elseif strcmp(gender{gg},'female')
    eqns(82) = Psi_AT2REA - ( 0.01  * (AT2R_eq - AT2R) + 1 );
end

%% Calculate Jacobian.

dfdy   = sym('dfdy'  , [num_vars,num_vars]);
dfdy_p = sym('dfdy_p', [num_vars,num_vars]);

tic
for i = 1:num_vars
    for j = 1:num_vars
        dfdy  (i,j) = diff(eqns(i), vars  (j));
        dfdy_p(i,j) = diff(eqns(i), vars_p(j));
    end
end
CPUtime = toc;

save_data_name = sprintf('%s_jac_anal.mat', gender{gg});
save_data_name = strcat('Data/', save_data_name);
save(save_data_name, 'dfdy', 'dfdy_p', 'CPUtime')

end % gender

end




























