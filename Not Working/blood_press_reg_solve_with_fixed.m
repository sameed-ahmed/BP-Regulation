% This is a long-term model of the cardiovascular system accounting for the
% effects of renal sympathetic nervous activity (rsna) on kidney functions.
% It is adopted from:
% "Long-Term Mathematical Model Involving Renal Sympathetic Nerve Activity,
% Arterial Pressure, and Sodium Excretion" - 2005 - Karaaslan, et. al.
% 
% A sex-specific submodel for the renin angiotension system is
% incorporated. It is adopted from:
% "Sex-specific Long-term Blood Pressure Regulation: Modeling and Analysis"
% - 2018 - Leete, Layton.

% Differential algebraic equation system f(t,x(t),x'(t);theta,P_ma) = 0.

function f = blood_press_reg_solve_with_fixed(t,x,x_p,pars,P_ma)

%% Retrieve parameters by name.

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
gen         = pars(39);
if     gen == 1
    gender = 'male';
elseif gen == 0
    gender = 'female';
end

%% Retrieve variables by name.

rsna          = x(1 ); rsna_p          = x_p(1 ); 
alpha_map     = x(2 ); alpha_map_p     = x_p(2 ); 
alpha_rap     = x(3 ); alpha_rap_p     = x_p(3 ); 
R_r           = x(4 ); R_r_p           = x_p(4 ); 
beta_rsna     = x(5 ); beta_rsna_p     = x_p(5 ); 
Phi_rb        = x(6 ); Phi_rb_p        = x_p(6 ); 
Phi_gfilt     = x(7 ); Phi_gfilt_p     = x_p(7 ); 
P_f           = x(8 ); P_f_p           = x_p(8 ); 
P_gh          = x(9 ); P_gh_p          = x_p(9 ); 
Sigma_tgf     = x(10); Sigma_tgf_p     = x_p(10); 
Phi_filsod    = x(11); Phi_filsod_p    = x_p(11); 
Phi_ptsodreab = x(12); Phi_ptsodreab_p = x_p(12); 
eta_ptsodreab = x(13); eta_ptsodreab_p = x_p(13); 
gamma_filsod  = x(14); gamma_filsod_p  = x_p(14); 
gamma_at      = x(15); gamma_at_p      = x_p(15); 
gamma_rsna    = x(16); gamma_rsna_p    = x_p(16); 
Phi_mdsod     = x(17); Phi_mdsod_p     = x_p(17); 
Phi_dtsodreab = x(18); Phi_dtsodreab_p = x_p(18); 
eta_dtsodreab = x(19); eta_dtsodreab_p = x_p(19); 
psi_al        = x(20); psi_al_p        = x_p(20); 
Phi_dtsod     = x(21); Phi_dtsod_p     = x_p(21); 
Phi_cdsodreab = x(22); Phi_cdsodreab_p = x_p(22); 
eta_cdsodreab = x(23); eta_cdsodreab_p = x_p(23); 
lambda_dt     = x(24); lambda_dt_p     = x_p(24); 
lambda_anp    = x(25); lambda_anp_p    = x_p(25); 
Phi_usod      = x(26); Phi_usod_p      = x_p(26); 
Phi_win       = x(27); Phi_win_p       = x_p(27); 
V_ecf         = x(28); V_ecf_p         = x_p(28); 
V_b           = x(29); V_b_p           = x_p(29); 
P_mf          = x(30); P_mf_p          = x_p(30); 
Phi_vr        = x(31); Phi_vr_p        = x_p(31); 
Phi_co        = x(32); Phi_co_p        = x_p(32); 
P_ra          = x(33); P_ra_p          = x_p(33); 
vas           = x(34); vas_p           = x_p(34); 
vas_f         = x(35); vas_f_p         = x_p(35); 
vas_d         = x(36); vas_d_p         = x_p(36); 
R_a           = x(37); R_a_p           = x_p(37); 
R_ba          = x(38); R_ba_p          = x_p(38); 
R_vr          = x(39); R_vr_p          = x_p(39); 
R_tp          = x(40); R_tp_p          = x_p(40); 
% P_ma          = x(41); P_ma_p          = x_p(41); 
epsilon_aum   = x(42-1); epsilon_aum_p   = x_p(42-1); 
a_auto        = x(43-1); a_auto_p        = x_p(43-1); 
a_chemo       = x(44-1); a_chemo_p       = x_p(44-1); 
a_baro        = x(45-1); a_baro_p        = x_p(45-1); 
C_adh         = x(46-1); C_adh_p         = x_p(46-1); 
N_adh         = x(47-1); N_adh_p         = x_p(47-1); 
N_adhs        = x(48-1); N_adhs_p        = x_p(48-1); 
delta_ra      = x(49-1); delta_ra_p      = x_p(49-1); 
Phi_twreab    = x(50-1); Phi_twreab_p    = x_p(50-1); 
mu_al         = x(51-1); mu_al_p         = x_p(51-1); 
mu_adh        = x(52-1); mu_adh_p        = x_p(52-1); 
Phi_u         = x(53-1); Phi_u_p         = x_p(53-1); 
M_sod         = x(54-1); M_sod_p         = x_p(54-1); 
C_sod         = x(55-1); C_sod_p         = x_p(55-1); 
nu_mdsod      = x(56-1); nu_mdsod_p      = x_p(56-1); 
nu_rsna       = x(57-1); nu_rsna_p       = x_p(57-1); 
C_al          = x(58-1); C_al_p          = x_p(58-1); 
N_al          = x(59-1); N_al_p          = x_p(59-1); 
N_als         = x(60-1); N_als_p         = x_p(60-1); 
xi_ksod       = x(61-1); xi_ksod_p       = x_p(61-1); 
xi_map        = x(62-1); xi_map_p        = x_p(62-1); 
xi_at         = x(63-1); xi_at_p         = x_p(63-1); 
hatC_anp      = x(64-1); hatC_anp_p      = x_p(64-1); 
AGT           = x(65-1); AGT_p           = x_p(65-1); 
nu_AT1        = x(66-1); nu_AT1_p        = x_p(66-1); 
R_sec         = x(67-1); R_sec_p         = x_p(67-1); 
PRC           = x(68-1); PRC_p           = x_p(68-1); 
PRA           = x(69-1); PRA_p           = x_p(69-1); 
AngI          = x(70-1); AngI_p          = x_p(70-1); 
AngII         = x(71-1); AngII_p         = x_p(71-1); 
AT1R          = x(72-1); AT1R_p          = x_p(72-1); 
AT2R          = x(73-1); AT2R_p          = x_p(73-1); 
Ang17         = x(74-1); Ang17_p         = x_p(74-1); 
AngIV         = x(75-1); AngIV_p         = x_p(75-1); 
R_aa          = x(76-1); R_aa_p          = x_p(76-1); 
R_ea          = x(77-1); R_ea_p          = x_p(77-1); 
Sigma_myo     = x(78-1); Sigma_myo_p     = x_p(78-1); 
Psi_AT1RAA    = x(79-1); Psi_AT1RAA_p    = x_p(79-1); 
Psi_AT1REA    = x(80-1); Psi_AT1REA_p    = x_p(80-1); 
Psi_AT2RAA    = x(81-1); Psi_AT2RAA_p    = x_p(81-1); 
Psi_AT2REA    = x(82-1); Psi_AT2REA_p    = x_p(82-1); 

%% Differential algebraic equation system f(t,x,x') = 0.

f = zeros(length(x),1);

% rsna
rsna0 = N_rsna * alpha_map * alpha_rap;
if     strcmp(gender,'male')
    f(1 ) = rsna - rsna0;
elseif strcmp(gender,'female')
    f(1 ) = rsna - rsna0^(1/rsna0);
end
% alpha_map
f(2 ) = alpha_map - ( 0.5 + 1 / (1 + exp((P_ma - 100) / 15)) );
% f(2 ) = alpha_map - ( 0.5 + 1.1 / (1 + exp((P_ma - 100) / 15)) );
% alpha_rap
f(3 ) = alpha_rap - ( 1 - 0.008 * P_ra );
% R_r
f(4 ) = R_r - ( R_aa + R_ea );
% beta_rsna
f(5 ) = beta_rsna - ( 1.5 * (rsna - 1) + 1 );
% Phi_rb
f(6 ) = Phi_rb - ( P_ma / R_r );
% Phi_gfilt
f(7 ) = Phi_gfilt - ( P_f * C_gcf );
% P_f
f(8 ) = P_f - ( P_gh - (P_B + P_go) );
% P_gh
f(9 ) = P_gh - ( P_ma - Phi_rb * R_aa );
% Sigma_tgf
f(10) = Sigma_tgf - ( 0.3408 + 3.449 / (3.88 + exp((Phi_mdsod - 3.859) / (-0.9617))) );
% Phi_filsod
f(11) = Phi_filsod - ( Phi_gfilt * C_sod );
% Phi_ptsodreab
f(12) = Phi_ptsodreab - ( Phi_filsod * eta_ptsodreab );
% eta_ptsodreab
f(13) = eta_ptsodreab - ( eta_etapt * gamma_filsod * gamma_at * gamma_rsna );
% gamma_filsod
f(14) = gamma_filsod - ( 0.85 + 0.3 / (1 + exp((Phi_filsod - 18)/138)) );
% gamma_at
f(15) = gamma_at - ( 0.95 + 0.12 / (1 + exp(2.6 - 1.8 * 1.301/20 * (AT1R*20/AT1R_eq))) );
% f(15) = gamma_at - ( 0.95 + 0.12 / (1 + exp(2.6 - 1.8 * log10(C_at))) );
% gamma_rsna
f(16) = gamma_rsna - ( 0.65+0.07 + 0.8*0.7 / (1 + exp((1 - rsna) / 2.18)) );
% f(16) = gamma_rsna - ( 0.5 + 0.7 / (1 + exp((1 - rsna) / 2.18)) );
% Phi_mdsod
f(17) = Phi_mdsod - ( Phi_filsod - Phi_ptsodreab );
% Phi_dtsodreab
f(18) = Phi_dtsodreab - ( Phi_mdsod * eta_dtsodreab );
% eta_dtsodreab
f(19) = eta_dtsodreab - ( eta_epsdt * psi_al );
% psi_al
f(20) = psi_al - ( 0.17 + 1.442*0.94 / (1 + 5*exp((0.48 - 1.2 * log10(C_al)) / 0.88)) );
% f(20) = psi_al - ( 0.17 + 0.94 / (1 + exp((0.48 - 1.2 * log10(C_al)) / 0.88)) );
% Phi_dtsod
f(21) = Phi_dtsod - ( Phi_mdsod - Phi_dtsodreab );
% Phi_cdsodreab
f(22) = Phi_cdsodreab - ( Phi_dtsod * eta_cdsodreab );
% eta_cdsodreab
f(23) = eta_cdsodreab - ( eta_etacd * lambda_dt * lambda_anp );
% lambda_dt
f(24) = lambda_dt - ( 0.82 + 0.39 / (1 + exp((Phi_dtsod - 1.7625) / 0.375)) );
% f(24) = lambda_dt - ( 0.82 + 0.39 / (1 + exp((Phi_dtsod - 1.6) / 2)) );
% lambda_anp
f(25) = lambda_anp - ( -0.1 * hatC_anp + 1.1 );
% f(25) = lambda_anp - ( -0.1 * hatC_anp + 1.1199 );
% Phi_usod
f(26) = Phi_usod - ( Phi_dtsod - Phi_cdsodreab );
% Phi_win
f(27) = Phi_win - ( max( 0, 0.008 / (1 + 86.1*1.822*(C_adh^(3*-1.607))) - 0.005 ) );
% f(27) = Phi_win - ( 0.008 / (1 + 1.822*(C_adh^(-1.607))) - 0.0053 );
% V_ecf
f(28) = V_ecf_p - ( Phi_win - Phi_u );
% V_b
f(29) = V_b - ( 4.5479392962 + 2.431217 / (1 + exp(-(V_ecf - 18.11278) * 0.47437)) );
% f(29) = V_b - ( 4.560227 + 2.431217 / (1 + exp(-(V_ecf - 18.11278) * 0.47437)) );
% P_mf
f(30) = P_mf - ( (7.436 * V_b - 30.18) * epsilon_aum );
% Phi_vr
f(31) = Phi_vr - ( (P_mf - P_ra) / R_vr );
% Phi_co
f(32) = Phi_co - ( Phi_vr );
% P_ra
if     strcmp(gender,'male')
    f(33) = P_ra - ( 0.2787 * exp(Phi_co * 0.2281) - 0.8268 );
%     f(33) = P_ra - ( max( 0, 0.2787 * exp(Phi_co * 0.2281) - 0.8268 ) );
elseif strcmp(gender,'female')
    f(33) = P_ra - ( 0.2787 * exp(Phi_co * 0.2281) - 0.8245 );
%     f(33) = P_ra - ( max( 0, 0.2787 * exp(Phi_co * 0.2281) - 0.8245 ) );
end
% vas
f(34) = vas_p - ( vas_f - vas_d );
% vas_f
f(35) = vas_f - ( (11.312 * exp(-Phi_co * 0.4799)) / 100000 );
% vas_d
f(36) = vas_d - ( vas * K_vd );
% R_a
f(37) = R_a - ( R_ba * epsilon_aum );
% R_ba
f(38) = R_ba - ( K_bar / vas );
% R_vr
f(39) = R_vr - ( (8 * R_bv + R_a) / 31 );
% R_tp
f(40) = R_tp - ( R_a + R_bv );
% % P_ma
% f(41) = P_ma - ( Phi_co * R_tp );
% epsilon_aum
f(42-1) = epsilon_aum - ( a_chemo + a_baro ); 
% a_auto
f(43-1) = a_auto - ( 3.0042 * exp(-P_ma * 0.011) );
% f(43-1) = a_auto - ( 3.079 * exp(-P_ma * 0.011) );
% a_chemo
f(44-1) = a_chemo - ( 1/4 * a_auto );
% a_baro
f(45-1) = a_baro_p - ( 3/4 * (a_auto_p - 0.0000667 * (a_baro - 1)) );
% C_adh
f(46-1) = C_adh - ( 4 * N_adh );
% N_adh
f(47-1) = N_adh_p - ( 1/T_adh * (N_adhs - N_adh) );
% N_adhs
f(48-1) = N_adhs - ( (C_sod - 141 + max( 1, epsilon_aum ) - 1 - delta_ra) / 3 );
% f(48-1) = N_adhs - ( (C_sod - 141 + epsilon_aum - 1 - delta_ra) / 3 );
% delta_ra
f(49-1) = delta_ra_p - ( 0.2 * P_ra_p - 0.0007 * delta_ra );
% Phi_twreab
f(50-1) = Phi_twreab - ( 0.025 - 0.001 / (mu_al * mu_adh) + 0.8 * Phi_gfilt );
% mu_al
f(51-1) = mu_al - ( 0.17 + 0.94 / (1 + exp((0.48 - 1.2 * log10(C_al)) / 0.88)) );
% mu_adh
f(52-1) = mu_adh - ( 0.3313 + 0.8 / (1 + exp(0.6 - 3.7 * log10(C_adh))) ); 
% f(52-1) = mu_adh - ( 0.37 + 0.8 / (1 + exp(0.6 - 3.7 * log10(C_adh))) );
% Phi_u
f(53-1) = Phi_u - ( max( 0.0003, Phi_gfilt - Phi_twreab ) );
% f(53-1) = Phi_u - ( Phi_gfilt - Phi_twreab );
% M_sod
f(54-1) = M_sod_p - ( Phi_sodin - Phi_usod );
% C_sod
f(55-1) = C_sod - ( max( 141, M_sod / V_ecf ) );
% f(55-1) = C_sod - ( M_sod / V_ecf );
% nu_mdsod
if     strcmp(gender,'male')
    f(56-1) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.731) / 0.6056)) );
elseif strcmp(gender,'female')
    f(56-1) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.637) / 0.6056)) );
end
% f(56-1) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.667) / 0.6056)) );
% nu_rsna
f(57-1) = nu_rsna - ( 1.822 - 2.056 / (1.358 + exp(rsna - 0.8667)) );
% f(57-1) = nu_rsna - ( 1.89 - 2.056 / (1.358 + exp(rsna - 0.8667)) );
% C_al
if     strcmp(gender,'male')
    f(58-1) = C_al - ( max( 1, N_al * 85      ) );
elseif strcmp(gender,'female')
    f(58-1) = C_al - ( max( 1, N_al * 69.1775 ) );
end
% N_al
f(59-1) = N_al_p - ( 1/T_al * (N_als - N_al) );
% N_als
f(60-1) = N_als - ( xi_ksod * xi_map * xi_at );
% xi_ksod
f(61-1) = xi_ksod - ( max( 0, (C_K / C_sod) / (5/144/(6+1)) - 6 ) ); 
% f(60-1) = xi_ksod - ( (C_K / C_sod) / 0.003525 - 9 );
% xi_map
if P_ma <= 100
    f(62-1) = xi_map - ( 70.1054 * exp(-0.0425 * P_ma) );
%     f(62-1) = xi_map - ( 69.03 * exp(-0.0425 * P_ma) );
else
    f(62-1) = xi_map - ( 1 );
end
% xi_at
f(63-1) = xi_at - ( 0.47 + 2.4 / (1 + exp((2.82 - 1.5 * 1.301/20 * (AT1R*20/AT1R_eq)) / 0.8)) );
% f(62-1) = xi_at - ( 0.4 + 2.4 / (1 + exp((2.82 - 1.5 * log10(C_at)) / 0.8)) );
% hatC_anp
f(64-1) = hatC_anp - ( 7.4052 - 6.554 / (1 + exp(P_ra - 3.762)) ); 
% f(63-1) = hatC_anp - ( 7.427 - 6.554 / (1 + exp(P_ra - 3.762)) );
% AGT
f(65-1) = AGT_p - ( k_AGT - PRA - log(2)/h_AGT * AGT );
% nu_AT1
f(66-1) = nu_AT1 - ( 10^(0.0102 - 0.95 * log10(AT1R / AT1R_eq)) );
% R_sec
f(67-1) = R_sec - ( N_rs * nu_mdsod * nu_rsna * nu_AT1 );
% PRC
f(68-1) = PRC_p - ( R_sec - log(2)/h_renin * PRC );
% PRA
f(69-1) = PRA - ( PRC * X_PRCPRA );
% AngI
f(70-1) = AngI_p - ( PRA - (c_ACE + c_Chym + c_NEP) * AngI - log(2)/h_AngI * AngI );
% AngII
f(71-1) = AngII_p - ( (c_ACE + c_Chym) * AngI - (c_ACE2 + c_IIIV + c_AT1R + c_AT2R) * AngII - log(2)/h_AngII * AngII );
% AT1R
f(72-1) = AT1R_p - ( c_AT1R * AngII - log(2)/h_AT1R * AT1R );
% AT2R
f(73-1) = AT2R_p - ( c_AT2R * AngII - log(2)/h_AT2R * AT2R );
% Ang17
f(74-1) = Ang17_p - ( c_NEP * AngI + c_ACE2 * AngII - log(2)/h_Ang17 * Ang17 );
% AngIV
f(75-1) = AngIV_p - ( c_IIIV * AngII - log(2)/h_AngIV * AngIV );
% R_aa
f(76-1) = R_aa - ( R_aass * beta_rsna * Sigma_tgf * Sigma_myo * Psi_AT1RAA * Psi_AT2RAA);
% R_ea
f(77-1) = R_ea - ( R_eass * Psi_AT1REA * Psi_AT2REA );
% Sigma_myo
f(78-1) = Sigma_myo - ( c_GPautoreg * (P_gh / P_ghnom - 1) + 1 );
% Psi_AT1RAA
f(79-1) = Psi_AT1RAA - ( 0.8   + 0.1902*0.055 * (AT1R*20 / AT1R_eq) - 0.185 / (AT1R*20 / AT1R_eq) );
% Psi_AT1REA
f(80-1) = Psi_AT1REA - ( 0.925 + 0.0835*0.05  * (AT1R*20 / AT1R_eq) - 0.17  / (AT1R*20 / AT1R_eq) );
% Psi_AT2RAA
if     strcmp(gender,'male')
    f(81-1) = Psi_AT2RAA - ( 1 );
elseif strcmp(gender,'female')
    f(81-1) = Psi_AT2RAA - ( 0.025 * (AT2R_eq - AT2R) + 1 );
end
% Psi_AT2REA
if     strcmp(gender,'male')
    f(82-1) = Psi_AT2REA - ( 1 );
elseif strcmp(gender,'female')
    f(82-1) = Psi_AT2REA - ( 0.01  * (AT2R_eq - AT2R) + 1 );
end

end





























