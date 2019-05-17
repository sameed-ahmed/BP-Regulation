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

% Differential algebraic equation system f(t,x(t),x'(t);theta) = 0.

function f = bp_reg_sim_RPP(t,x,x_p,pars,fixed_var_pars,Phi_win_input,...
                            tchange,drugs,RPP,RPP_per,SSdata,scenario)

%% Scenarios

% Normal          - normal conditions
% Denerve         - cut off rsna from kidney
% Denerve & AT2R- - cut off rsna from kidney and block AT2R
% scenario = {'Normal', 'Denerve', 'Denerve & AT2R-'};

%% Retrieve drugs by name.

% drugs = [Ang II inf rate fmol/(ml min), ACEi target level, ARB target level]
if     t < tchange
    k_AngII   = 0; 
    gamma_ace = 0;
    gamma_arb = 0;
elseif t >= tchange
    k_AngII   = drugs(1);
    gamma_ace = drugs(2);
    gamma_arb = drugs(3);
end

%% Renal perfusion pressure.

% if     t < tchange
%     RPP = RPP;
% elseif t >= 1 *tchange && t < 2 *tchange
%     RPP = RPP + 1 *0.1;
% elseif t >= 2 *tchange && t < 3 *tchange
%     RPP = RPP + 2 *0.1;
% elseif t >= 3 *tchange && t < 4 *tchange
%     RPP = RPP + 3 *0.1;
% elseif t >= 4 *tchange && t < 5 *tchange
%     RPP = RPP + 4 *0.1;
% elseif t >= 5 *tchange && t < 6 *tchange
%     RPP = RPP + 5 *0.1;
% elseif t >= 6 *tchange && t < 7 *tchange
%     RPP = RPP + 6 *0.1;
% elseif t >= 7 *tchange && t < 8 *tchange
%     RPP = RPP + 7 *0.1;
% elseif t >= 8 *tchange && t < 9 *tchange
%     RPP = RPP + 8 *0.1;
% elseif t >= 9 *tchange && t < 10*tchange
%     RPP = RPP + 9 *0.1;
% elseif t >= 10*tchange 
%     RPP = RPP + 10*0.1;
% end

% if     t < tchange
%     RPP = RPP;
% elseif t >= 1 *tchange && t < 2 *tchange
%     RPP = RPP + 1 *10;
% elseif t >= 2 *tchange && t < 3 *tchange
%     RPP = RPP + 2 *10;
% elseif t >= 3 *tchange && t < 4 *tchange
%     RPP = RPP + 3 *10;
% elseif t >= 4 *tchange && t < 5 *tchange
%     RPP = RPP + 4 *10;
% elseif t >= 5 *tchange
%     RPP = RPP + 5 *10;
% end

if     t < tchange
    RPP = RPP;
elseif t >= 1 *tchange 
%     RPP = RPP;
%     RPP = RPP + RPP_per;
    RPP = RPP_per * tanh(5 * (t-tchange)) + RPP;
end

%% Retrieve parameters by name.

% Scaling factor
% Rat flow = Human flow x SF
SF = pars(end-1);
% Rat resistance = Human resistance x SF
SF_R = pars(end);

N_rsna           = pars(1 );
R_aass           = pars(2 );
R_eass           = pars(3 );
P_B              = pars(4 );
P_go             = pars(5 );
C_gcf            = pars(6 );
eta_ptsodreab_eq = pars(7 );
eta_dtsodreab_eq = pars(8 );
eta_cdsodreab_eq = pars(9 );
eta_ptwreab_eq   = pars(10);
eta_dtwreab_eq   = pars(11);
eta_cdwreab_eq   = pars(12);
K_vd             = pars(13);
K_bar            = pars(14);
R_bv             = pars(15);
T_adh            = pars(16);
Phi_sodin        = pars(17);
C_K              = pars(18);
T_al             = pars(19);
N_rs             = pars(20);
X_PRCPRA         = pars(21);
h_renin          = pars(22);
h_AGT            = pars(23);
h_AngI           = pars(24);
h_AngII          = pars(25);
h_Ang17          = pars(26);
h_AngIV          = pars(27);
h_AT1R           = pars(28);
h_AT2R           = pars(29);
k_AGT            = pars(30);
c_ACE            = pars(31)*(1-gamma_ace);
c_Chym           = pars(32);
c_NEP            = pars(33);
c_ACE2           = pars(34);
c_IIIV           = pars(35);
c_AT1R           = pars(36)*(1-gamma_arb);
c_AT2R           = pars(37);
AT1R_eq          = pars(38);
AT2R_eq          = pars(39);
gen              = pars(40);
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
lambda_al     = x(26); lambda_al_p     = x_p(26); 
Phi_usod      = x(27); Phi_usod_p      = x_p(27); 
% Phi_win       = x(28); Phi_win_p       = x_p(28); 
V_ecf         = x(29-1); V_ecf_p         = x_p(29-1); 
V_b           = x(30-1); V_b_p           = x_p(30-1); 
P_mf          = x(31-1); P_mf_p          = x_p(31-1); 
Phi_vr        = x(32-1); Phi_vr_p        = x_p(32-1); 
Phi_co        = x(33-1); Phi_co_p        = x_p(33-1); 
P_ra          = x(34-1); P_ra_p          = x_p(34-1); 
vas           = x(35-1); vas_p           = x_p(35-1); 
vas_f         = x(36-1); vas_f_p         = x_p(36-1); 
vas_d         = x(37-1); vas_d_p         = x_p(37-1); 
R_a           = x(38-1); R_a_p           = x_p(38-1); 
R_ba          = x(39-1); R_ba_p          = x_p(39-1); 
R_vr          = x(40-1); R_vr_p          = x_p(40-1); 
R_tp          = x(41-1); R_tp_p          = x_p(41-1); 
P_ma          = x(42-1); P_ma_p          = x_p(42-1); 
epsilon_aum   = x(43-1); epsilon_aum_p   = x_p(43-1); 
a_auto        = x(44-1); a_auto_p        = x_p(44-1); 
a_chemo       = x(45-1); a_chemo_p       = x_p(45-1); 
a_baro        = x(46-1); a_baro_p        = x_p(46-1); 
C_adh         = x(47-1); C_adh_p         = x_p(47-1); 
N_adh         = x(48-1); N_adh_p         = x_p(48-1); 
N_adhs        = x(49-1); N_adhs_p        = x_p(49-1); 
delta_ra      = x(50-1); delta_ra_p      = x_p(50-1); 
Phi_ptwreab   = x(51-1); Phi_ptwreab_p   = x_p(51-1); 
eta_ptwreab   = x(52-1); eta_ptwreab_p   = x_p(52-1); 
mu_ptsodreab  = x(53-1); mu_ptsodreab_p  = x_p(53-1); 
Phi_mdu       = x(54-1); Phi_mdu_p       = x_p(54-1); 
Phi_dtwreab   = x(55-1); Phi_dtwreab_p   = x_p(55-1); 
eta_dtwreab   = x(56-1); eta_dtwreab_p   = x_p(56-1); 
mu_dtsodreab  = x(57-1); mu_dtsodreab_p  = x_p(57-1); 
Phi_dtu       = x(58-1); Phi_dtu_p       = x_p(58-1); 
Phi_cdwreab   = x(59-1); Phi_cdwreab_p   = x_p(59-1); 
eta_cdwreab   = x(60-1); eta_cdwreab_p   = x_p(60-1); 
mu_cdsodreab  = x(61-1); mu_cdsodreab_p  = x_p(61-1); 
mu_adh        = x(62-1); mu_adh_p        = x_p(62-1); 
Phi_u         = x(63-1); Phi_u_p         = x_p(63-1); 
M_sod         = x(64-1); M_sod_p         = x_p(64-1); 
C_sod         = x(65-1); C_sod_p         = x_p(65-1); 
nu_mdsod      = x(66-1); nu_mdsod_p      = x_p(66-1); 
nu_rsna       = x(67-1); nu_rsna_p       = x_p(67-1); 
C_al          = x(68-1); C_al_p          = x_p(68-1); 
N_al          = x(69-1); N_al_p          = x_p(69-1); 
N_als         = x(70-1); N_als_p         = x_p(70-1); 
xi_ksod       = x(71-1); xi_ksod_p       = x_p(71-1); 
xi_map        = x(72-1); xi_map_p        = x_p(72-1); 
xi_at         = x(73-1); xi_at_p         = x_p(73-1); 
hatC_anp      = x(74-1); hatC_anp_p      = x_p(74-1); 
AGT           = x(75-1); AGT_p           = x_p(75-1); 
nu_AT1        = x(76-1); nu_AT1_p        = x_p(76-1); 
R_sec         = x(77-1); R_sec_p         = x_p(77-1); 
PRC           = x(78-1); PRC_p           = x_p(78-1); 
PRA           = x(79-1); PRA_p           = x_p(79-1); 
AngI          = x(80-1); AngI_p          = x_p(80-1); 
AngII         = x(81-1); AngII_p         = x_p(81-1); 
AT1R          = x(82-1); AT1R_p          = x_p(82-1); 
AT2R          = x(83-1); AT2R_p          = x_p(83-1); 
Ang17         = x(84-1); Ang17_p         = x_p(84-1); 
AngIV         = x(85-1); AngIV_p         = x_p(85-1); 
R_aa          = x(86-1); R_aa_p          = x_p(86-1); 
R_ea          = x(87-1); R_ea_p          = x_p(87-1); 
Sigma_myo     = x(88-1); Sigma_myo_p     = x_p(88-1); 
Psi_AT1RAA    = x(89-1); Psi_AT1RAA_p    = x_p(89-1); 
Psi_AT1REA    = x(90-1); Psi_AT1REA_p    = x_p(90-1); 
Psi_AT2RAA    = x(91-1); Psi_AT2RAA_p    = x_p(91-1); 
Psi_AT2REA    = x(92-1); Psi_AT2REA_p    = x_p(92-1); 

%% Differential algebraic equation system f(t,x,x') = 0.

f = zeros(length(x),1);

% rsna
rsna0 = N_rsna * alpha_map * alpha_rap;
if     strcmp(gender,'male')
    f(1 ) = rsna - rsna0;
elseif strcmp(gender,'female')
    f(1 ) = rsna - rsna0^(1/rsna0);
%     f(1 ) = rsna - rsna0; % male
end
% alpha_map
f(2 ) = alpha_map - ( 0.5 + 1 / (1 + exp((P_ma - fixed_var_pars(1)) / 15)) );
% alpha_rap
f(3 ) = alpha_rap - ( 1 - 0.008 * P_ra );
% R_r
f(4 ) = R_r - ( R_aa + R_ea );
% MODIFY ------------------------------------------------------------------
% beta_rsna
if     t < tchange
    f(5 ) = beta_rsna - ( 2 / (1 + exp(-3.16 * (rsna - 1))) );
elseif t >= tchange
    if   strcmp(scenario,'Denerve') || strcmp(scenario,'Denerve & AT2R-')
%         f(5 ) = beta_rsna - ( SSdata(5) );
        f(5 ) = beta_rsna - ( 1 );
    else
        f(5 ) = beta_rsna - ( 2 / (1 + exp(-3.16 * (rsna - 1))) );
    end
end
% MODIFY ------------------------------------------------------------------ - 1) + 1 );
% VARY --------------------------------------------------------------------
% Phi_rb
if     t < tchange
    f(6 ) = Phi_rb - ( P_ma / R_r );
elseif t >= tchange
    f(6 ) = Phi_rb - ( RPP / R_r );
end
% VARY --------------------------------------------------------------------
% Phi_gfilt
f(7 ) = Phi_gfilt - ( P_f * C_gcf );
% P_f
f(8 ) = P_f - ( P_gh - (P_B + P_go) );
% VARY --------------------------------------------------------------------
% P_gh
if     t < tchange
    f(9 ) = P_gh - ( P_ma - Phi_rb * R_aa );
elseif t >= tchange
    f(9 ) = P_gh - ( RPP - Phi_rb * R_aa );
end
% VARY --------------------------------------------------------------------
% Sigma_tgf - rat - female reabsorption
% f(10) = Sigma_tgf - ( 0.3408 + 3.449 / (3.88 + exp((Phi_mdsod - 3.859) / (-0.9617))) );
if     strcmp(gender,'male')
    f(10) = Sigma_tgf - ( 0.3408 + 3.449 / (3.88 + exp((Phi_mdsod - fixed_var_pars(2)) / (-0.9617 * SF) )) );
%     f(10) = Sigma_tgf - ( 0.3408 + 3.449 / (3.88 + exp((Phi_mdsod - 3.859 * SF * 2.500) / (-0.9617 * SF * 2.500) )) ); % female
elseif strcmp(gender,'female')
    f(10) = Sigma_tgf - ( 0.3408 + 3.449 / (3.88 + exp((Phi_mdsod - 3.859 * SF * 2.500) / (-0.9617 * SF * 2.500) )) );
%     f(10) = Sigma_tgf - ( 0.3408 + 3.449 / (3.88 + exp((Phi_mdsod - 3.859 * SF) / (-0.9617 * SF) )) ); % male
end
% Phi_filsod
f(11) = Phi_filsod - ( Phi_gfilt * C_sod );
% Phi_ptsodreab
f(12) = Phi_ptsodreab - ( Phi_filsod * eta_ptsodreab );
% eta_ptsodreab
f(13) = eta_ptsodreab - ( eta_ptsodreab_eq * gamma_filsod * gamma_at * gamma_rsna );
% gamma_filsod - rat
% f(14) = gamma_filsod - ( 0.85 + 0.3 / (1 + exp((Phi_filsod - 18)/138)) );
f(14) = gamma_filsod - ( 0.85 + 0.3 / (1 + exp((Phi_filsod - fixed_var_pars(3))/(138 * SF) )) );
% gamma_at
f(15) = gamma_at - ( 0.95 + 0.12 / (1 + exp(2.6785 - 2.342 * (AT1R/AT1R_eq))) );
% MODIFY ------------------------------------------------------------------
% gamma_rsna
if     t < tchange
    f(16) = gamma_rsna - ( 0.72 + 0.56 / (1 + exp((1 - rsna) / 2.18)) );
elseif t >= tchange
    if   strcmp(scenario,'Denerve') || strcmp(scenario,'Denerve & AT2R-')
%         f(16) = gamma_rsna - ( SSdata(16) );
        f(16) = gamma_rsna - ( 1 );
    else
        f(16) = gamma_rsna - ( 0.72 + 0.56 / (1 + exp((1 - rsna) / 2.18)) );
    end
end
% MODIFY ------------------------------------------------------------------
% Phi_mdsod
f(17) = Phi_mdsod - ( Phi_filsod - Phi_ptsodreab );
% Phi_dtsodreab
f(18) = Phi_dtsodreab - ( Phi_mdsod * eta_dtsodreab );
% eta_dtsodreab
f(19) = eta_dtsodreab - ( eta_dtsodreab_eq * psi_al );
% psi_al - rat
% f(20) = psi_al - ( 0.17 + 0.94 / (1 + exp((0.48 - 1.2 * log10(C_al)) / 0.88)) );
if     strcmp(gender,  'male')
%     f(20) = psi_al - ( 1/(395^0.3) * C_al^0.3 );
%     % ------------------------------------------------------
%     h1 = 2 / eta_dtsodreab_eq;
%     g1 = h1 / 2;
%     f1 = -1/395 * log(h1 / (1 + g1) - 1);
%     f(20) = psi_al - ( -g1 + h1 / (1 + exp(-f1 * C_al)) );
%     % ------------------------------------------------------
    % ------------------------------------------------------
    dd = 0.3 / eta_dtsodreab_eq;
    aa = 1 / eta_dtsodreab_eq - dd;
    bb = 1.34;
    cc = (aa * 395^bb) / (1 - dd) - 395^bb;
    f(20) = psi_al - ( (aa * C_al^bb) / (cc + C_al^bb) + dd );
    % ------------------------------------------------------
%     f(20) = psi_al - 1;
elseif strcmp(gender,'female')
    f(20) = psi_al - ( 1/(379^0.3) * C_al^0.3 );
end
% Phi_dtsod
f(21) = Phi_dtsod - ( Phi_mdsod - Phi_dtsodreab );
% Phi_cdsodreab
f(22) = Phi_cdsodreab - ( Phi_dtsod * eta_cdsodreab );
% eta_cdsodreab
f(23) = eta_cdsodreab - ( eta_cdsodreab_eq * lambda_dt * lambda_anp );
% lambda_dt - rat - female reabsorption
% f(24) = lambda_dt - ( 0.82 + 0.39 / (1 + exp((Phi_dtsod - 1.7625) / 0.375)) );
if     strcmp(gender,'male')
%     f(24) = lambda_dt - ( 0.82 + 0.39 / (1 + exp((Phi_dtsod - 1.7625 * SF) / (0.375 * SF) )) );
    f(24) = lambda_dt - ( 0.82 + 0.2553 / (1 + exp((Phi_dtsod - fixed_var_pars(4)) / (0.245 * SF) )) );
%     f(24) = lambda_dt - ( 0.82 + 0.2109 / (1 + exp((Phi_dtsod - 2.22 * SF * 2.504) / (0.224 * SF * 2.504) )) ); % female
elseif strcmp(gender,'female')
%     f(24) = lambda_dt - ( 0.82 + 0.39 / (1 + exp((Phi_dtsod - 1.7625 * SF * 2.504) / (0.375 * SF * 2.504) )) );
    f(24) = lambda_dt - ( 0.82 + 0.2109 / (1 + exp((Phi_dtsod - 2.22 * SF * 2.504) / (0.224 * SF * 2.504) )) );
%     f(24) = lambda_dt - ( 0.82 + 0.2553 / (1 + exp((Phi_dtsod - 2.03 * SF) / (0.245 * SF) )) ); % male
end
% lambda_anp
f(25) = lambda_anp - ( -0.1 * hatC_anp + 1.1 );
% lambda_al
if     strcmp(gender,  'male')
    f(26) = lambda_al - ( 1/(395^0.06) * C_al^0.06 );
%     f(26) = lambda_al - ( 1/(395^0.2) * C_al^0.2 );
%     f(26) = lambda_al - ( 1 );
elseif strcmp(gender,'female')
    f(26) = lambda_al - ( 1/(379^0.06) * C_al^0.06 );
end
% Phi_usod
f(27) = Phi_usod - ( Phi_dtsod - Phi_cdsodreab );
% FIXED -------------------------------------------------------------------
% % Phi_win - rat
% % f(28) = Phi_win - ( 0.003 / (1 + exp(-2.25 * (C_adh - 3.87))) );
% f(28) = Phi_win - ( 0.003 * SF / (1 + exp(-2.25 * (C_adh - 3.87))) );
% FIXED -------------------------------------------------------------------
% V_ecf
f(29-1) = V_ecf_p - ( Phi_win_input - Phi_u );
% V_b - rat
% f(30-1) = V_b - ( 4.5479 + 2.4312 / (1 + exp(-(V_ecf - 18.1128) * 0.4744)) );
f(30-1) = V_b - ( 4.5479+10 + 2.4312 / (1 + exp(-(V_ecf-30 - 18.1128) * (0.4744) )) );
% P_mf - rat
% f(31-1) = P_mf - ( (7.436 * V_b - 30.18) * epsilon_aum );
pmfpmf = (7.4360/3);
f(31-1) = P_mf - ( ( pmfpmf * V_b - 30.18) * epsilon_aum );
% Phi_vr
f(32-1) = Phi_vr - ( (P_mf - P_ra) / R_vr );
% Phi_co
f(33-1) = Phi_co - ( Phi_vr );
% P_ra - rat
% f(34-1) = P_ra - ( max( 0, 0.2787 * exp(Phi_co * 0.2281) - 0.8256 ) );
prapra = 0.2787 * exp(fixed_var_pars(end) * 0.2281 * SF_R);
f(34-1) = P_ra - ( max( 0, 0.2787 * exp(Phi_co * 0.2281 * SF_R) - prapra ) );
% vas
f(35-1) = vas_p - ( vas_f - vas_d );
% vas_f - rat
% f(36-1) = vas_f - ( (11.312 * exp(-Phi_co * 0.4799)) / 100000 );
vvv = -1/fixed_var_pars(end) * log(1/11.312);
f(36-1) = vas_f - ( (11.312 * exp(-Phi_co * vvv)) / 100000 );
% vas_d
f(37-1) = vas_d - ( vas * K_vd );
% R_a
f(38-1) = R_a - ( R_ba * epsilon_aum );
% R_ba
f(39-1) = R_ba - ( K_bar / vas );
% R_vr
f(40-1) = R_vr - ( (8 * R_bv + R_a) / 31 );
% R_tp
f(41-1) = R_tp - ( R_a + R_bv );
% P_ma
f(42-1) = P_ma - ( Phi_co * R_tp );
% epsilon_aum
f(43-1) = epsilon_aum - ( 4/5 * (a_chemo + a_baro) );
% a_auto
f(44-1) = a_auto - ( 3.0042 * exp(-fixed_var_pars(5) * P_ma) );
% a_chemo
f(45-1) = a_chemo - ( 1/4 * a_auto );
% a_baro
f(46-1) = a_baro_p - ( 3/4 * (a_auto_p - 0.0000667 * (a_baro - 1)) );
% C_adh
f(47-1) = C_adh - ( 4 * N_adh );
% N_adh
f(48-1) = N_adh_p - ( 1/T_adh * (N_adhs - N_adh) );
% N_adhs
% f(49-1) = N_adhs - ( (C_sod - 141 + max( 0, epsilon_aum - 1 ) - delta_ra) / 3 );
f(49-1) = N_adhs - ( (max( 0, C_sod - fixed_var_pars(6)) + max( 0, epsilon_aum - 1 ) - delta_ra) / 3 );
% delta_ra
f(50-1) = delta_ra_p - ( 0.2 * P_ra_p - 0.0007 * delta_ra );

% Phi_ptwreab
f(51-1) = Phi_ptwreab - ( Phi_gfilt * eta_ptwreab );
% eta_ptwreab
f(52-1) = eta_ptwreab - ( eta_ptwreab_eq * mu_ptsodreab );
% mu_ptsodreab
f(53-1) = mu_ptsodreab - ( 0.5 * 7/43 * tanh(13 * (eta_ptsodreab/eta_ptsodreab_eq - 1)) + 1 );
% f(53-1) = mu_ptsodreab - ( 1 );
% Phi_mdu
f(54-1) = Phi_mdu - ( Phi_gfilt - Phi_ptwreab );
% Phi_dtwreab
f(55-1) = Phi_dtwreab - ( Phi_mdu * eta_dtwreab );
% eta_dtwreab
f(56-1) = eta_dtwreab - ( eta_dtwreab_eq * mu_dtsodreab );
% mu_dtsodreab
f(57-1) = mu_dtsodreab - ( 0.5 * 2/3 * tanh(3.2 * (eta_dtsodreab/eta_dtsodreab_eq - 1)) + 1 );
% f(57-1) = mu_dtsodreab - ( 1 );
% Phi_dtu
f(58-1) = Phi_dtu - ( Phi_mdu - Phi_dtwreab );
% Phi_cdwreab
f(59-1) = Phi_cdwreab - ( Phi_dtu * eta_cdwreab );
% eta_cdwreab
f(60-1) = eta_cdwreab - ( eta_cdwreab_eq * mu_cdsodreab * mu_adh );
% mu_cdsodreab
f(61-1) = mu_cdsodreab - ( 0.5 * 11/39 * tanh(9.7 * (eta_cdsodreab/eta_cdsodreab_eq - 1)) + 1 );
% f(61-1) = mu_cdsodreab - ( 1 );
% mu_adh
aaa = 1.0328;
bbb = 0.1938;
ccc = -1/4 * log((aaa - 1) / bbb);
f(62-1) = mu_adh - ( aaa - bbb * exp(-ccc * C_adh) );
% Phi_u - rat
% f(63-1) = Phi_u - ( max( 0.0003, Phi_gfilt - Phi_twreab ) );
f(63-1) = Phi_u - ( Phi_dtu - Phi_cdwreab );

% M_sod
f(64-1) = M_sod_p - ( Phi_sodin - Phi_usod );
% C_sod
f(65-1) = C_sod - ( M_sod / V_ecf );
% nu_mdsod - rat - female reabsorption
% if     strcmp(gender,'male')
%     f(66-1) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.731) / 0.6056)) );
% elseif strcmp(gender,'female')
%     f(66-1) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.637) / 0.6056)) );
% end
% % f(66-1) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.667) / 0.6056)) );
if     strcmp(gender,'male')
    f(66-1) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - fixed_var_pars(7)) / (0.6056 * SF) )) );
%     f(66-1) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.637 * SF * 2.500) / (0.6056 * SF * 2.500) )) ); % female
elseif strcmp(gender,'female')
    f(66-1) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.637 * SF * 2.500) / (0.6056 * SF * 2.500) )) );
%     f(66-1) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.731 * SF) / (0.6056 * SF) )) ); % male
end
% MODIFY ------------------------------------------------------------------
% nu_rsna
if     t < tchange
    f(67-1) = nu_rsna - ( 1.822 - 2.056 / (1.358 + exp(rsna - 0.8667)) );
elseif t >= tchange
    if   strcmp(scenario,'Denerve') || strcmp(scenario,'Denerve & AT2R-')
%         f(67-1) = nu_rsna - ( SSdata(66-1) );
        f(67-1) = nu_rsna - ( 1 );
    else
        f(67-1) = nu_rsna - ( 1.822 - 2.056 / (1.358 + exp(rsna - 0.8662)) );
    end
end
% MODIFY ------------------------------------------------------------------
% C_al - rat
% if     strcmp(gender,  'male')
%     f(68-1) = C_al - ( N_al * 85      );
% elseif strcmp(gender,'female')
%     f(68-1) = C_al - ( N_al * 69.1775 );
% end
if     strcmp(gender,  'male')
    f(68-1) = C_al - ( N_al * 395 );
elseif strcmp(gender,'female')
    f(68-1) = C_al - ( N_al * 379 );
%     f(68-1) = C_al - ( N_al * 395.3 ); % male
end
% N_al
f(69-1) = N_al_p - ( 1/T_al * (N_als - N_al) );
% N_als
f(70-1) = N_als - ( xi_ksod * xi_map * xi_at );
% xi_ksod
f(71-1) = xi_ksod - ( 5 / ( 1 + exp(0.265 * (C_sod/C_K - fixed_var_pars(8))) ) ); 
% xi_map
if P_ma <= 100
    f(72-1) = xi_map - ( (1/exp(-0.0425 * 100)) * exp(-0.0425 * P_ma) );
else
    f(72-1) = xi_map - ( 1 );
end
% xi_at
f(73-1) = xi_at - ( 0.47 + 2.4 / (1 + exp(3.525 - 2.2642 * (AT1R/AT1R_eq))) );
% hatC_anp
f(74-1) = hatC_anp - ( 7.4052 - 6.554 / (1 + exp(P_ra - 3.762)) ); 
% AGT
f(75-1) = AGT_p - ( k_AGT - PRA - log(2)/h_AGT * AGT );
% nu_AT1
f(76-1) = nu_AT1 - ( (AT1R / AT1R_eq)^(-0.95) );
% R_sec
f(77-1) = R_sec - ( N_rs * nu_mdsod * nu_rsna * nu_AT1 );
% PRC
f(78-1) = PRC_p - ( R_sec - log(2)/h_renin * PRC );
% PRA
f(79-1) = PRA - ( PRC * X_PRCPRA );
% AngI
f(80-1) = AngI_p - ( PRA - (c_ACE + c_Chym + c_NEP) * AngI - log(2)/h_AngI * AngI );
% AngII
f(81-1) = AngII_p - ( k_AngII + (c_ACE + c_Chym) * AngI - (c_ACE2 + c_IIIV + c_AT1R + c_AT2R) * AngII - log(2)/h_AngII * AngII );
% AT1R
f(82-1) = AT1R_p - ( c_AT1R * AngII - log(2)/h_AT1R * AT1R );
% MODIFY ------------------------------------------------------------------
% AT2R
if   strcmp(scenario,'Denerve & AT2R-')
    alpha = 10;
    f(83-1) = AT2R_p - ( c_AT2R * AngII - log(2)/h_AT2R * AT2R - alpha * AT2R );
else
    f(83-1) = AT2R_p - ( c_AT2R * AngII - log(2)/h_AT2R * AT2R );
end
% MODIFY ------------------------------------------------------------------
% Ang17
f(84-1) = Ang17_p - ( c_NEP * AngI + c_ACE2 * AngII - log(2)/h_Ang17 * Ang17 );
% AngIV
f(85-1) = AngIV_p - ( c_IIIV * AngII - log(2)/h_AngIV * AngIV );
% R_aa
f(86-1) = R_aa - ( R_aass * beta_rsna * Sigma_tgf * Sigma_myo * Psi_AT1RAA * Psi_AT2RAA);
% R_ea
f(87-1) = R_ea - ( R_eass * Psi_AT1REA * Psi_AT2REA );
% Sigma_myo
f(88-1) = Sigma_myo - ( 0.8 + 1.2 / ( 1 + (5/1) * exp(-0.4 * (P_gh - fixed_var_pars(9))) ) );
% f(88-1) = Sigma_myo - ( 5 * (P_gh / 62 - 1) + 1 );
% Psi_AT1RAA
f(89-1) = Psi_AT1RAA - ( 0.8   + 0.2092 * (AT1R / AT1R_eq) - 0.0092 / (AT1R / AT1R_eq) );
% Psi_AT1REA
f(90-1) = Psi_AT1REA - ( 0.925 + 0.0835 * (AT1R / AT1R_eq) - 0.0085 / (AT1R / AT1R_eq) );
% Psi_AT2RAA
if     strcmp(gender,'male')
    f(91-1) = Psi_AT2RAA - ( 1 );
elseif strcmp(gender,'female')
%     f(91-1) = Psi_AT2RAA - ( 0.025 * (AT2R_eq - AT2R) + 1 );
    f(91-1) = Psi_AT2RAA - ( 0.9 + 0.1 * exp(-(AT2R/AT2R_eq - 1)) );
%     f(91-1) = Psi_AT2RAA - ( 0.7 + 0.3 * exp(-1.3 * (AT2R/AT2R_eq - 1)) );
%     f(91-1) = Psi_AT2REA - ( 0.93 + 0.07 * exp(-0.9 * (AT2R/AT2R_eq - 1)) );
%     f(91-1) = Psi_AT2RAA - ( 1 );
end
% Psi_AT2REA
if     strcmp(gender,'male')
    f(92-1) = Psi_AT2REA - ( 1 );
elseif strcmp(gender,'female')
%     f(92-1) = Psi_AT2REA - ( 0.01  * (AT2R_eq - AT2R) + 1 );
    f(92-1) = Psi_AT2REA - ( 0.9 + 0.1 * exp(-(AT2R/AT2R_eq - 1)) );
%     f(92-1) = Psi_AT2REA - ( 0.8 + 0.2 * exp(-1.4 * (AT2R/AT2R_eq - 1)) );
%     f(92-1) = Psi_AT2RAA - ( 0.85 + 0.15 * exp(-0.85 * (AT2R/AT2R_eq - 1)) );
%     f(92-1) = Psi_AT2REA - ( 1 );
end

end






























