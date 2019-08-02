% This is just a list of all the variables for reference.

% fixed_ind = [2, 10, 14, 24, 44, 49, 66, 71, 88];

% Order in which they appear in the code.
rsna          = x(1 ); rsna_p          = x_p(1 ); 
alpha_map     = x(2 ); alpha_map_p     = x_p(2 ); % fi
alpha_rap     = x(3 ); alpha_rap_p     = x_p(3 ); 
R_r           = x(4 ); R_r_p           = x_p(4 ); 
beta_rsna     = x(5 ); beta_rsna_p     = x_p(5 ); 
Phi_rb        = x(6 ); Phi_rb_p        = x_p(6 ); 
Phi_gfilt     = x(7 ); Phi_gfilt_p     = x_p(7 ); 
P_f           = x(8 ); P_f_p           = x_p(8 ); 
P_gh          = x(9 ); P_gh_p          = x_p(9 ); 
Sigma_tgf     = x(10); Sigma_tgf_p     = x_p(10); % fi % sfs
Phi_filsod    = x(11); Phi_filsod_p    = x_p(11); 
Phi_ptsodreab = x(12); Phi_ptsodreab_p = x_p(12); 
eta_ptsodreab = x(13); eta_ptsodreab_p = x_p(13); 
gamma_filsod  = x(14); gamma_filsod_p  = x_p(14); % fi % sfs
gamma_at      = x(15); gamma_at_p      = x_p(15); 
gamma_rsna    = x(16); gamma_rsna_p    = x_p(16); 
Phi_mdsod     = x(17); Phi_mdsod_p     = x_p(17); 
Phi_dtsodreab = x(18); Phi_dtsodreab_p = x_p(18); 
eta_dtsodreab = x(19); eta_dtsodreab_p = x_p(19); 
psi_al        = x(20); psi_al_p        = x_p(20); 
Phi_dtsod     = x(21); Phi_dtsod_p     = x_p(21); 
Phi_cdsodreab = x(22); Phi_cdsodreab_p = x_p(22); 
eta_cdsodreab = x(23); eta_cdsodreab_p = x_p(23); 
lambda_dt     = x(24); lambda_dt_p     = x_p(24); % fi % sfs
lambda_anp    = x(25); lambda_anp_p    = x_p(25); 
lambda_al     = x(26); lambda_al_p     = x_p(26); 
Phi_usod      = x(27); Phi_usod_p      = x_p(27); 
Phi_win       = x(28); Phi_win_p       = x_p(28); 
V_ecf         = x(29); V_ecf_p         = x_p(29); 
V_b           = x(30); V_b_p           = x_p(30); % sfv 
P_mf          = x(31); P_mf_p          = x_p(31); % sfv 
Phi_vr        = x(32); Phi_vr_p        = x_p(32); 
Phi_co        = x(33); Phi_co_p        = x_p(33); 
P_ra          = x(34); P_ra_p          = x_p(34); % sfr 
vas           = x(35); vas_p           = x_p(35); 
vas_f         = x(36); vas_f_p         = x_p(36); 
vas_d         = x(37); vas_d_p         = x_p(37); 
R_a           = x(38); R_a_p           = x_p(38); 
R_ba          = x(39); R_ba_p          = x_p(39); 
R_vr          = x(40); R_vr_p          = x_p(40); 
R_tp          = x(41); R_tp_p          = x_p(41); 
P_ma          = x(42); P_ma_p          = x_p(42); 
epsilon_aum   = x(43); epsilon_aum_p   = x_p(43); 
a_auto        = x(44); a_auto_p        = x_p(44); % fi
a_chemo       = x(45); a_chemo_p       = x_p(45); 
a_baro        = x(46); a_baro_p        = x_p(46); 
C_adh         = x(47); C_adh_p         = x_p(47); 
N_adh         = x(48); N_adh_p         = x_p(48); 
N_adhs        = x(49); N_adhs_p        = x_p(49); % fi
delta_ra      = x(50); delta_ra_p      = x_p(50); 
Phi_ptwreab   = x(51); Phi_ptwreab_p   = x_p(51);                           %#ok<*NASGU>
eta_ptwreab   = x(52); eta_ptwreab_p   = x_p(52); 
mu_ptsodreab  = x(53); mu_ptsodreab_p  = x_p(53); 
Phi_mdu       = x(54); Phi_mdu_p       = x_p(54); 
Phi_dtwreab   = x(55); Phi_dtwreab_p   = x_p(55); 
eta_dtwreab   = x(56); eta_dtwreab_p   = x_p(56); 
mu_dtsodreab  = x(57); mu_dtsodreab_p  = x_p(57); 
Phi_dtu       = x(58); Phi_dtu_p       = x_p(58); 
Phi_cdwreab   = x(59); Phi_cdwreab_p   = x_p(59); 
eta_cdwreab   = x(60); eta_cdwreab_p   = x_p(60); 
mu_cdsodreab  = x(61); mu_cdsodreab_p  = x_p(61); 
mu_adh        = x(62); mu_adh_p        = x_p(62); 
Phi_u         = x(63); Phi_u_p         = x_p(63); 
M_sod         = x(64); M_sod_p         = x_p(64); 
C_sod         = x(65); C_sod_p         = x_p(65); 
nu_mdsod      = x(66); nu_mdsod_p      = x_p(66); % fi % sfs
nu_rsna       = x(67); nu_rsna_p       = x_p(67); 
C_al          = x(68); C_al_p          = x_p(68); 
N_al          = x(69); N_al_p          = x_p(69); 
N_als         = x(70); N_als_p         = x_p(70); 
xi_ksod       = x(71); xi_ksod_p       = x_p(71); % fi
xi_map        = x(72); xi_map_p        = x_p(72); 
xi_at         = x(73); xi_at_p         = x_p(73); 
hatC_anp      = x(74); hatC_anp_p      = x_p(74); 
AGT           = x(75); AGT_p           = x_p(75); 
nu_AT1        = x(76); nu_AT1_p        = x_p(76); 
R_sec         = x(77); R_sec_p         = x_p(77); 
PRC           = x(78); PRC_p           = x_p(78); 
PRA           = x(79); PRA_p           = x_p(79); 
AngI          = x(80); AngI_p          = x_p(80); 
AngII         = x(81); AngII_p         = x_p(81); 
AT1R          = x(82); AT1R_p          = x_p(82); 
AT2R          = x(83); AT2R_p          = x_p(83); 
Ang17         = x(84); Ang17_p         = x_p(84); 
AngIV         = x(85); AngIV_p         = x_p(85); 
R_aa          = x(86); R_aa_p          = x_p(86); 
R_ea          = x(87); R_ea_p          = x_p(87); 
Sigma_myo     = x(88); Sigma_myo_p     = x_p(88); % fi
Psi_AT1RAA    = x(89); Psi_AT1RAA_p    = x_p(89); 
Psi_AT1REA    = x(90); Psi_AT1REA_p    = x_p(90); 
Psi_AT2RAA    = x(91); Psi_AT2RAA_p    = x_p(91); 
Psi_AT2REA    = x(92); Psi_AT2REA_p    = x_p(92); 



% fixed_ind = [2, 10, 14, 24, 44, 49, 66, 71, 88];

% Order in which they appear in the manuscript.
% Cardiovascular function
P_mf          = x(31); P_mf_p          = x_p(31); % sfv 
P_ra          = x(34); P_ra_p          = x_p(34); % sfr 
P_ma          = x(42); P_ma_p          = x_p(42); 
Phi_vr        = x(32); Phi_vr_p        = x_p(32); 
Phi_co        = x(33); Phi_co_p        = x_p(33); 
vas           = x(35); vas_p           = x_p(35); 
vas_f         = x(36); vas_f_p         = x_p(36); 
vas_d         = x(37); vas_d_p         = x_p(37); 
R_a           = x(38); R_a_p           = x_p(38); 
R_ba          = x(39); R_ba_p          = x_p(39); 
R_vr          = x(40); R_vr_p          = x_p(40); 
R_tp          = x(41); R_tp_p          = x_p(41); 
epsilon_aum   = x(43); epsilon_aum_p   = x_p(43); 
a_chemo       = x(45); a_chemo_p       = x_p(45); 
a_auto        = x(44); a_auto_p        = x_p(44); % fi
a_baro        = x(46); a_baro_p        = x_p(46); 
% Renal hemodynamics
Phi_rb        = x(6 ); Phi_rb_p        = x_p(6 ); 
Phi_gfilt     = x(7 ); Phi_gfilt_p     = x_p(7 ); 
P_f           = x(8 ); P_f_p           = x_p(8 ); 
P_gh          = x(9 ); P_gh_p          = x_p(9 ); 
R_r           = x(4 ); R_r_p           = x_p(4 ); 
R_aa          = x(86); R_aa_p          = x_p(86); 
R_ea          = x(87); R_ea_p          = x_p(87); 
beta_rsna     = x(5 ); beta_rsna_p     = x_p(5 ); 
Sigma_tgf     = x(10); Sigma_tgf_p     = x_p(10); % fi % sfs
Sigma_myo     = x(88); Sigma_myo_p     = x_p(88); % fi
Psi_AT1RAA    = x(89); Psi_AT1RAA_p    = x_p(89); 
Psi_AT1REA    = x(90); Psi_AT1REA_p    = x_p(90); 
Psi_AT2RAA    = x(91); Psi_AT2RAA_p    = x_p(91); 
Psi_AT2REA    = x(92); Psi_AT2REA_p    = x_p(92); 
% Renal function
Phi_filsod    = x(11); Phi_filsod_p    = x_p(11); 
Phi_ptsodreab = x(12); Phi_ptsodreab_p = x_p(12); 
Phi_mdsod     = x(17); Phi_mdsod_p     = x_p(17); 
Phi_dtsodreab = x(18); Phi_dtsodreab_p = x_p(18); 
Phi_dtsod     = x(21); Phi_dtsod_p     = x_p(21); 
Phi_cdsodreab = x(22); Phi_cdsodreab_p = x_p(22); 
Phi_usod      = x(27); Phi_usod_p      = x_p(27); 
eta_ptsodreab = x(13); eta_ptsodreab_p = x_p(13); 
eta_dtsodreab = x(19); eta_dtsodreab_p = x_p(19); 
eta_cdsodreab = x(23); eta_cdsodreab_p = x_p(23); 
gamma_filsod  = x(14); gamma_filsod_p  = x_p(14); % fi % sfs
gamma_at      = x(15); gamma_at_p      = x_p(15); 
gamma_rsna    = x(16); gamma_rsna_p    = x_p(16); 
psi_al        = x(20); psi_al_p        = x_p(20); 
lambda_dt     = x(24); lambda_dt_p     = x_p(24); % fi % sfs
lambda_anp    = x(25); lambda_anp_p    = x_p(25); 
lambda_al     = x(26); lambda_al_p     = x_p(26); 
Phi_ptwreab   = x(51); Phi_ptwreab_p   = x_p(51); 
Phi_mdu       = x(54); Phi_mdu_p       = x_p(54); 
Phi_dtwreab   = x(55); Phi_dtwreab_p   = x_p(55); 
Phi_dtu       = x(58); Phi_dtu_p       = x_p(58); 
Phi_cdwreab   = x(59); Phi_cdwreab_p   = x_p(59); 
Phi_u         = x(63); Phi_u_p         = x_p(63); 
eta_ptwreab   = x(52); eta_ptwreab_p   = x_p(52); 
eta_dtwreab   = x(56); eta_dtwreab_p   = x_p(56); 
eta_cdwreab   = x(60); eta_cdwreab_p   = x_p(60); 
mu_ptsodreab  = x(53); mu_ptsodreab_p  = x_p(53); 
mu_dtsodreab  = x(57); mu_dtsodreab_p  = x_p(57); 
mu_cdsodreab  = x(61); mu_cdsodreab_p  = x_p(61); 
mu_adh        = x(62); mu_adh_p        = x_p(62); 
% Renin-angiotensin-aldosterone system
R_sec         = x(77); R_sec_p         = x_p(77); 
PRC           = x(78); PRC_p           = x_p(78); 
PRA           = x(79); PRA_p           = x_p(79); 
AGT           = x(75); AGT_p           = x_p(75); 
AngI          = x(80); AngI_p          = x_p(80); 
AngII         = x(81); AngII_p         = x_p(81); 
Ang17         = x(84); Ang17_p         = x_p(84); 
AngIV         = x(85); AngIV_p         = x_p(85); 
AT1R          = x(82); AT1R_p          = x_p(82); 
AT2R          = x(83); AT2R_p          = x_p(83); 
nu_mdsod      = x(66); nu_mdsod_p      = x_p(66); % fi % sfs
nu_rsna       = x(67); nu_rsna_p       = x_p(67); 
nu_AT1        = x(76); nu_AT1_p        = x_p(76); 
N_als         = x(70); N_als_p         = x_p(70); 
N_al          = x(69); N_al_p          = x_p(69); 
C_al          = x(68); C_al_p          = x_p(68); 
xi_ksod       = x(71); xi_ksod_p       = x_p(71); % fi
xi_map        = x(72); xi_map_p        = x_p(72); 
xi_at         = x(73); xi_at_p         = x_p(73); 
% Miscellaneous
rsna          = x(1 ); rsna_p          = x_p(1 ); 
alpha_map     = x(2 ); alpha_map_p     = x_p(2 ); % fi
alpha_rap     = x(3 ); alpha_rap_p     = x_p(3 ); 
Phi_win       = x(28); Phi_win_p       = x_p(28); 
V_ecf         = x(29); V_ecf_p         = x_p(29); 
V_b           = x(30); V_b_p           = x_p(30); % sfv 
N_adhs        = x(49); N_adhs_p        = x_p(49); % fi
N_adh         = x(48); N_adh_p         = x_p(48); 
C_adh         = x(47); C_adh_p         = x_p(47); 
delta_ra      = x(50); delta_ra_p      = x_p(50); 
M_sod         = x(64); M_sod_p         = x_p(64); 
C_sod         = x(65); C_sod_p         = x_p(65); 
hatC_anp      = x(74); hatC_anp_p      = x_p(74); 




























