% This is a model of long-term model blood pressure regulation.
% It is adopted with modifications from Karaaslan 2005 and Leete 2018.

% This function file is to solve for the steady state solution in different scenarios.

% Input
% t        - time
% x        - variables
% x_p      - variable derivatives
% pars     - parameters
% tchange  - time at which to change RPP in simulation
% varargin - cell containing different scenarios and the corresponding parameters 
%            it is a variable length input argument

% Output
% f        - left hand side of f(t,x(t),x'(t);theta) = 0.

function f = bp_reg_mod(t,x,x_p,pars,tchange,varargin)

%% Retrieve species and sex identifier. 

spe_par = pars(44);
sex_par = pars(45);
if     spe_par == 1
    species = 'human';
elseif spe_par == 0
    species = 'rat';
end
if     sex_par == 1
    sex = 'male';
elseif sex_par == 0
    sex = 'female';
end

%% Retrieve species specific things.

if     strcmp(species, 'human')
    % Scaling factors
    % Rat value = Human value x SF
    % Note: This includes conversion of units.
    SF_S =                    1; % sodium flow
    SF_U =                    1; % urine flow
    SF_R =                    1; % resistance
    SF_V =                    1; % volume
elseif strcmp(species, 'rat')
    % Steady state variable values
    SSdata_input     = pars(end-92:end);
    pars(end-92:end) = '';
    
    % Fixed variable parameters
    fixed_var_pars  = pars(end-8:end);
    pars(end-8:end) = '';

    % Physiological variables which determine scaling factors.
    Phi_usod_new = pars(end-3)      ; % Munger 1988, Karaaslan 2005
    Phi_u_new    = pars(end-2)      ; % Munger 1988, Layton 2016
    R_r_new      = pars(end-1)      ; % Munger 1988
    W_b          = pars(end  )      ; % Munger 1988
    V_b_new      = 0.06 * W_b + 0.77; % Lee 1985
    % Scaling factors
    % Rat value = Human value x SF
    % Note: This includes conversion of units.
    SF_S = Phi_usod_new / 0.126; % sodium flow
    SF_U = Phi_u_new    / 0.001; % urine flow
    SF_R = R_r_new      / 83.3 ; % resistance
    SF_V = V_b_new      / 5    ; % volume
end

%% Default parameter inputs for changes in simulation.

kappa_AngII = 0; % Ang II infusion rate fmol/ml/min
kappa_ACEi  = 0; % ACE inhibitor %
kappa_ARB   = 0; % Angiotensin receptor blocker %

RPP_ind = false; % Boolean for controlling renal perfusion pressure
RPP_per = 0    ; % Renal perfusion pressure mmHg
denerve = false; % Boolean for unilateral renal denervation by fixing rsna = 1, which is baseline.
fix_win = false; % Boolean for fixing water intake.

no_TGF  = false; % Boolean for canceling the tubuloglomerular feedback.
no_Myo  = false; % Boolean for canceling the myogenic response.
lin_Myo = false; % Boolean for having linear myogenic response.

m_RSNA = false; % Boolean for having male RSNA in the female model.
m_AT2R = false; % Boolean for having male effect of AT2R in the female model.

%% Read and assign optional parameters.

% The odd inputs of varargin are strings for each scenario. The
% corresponding even inputs are the values for the effect parameters to
% modify something.

varargin = varargin{:};
for i = 1:2:length(varargin)
    
    if     strcmp(varargin{i},'AngII_')
        kappa_AngII = varargin{i + 1};
    elseif strcmp(varargin{i},'ACEi_')
        kappa_ACEi  = varargin{i + 1};
    elseif strcmp(varargin{i},'ARB_')
        kappa_ARB  = varargin{i + 1};
        
    elseif strcmp(varargin{i},'RPP_')
        RPP_ind = true;
        RPP_per = varargin{i + 1};
    elseif strcmp(varargin{i},'Denerve_')
        denerve = varargin{i + 1};
    elseif strcmp(varargin{i},'Fixed_WatIn_')
        fix_win = varargin{i + 1};
        
    elseif strcmp(varargin{i},'no_TGF_')
        no_TGF  = varargin{i + 1};
    elseif strcmp(varargin{i},'no_Myo_')
        no_Myo  = varargin{i + 1};
    elseif strcmp(varargin{i},'lin_Myo_')
        lin_Myo = varargin{i + 1};
        
    elseif strcmp(varargin{i},'m_RSNA_') || strcmp(varargin{i},'m_RSNA_m_Reab_')
        m_RSNA = varargin{i + 1};
    elseif strcmp(varargin{i},'m_AT2R_')
        m_AT2R = varargin{i + 1};
    end
end

if RPP_ind
    if     t <  tchange
        RPP = SSdata_input(42);
    elseif t >= tchange 
        RPP = RPP_per * tanh(5 * (t-tchange)) + SSdata_input(42);
    end
end

if     t <  tchange
    kappa_AngII = 0; 
    kappa_ACEi  = 0; 
    kappa_ARB   = 0;
elseif t >= tchange 
    kappa_AngII = kappa_AngII; 
    kappa_ACEi  = kappa_ACEi; 
    kappa_ARB   = kappa_ARB;
end

%% Formulate human to rat scaling factors.

% Rat value = Human value x SF
% Note: This includes conversion of units.
if     strcmp(species, 'human')
    SF_S =                    1; % sodium flow
    SF_U =                    1; % urine flow
    SF_R =                    1; % resistance
    SF_V =                    1; % volume
elseif strcmp(species, 'rat')
    % Physiological variables which determine scaling factors.
    Phi_usod_new = pars(end-3)      ; % Munger 1988, Karaaslan 2005
    Phi_u_new    = pars(end-2)      ; % Munger 1988, Layton 2016
    R_r_new      = pars(end-1)      ; % Munger 1988
    W_b          = pars(end  )      ; % Munger 1988
    V_b_new      = 0.06 * W_b + 0.77; % Lee 1985
    
    SF_S = Phi_usod_new / 0.126; % sodium flow
    SF_U = Phi_u_new    / 0.001; % urine flow
    SF_R = R_r_new      / 83.3 ; % resistance
    SF_V = V_b_new      / 5    ; % volume
end

%% Retrieve parameters by name.

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
K_bar            = pars(14) * SF_R;
R_bv             = pars(15) * SF_R;
N_adhs_eq        = pars(16);
T_adh            = pars(17);
Phi_sodin_const  = pars(18);
N_als_eq         = pars(19);
C_K              = pars(20);
T_al             = pars(21);
N_rs             = pars(22);
X_PRCPRA         = pars(23);
h_renin          = pars(24);
h_AGT            = pars(25);
h_AngI           = pars(26);
h_AngII          = pars(27);
h_Ang17          = pars(28);
h_AngIV          = pars(29);
h_AT1R           = pars(30);
h_AT2R           = pars(31);
k_AGT            = pars(32);
c_ACE            = pars(33) * (1-kappa_ACEi);
c_Chym           = pars(34);
c_NEP            = pars(35);
c_ACE2           = pars(36);
c_IIIV           = pars(37);
c_AT1R           = pars(38) * (1-kappa_ARB);
c_AT2R           = pars(39);
AT1R_eq          = pars(40);
AT2R_eq          = pars(41);
Psi_AT2RAA_eq    = pars(42);
Psi_AT2REA_eq    = pars(43);

%% Retrieve variables by name.

rsna          = x(1 ); rsna_p          = x_p(1 ); 
alpha_map     = x(2 ); alpha_map_p     = x_p(2 ); %v
alpha_rap     = x(3 ); alpha_rap_p     = x_p(3 ); 
R_r           = x(4 ); R_r_p           = x_p(4 ); 
beta_rsna     = x(5 ); beta_rsna_p     = x_p(5 ); 
Phi_rb        = x(6 ); Phi_rb_p        = x_p(6 ); 
Phi_gfilt     = x(7 ); Phi_gfilt_p     = x_p(7 ); 
P_f           = x(8 ); P_f_p           = x_p(8 ); 
P_gh          = x(9 ); P_gh_p          = x_p(9 ); 
Sigma_tgf     = x(10); Sigma_tgf_p     = x_p(10); %v
Phi_filsod    = x(11); Phi_filsod_p    = x_p(11); 
Phi_ptsodreab = x(12); Phi_ptsodreab_p = x_p(12); 
eta_ptsodreab = x(13); eta_ptsodreab_p = x_p(13); 
gamma_filsod  = x(14); gamma_filsod_p  = x_p(14); %v
gamma_at      = x(15); gamma_at_p      = x_p(15); 
gamma_rsna    = x(16); gamma_rsna_p    = x_p(16); 
Phi_mdsod     = x(17); Phi_mdsod_p     = x_p(17); 
Phi_dtsodreab = x(18); Phi_dtsodreab_p = x_p(18); 
eta_dtsodreab = x(19); eta_dtsodreab_p = x_p(19); 
psi_al        = x(20); psi_al_p        = x_p(20); 
Phi_dtsod     = x(21); Phi_dtsod_p     = x_p(21); 
Phi_cdsodreab = x(22); Phi_cdsodreab_p = x_p(22); 
eta_cdsodreab = x(23); eta_cdsodreab_p = x_p(23); 
lambda_dt     = x(24); lambda_dt_p     = x_p(24); %v
lambda_anp    = x(25); lambda_anp_p    = x_p(25); 
lambda_al     = x(26); lambda_al_p     = x_p(26); 
Phi_usod      = x(27); Phi_usod_p      = x_p(27); 
Phi_sodin     = x(28); Phi_sodin_p     = x_p(28); 
V_ecf         = x(29); V_ecf_p         = x_p(29); 
V_b           = x(30); V_b_p           = x_p(30); 
P_mf          = x(31); P_mf_p          = x_p(31); 
Phi_vr        = x(32); Phi_vr_p        = x_p(32); 
Phi_co        = x(33); Phi_co_p        = x_p(33); 
P_ra          = x(34); P_ra_p          = x_p(34); 
vas           = x(35); vas_p           = x_p(35); 
vas_f         = x(36); vas_f_p         = x_p(36); 
vas_d         = x(37); vas_d_p         = x_p(37); 
R_a           = x(38); R_a_p           = x_p(38); 
R_ba          = x(39); R_ba_p          = x_p(39); 
R_vr          = x(40); R_vr_p          = x_p(40); 
R_tp          = x(41); R_tp_p          = x_p(41); 
P_ma          = x(42); P_ma_p          = x_p(42); 
epsilon_aum   = x(43); epsilon_aum_p   = x_p(43); 
a_auto        = x(44); a_auto_p        = x_p(44); %v
a_chemo       = x(45); a_chemo_p       = x_p(45); 
a_baro        = x(46); a_baro_p        = x_p(46); 
C_adh         = x(47); C_adh_p         = x_p(47); 
N_adh         = x(48); N_adh_p         = x_p(48); 
N_adhs        = x(49); N_adhs_p        = x_p(49); %v
delta_ra      = x(50); delta_ra_p      = x_p(50); 

M_sod         = x(51); M_sod_p         = x_p(51); 
C_sod         = x(52); C_sod_p         = x_p(52); 
nu_mdsod      = x(53); nu_mdsod_p      = x_p(53); % _v
nu_rsna       = x(54); nu_rsna_p       = x_p(54); 
C_al          = x(55); C_al_p          = x_p(55); 
N_al          = x(56); N_al_p          = x_p(56); 
N_als         = x(57); N_als_p         = x_p(57); 
xi_ksod       = x(58); xi_ksod_p       = x_p(58); % _v
xi_map        = x(59); xi_map_p        = x_p(59); 
xi_at         = x(60); xi_at_p         = x_p(60); 
hatC_anp      = x(61); hatC_anp_p      = x_p(61); 
AGT           = x(62); AGT_p           = x_p(62); 
nu_AT1        = x(63); nu_AT1_p        = x_p(63); 
R_sec         = x(64); R_sec_p         = x_p(64); 
PRC           = x(65); PRC_p           = x_p(65); 
PRA           = x(66); PRA_p           = x_p(66); 
AngI          = x(67); AngI_p          = x_p(67); 
AngII         = x(68); AngII_p         = x_p(68); 
AT1R          = x(69); AT1R_p          = x_p(69); 
AT2R          = x(70); AT2R_p          = x_p(70); 
Ang17         = x(71); Ang17_p         = x_p(71); 
AngIV         = x(72); AngIV_p         = x_p(72); 
R_aa          = x(73); R_aa_p          = x_p(73); 
R_ea          = x(74); R_ea_p          = x_p(74); 
Sigma_myo     = x(75); Sigma_myo_p     = x_p(75); % _v
Psi_AT1RAA    = x(76); Psi_AT1RAA_p    = x_p(76); 
Psi_AT1REA    = x(77); Psi_AT1REA_p    = x_p(77); 
Psi_AT2RAA    = x(78); Psi_AT2RAA_p    = x_p(78); 
Psi_AT2REA    = x(79); Psi_AT2REA_p    = x_p(79); 

Phi_ptwreab   = x(80); Phi_ptwreab_p   = x_p(80); 
eta_ptwreab   = x(81); eta_ptwreab_p   = x_p(81); 
mu_ptsodreab  = x(82); mu_ptsodreab_p  = x_p(82); 
Phi_mdu       = x(83); Phi_mdu_p       = x_p(83); 
Phi_dtwreab   = x(84); Phi_dtwreab_p   = x_p(84); 
eta_dtwreab   = x(85); eta_dtwreab_p   = x_p(85); 
mu_dtsodreab  = x(86); mu_dtsodreab_p  = x_p(86); 
Phi_dtu       = x(87); Phi_dtu_p       = x_p(87); 
Phi_cdwreab   = x(88); Phi_cdwreab_p   = x_p(88); 
eta_cdwreab   = x(89); eta_cdwreab_p   = x_p(89); 
mu_cdsodreab  = x(90); mu_cdsodreab_p  = x_p(90); 
mu_adh        = x(91); mu_adh_p        = x_p(91); 
Phi_u         = x(92); Phi_u_p         = x_p(92); 
Phi_win       = x(93); Phi_win_p       = x_p(93); 

%% Differential algebraic equation system f(t,x,x') = 0.

f = zeros(length(x),1);

% rsna
rsna0 = N_rsna * alpha_map * alpha_rap;
if     strcmp(sex,'male')
        f(1 ) = rsna - rsna0;
elseif strcmp(sex,'female')
    if   m_RSNA
        f(1 ) = rsna - rsna0;
    else
        f(1 ) = rsna - rsna0^(1/rsna0);
    end
end
% alpha_map
f(2 ) = alpha_map - ( 0.5 + 1 / (1 + exp((P_ma - fixed_var_pars(1)) / 15)) );
% alpha_rap
f(3 ) = alpha_rap - ( 1 - 0.008 * P_ra );
% R_r
f(4 ) = R_r - ( R_aa + R_ea );
% beta_rsna
if   denerve
    f(5 ) = beta_rsna - ( 1 );
else
    f(5 ) = beta_rsna - ( 2 / (1 + exp(-3.16 * (rsna - 1))) );
end
% Phi_rb
if     RPP_ind
    f(6 ) = Phi_rb - ( RPP / R_r );
else
    f(6 ) = Phi_rb - ( P_ma / R_r );
end
% Phi_gfilt
f(7 ) = Phi_gfilt - ( P_f * C_gcf );
% P_f
f(8 ) = P_f - ( P_gh - (P_B + P_go) );
% P_gh
if     RPP_ind
    f(9 ) = P_gh - ( RPP - Phi_rb * R_aa );
else
    f(9 ) = P_gh - ( P_ma - Phi_rb * R_aa );
end
% Sigma_tgf
if   no_TGF
    f(10) = Sigma_tgf - ( 1 );
else
    f(10) = Sigma_tgf - ( 0.3408 + 3.449 / (3.88 + exp((Phi_mdsod - fixed_var_pars(2)) / (-0.9617 * SF_S) )) );
end
% Phi_filsod
f(11) = Phi_filsod - ( Phi_gfilt * C_sod );
% Phi_ptsodreab
f(12) = Phi_ptsodreab - ( Phi_filsod * eta_ptsodreab );
% eta_ptsodreab
f(13) = eta_ptsodreab - ( eta_ptsodreab_eq * gamma_filsod * gamma_at * gamma_rsna );
% gamma_filsod
f(14) = gamma_filsod - ( 0.8 + 0.3 / (1 + exp((Phi_filsod - fixed_var_pars(3))/(138 * SF_S) )) );
% gamma_at
gammaat_c = 0.8017; gammaat_d = 0.92; gammaat_a = 0.136;
gammaat_b = -1/(1-gammaat_c) * log(gammaat_a/(1-gammaat_d) - 1);
f(15) = gamma_at - ( gammaat_d + gammaat_a / (1 + exp(-gammaat_b * ((AT1R/AT1R_eq) - gammaat_c)) ) );
% gamma_rsna
if   denerve
    f(16) = gamma_rsna - ( 1 );
else
    f(16) = gamma_rsna - ( 0.72 + 0.56 / (1 + exp((1 - rsna) / 2.18)) );
end
% Phi_mdsod
f(17) = Phi_mdsod - ( Phi_filsod - Phi_ptsodreab );
% Phi_dtsodreab
f(18) = Phi_dtsodreab - ( Phi_mdsod * eta_dtsodreab );
% eta_dtsodreab
f(19) = eta_dtsodreab - ( eta_dtsodreab_eq * psi_al );
% psi_al
psial_b = 0.1; psial_d = 1.05 / psial_b; psial_a = (1 + psial_b) * psial_d;
psial_c = -1/387 * log((psial_a / (1 + psial_d) - 1) / psial_b);
f(20) = psi_al - ( psial_a / (1 + psial_b * exp(-psial_c * C_al)) - psial_d );
% Phi_dtsod
f(21) = Phi_dtsod - ( Phi_mdsod - Phi_dtsodreab );
% Phi_cdsodreab
f(22) = Phi_cdsodreab - ( Phi_dtsod * eta_cdsodreab );
% eta_cdsodreab
f(23) = eta_cdsodreab - ( eta_cdsodreab_eq * lambda_dt * lambda_anp );
% lambda_dt
if     strcmp(sex,  'male')
    lambdadt_b = 0.275; lambdadt_c = 2.3140;
elseif strcmp(sex,'female')
    lambdadt_b = 1/eta_cdsodreab_eq - 0.8; lambdadt_c = 2.86;
end
f(24) = lambda_dt - ( 0.8 + lambdadt_b/ (1 + exp(lambdadt_c/SF_S * (Phi_dtsod - fixed_var_pars(4))) ) );
% lambda_anp
f(25) = lambda_anp - ( -0.1 * hatC_anp + 1.1 );
% lambda_al
f(26) = lambda_al - ( 1/(387^0.06) * C_al^0.06 );
% Phi_usod
f(27) = Phi_usod - ( Phi_dtsod - Phi_cdsodreab );

% Phi_sodin
f(28) = Phi_sodin - ( Phi_sodin_const );

% V_ecf
f(29) = V_ecf_p - ( Phi_win - Phi_u );
% V_b
f(30) = V_b - ( SF_V*( 4.5479 + 2.4312 / (1 + exp(-(V_ecf - 18.1128*SF_V) * (0.4744/SF_V) )) ) );
% P_mf
pmfpmf = (7.4360/SF_V);
f(31) = P_mf - ( ( pmfpmf * V_b - 30.18) * epsilon_aum );
% Phi_vr
f(32) = Phi_vr - ( (P_mf - P_ra) / R_vr );
% Phi_co
f(33) = Phi_co - ( Phi_vr );
% P_ra
prapra = 0.2787 * exp(SSdata_input(33) * 0.2281 * SF_R);
prapra = prapra + 0.01 * prapra;
f(34) = P_ra - ( max( 0, 0.2787 * exp(Phi_co * 0.2281 * SF_R) - prapra ) );
% vas
f(35) = vas_p - ( 1 / 1000 * (vas_f - vas_d) );
% vas_f
vvv = -1/SSdata_input(33) * log(1/11.312);
f(36) = vas_f - ( (11.312 * exp(-Phi_co * vvv)) / 100 );
% vas_d
f(37) = vas_d - ( vas * K_vd );
% R_a
f(38) = R_a - ( R_ba * epsilon_aum );
% R_ba
f(39) = R_ba - ( K_bar / vas );
% R_vr
f(40) = R_vr - ( (8 * R_bv + R_a) / 31 );
% R_tp
f(41) = R_tp - ( R_a + R_bv );
% P_ma
f(42) = P_ma - ( Phi_co * R_tp );
% epsilon_aum
f(43) = epsilon_aum - ( 4/5 * (a_chemo + a_baro) );
% a_auto
f(44) = a_auto - ( 3.0042 * exp(-fixed_var_pars(5) * P_ma) );
% a_chemo
f(45) = a_chemo - ( 1/4 * a_auto );
% a_baro
f(46) = a_baro_p - ( 3/4 * (a_auto_p - 0.0000667 * (a_baro - 1)) );
% C_adh
f(47) = C_adh - ( 4 * N_adh );
% N_adh
f(48) = N_adh_p - ( 1/T_adh * (N_adhs - N_adh) );
% N_adhs
f(49) = N_adhs - ( N_adhs_eq * (max( 0, C_sod - fixed_var_pars(6)) + max( 0, epsilon_aum - 1 ) - delta_ra) / 3 );
% delta_ra
f(50) = delta_ra_p - ( 0.2 * P_ra_p - 0.0007 * delta_ra );

% M_sod
f(51) = M_sod_p - ( Phi_sodin - Phi_usod );
% C_sod
f(52) = C_sod - ( M_sod / V_ecf );
% nu_mdsod
f(53) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - fixed_var_pars(7)) / (0.6056 * SF_S) )) );
% nu_rsna
if   denerve
    f(54) = nu_rsna - ( 1 );
else
    f(54) = nu_rsna - ( 1.822 - 2.056 / (1.358 + exp(rsna - 0.8662)) );
end
% C_al
f(55) = C_al - ( N_al * 387 );
% N_al
f(56) = N_al_p - ( 1/T_al * (N_als - N_al) );
% N_als
f(57) = N_als - ( N_als_eq * xi_ksod * xi_map * xi_at );
% xi_ksod
f(58) = xi_ksod - ( 5 / ( 1 + exp(0.265 * (C_sod/C_K - fixed_var_pars(8))) ) ); 
% xi_map
if P_ma <= 100
    f(59) = xi_map - ( (1/exp(-0.0425 * 100)) * exp(-0.0425 * P_ma) );
else
    f(59) = xi_map - ( 1 );
end
% xi_at
xiat_a = 0.1; xiat_b = 2.9; xiat_c = 2.0;
xiat_d = 1 + 1/xiat_c * log(xiat_b/(1-xiat_a) - 1);
f(60) = xi_at - ( xiat_a + xiat_b / (1 + exp(-xiat_c * ((AT1R/AT1R_eq) - xiat_d)) ) );
% hatC_anp
f(61) = hatC_anp - ( 7.4052 - 6.554 / (1 + exp(P_ra - 3.762)) ); 
% AGT
f(62) = AGT_p - ( k_AGT - PRA - log(2)/h_AGT * AGT );
% nu_AT1
f(63) = nu_AT1 - ( (AT1R / AT1R_eq)^(-0.95) );
% R_sec
f(64) = R_sec - ( N_rs * nu_mdsod * nu_rsna * nu_AT1 );
% PRC
f(65) = PRC_p - ( R_sec - log(2)/h_renin * PRC );
% PRA
f(66) = PRA - ( PRC * X_PRCPRA );
% AngI
f(67) = AngI_p - ( PRA - (c_ACE + c_Chym + c_NEP) * AngI - log(2)/h_AngI * AngI );
% AngII
f(68) = AngII_p - ( kappa_AngII + (c_ACE + c_Chym) * AngI - (c_ACE2 + c_IIIV + c_AT1R + c_AT2R) * AngII - log(2)/h_AngII * AngII );
% AT1R
f(69) = AT1R_p - ( c_AT1R * AngII - log(2)/h_AT1R * AT1R );
% AT2R
f(70) = AT2R_p - ( c_AT2R * AngII - log(2)/h_AT2R * AT2R );
% Ang17
f(71) = Ang17_p - ( c_NEP * AngI + c_ACE2 * AngII - log(2)/h_Ang17 * Ang17 );
% AngIV
f(72) = AngIV_p - ( c_IIIV * AngII - log(2)/h_AngIV * AngIV );
% R_aa
f(73) = R_aa - ( R_aass * beta_rsna * Sigma_tgf * Sigma_myo * Psi_AT1RAA * Psi_AT2RAA);
% R_ea
f(74) = R_ea - ( R_eass * Psi_AT1REA * Psi_AT2REA );
% Sigma_myo
sigmamyo_a = 0.75; sigmamyo_b = 1.2; sigmamyo_d = 0.6;
sigmamyo_c = sigmamyo_b / (1-sigmamyo_a) - 1;
if     no_Myo
    f(75) = Sigma_myo - ( 1 );
elseif lin_Myo
    f(75) = Sigma_myo - ( 5 * (P_gh / fixed_var_pars(9) - 1) + 1 );
else
    f(75) = Sigma_myo - ( sigmamyo_a + sigmamyo_b / ( 1 + sigmamyo_c * exp(-sigmamyo_d * (P_gh - fixed_var_pars(9))) ) );
end
% Psi_AT1RAA
f(76) = Psi_AT1RAA - ( 0.8   + 0.2092 * (AT1R / AT1R_eq) - 0.0092 / (AT1R / AT1R_eq) );
% Psi_AT1REA
f(77) = Psi_AT1REA - ( 0.925 + 0.0835 * (AT1R / AT1R_eq) - 0.0085 / (AT1R / AT1R_eq) );
% Psi_AT2RAA
if     strcmp(sex,'male')
        f(78) = Psi_AT2RAA - ( 1 );
elseif strcmp(sex,'female')
    if   m_AT2R
        f(78) = Psi_AT2RAA - ( 1 );
    else
        f(78) = Psi_AT2RAA - ( Psi_AT2RAA_eq * (0.9 + 0.1 * exp(-(AT2R/AT2R_eq - 1))) );
    end
end
% Psi_AT2REA
if     strcmp(sex,'male')
        f(79) = Psi_AT2REA - ( 1 );
elseif strcmp(sex,'female')
    if   m_AT2R
        f(79) = Psi_AT2REA - ( 1 );
    else
        f(79) = Psi_AT2REA - ( Psi_AT2REA_eq * (0.9 + 0.1 * exp(-(AT2R/AT2R_eq - 1))) );
    end    
end

% Phi_ptwreab
f(80) = Phi_ptwreab - ( Phi_gfilt * eta_ptwreab );
% eta_ptwreab
f(81) = eta_ptwreab - ( eta_ptwreab_eq * mu_ptsodreab );
% mu_ptsodreab
musodreab_a = 0.12; musodreab_b = 10;
f(82) = mu_ptsodreab - ( musodreab_a * tanh(musodreab_b * (eta_ptsodreab/eta_ptsodreab_eq - 1)) + 1 );
% Phi_mdu
f(83) = Phi_mdu - ( Phi_gfilt - Phi_ptwreab );
% Phi_dtwreab
f(84) = Phi_dtwreab - ( Phi_mdu * eta_dtwreab );
% eta_dtwreab
f(85) = eta_dtwreab - ( eta_dtwreab_eq * mu_dtsodreab );
% mu_dtsodreab
f(86) = mu_dtsodreab - ( musodreab_a * tanh(musodreab_b * (eta_dtsodreab/eta_dtsodreab_eq - 1)) + 1 );
% Phi_dtu
f(87) = Phi_dtu - ( Phi_mdu - Phi_dtwreab );
% Phi_cdwreab
f(88) = Phi_cdwreab - ( Phi_dtu * eta_cdwreab );
% eta_cdwreab
f(89) = eta_cdwreab - ( eta_cdwreab_eq * mu_cdsodreab * mu_adh );
% mu_cdsodreab
f(90) = mu_cdsodreab - ( musodreab_a * tanh(musodreab_b * (eta_cdsodreab/eta_cdsodreab_eq - 1)) + 1 );
% mu_adh
muadh_a = 1.0328; muadh_b = 0.1938;
muadh_c = -1/4 * log((muadh_a - 1) / muadh_b);
f(91) = mu_adh - ( muadh_a - muadh_b * exp(-muadh_c * C_adh) );
% Phi_u
f(92) = Phi_u - ( Phi_dtu - Phi_cdwreab );
% Phi_win
if     fix_win
    f(93) = Phi_win - ( SSdata_input(93) );
else
    phiwin_a = 0.8; phiwin_c = 0.002313;
    phiwin_b = SSdata_input(47) + 1 / phiwin_a * log(phiwin_c*SF_U / 0.0150 - 1);
    f(93) = Phi_win - ( phiwin_c * SF_U / (1 + exp(-phiwin_a * (C_adh - phiwin_b))) );
end

end






























