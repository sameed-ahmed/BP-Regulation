% This sets an initial guess for all of the variables that is close to the
% steady state value and then saves it as a data file.
% 
% Values are from Karaaslan 2005, Leete 2018, and Munger 1988.

function solver_initial_guess_data

species = {'human', 'rat'   };
sex     = {'male' , 'female'};

for sex_ind = 1:2 % sex

% Species
spe_ind = 2;

%% Scaling factors

% Physiological variables which determine scaling factors.
if     strcmp(sex{sex_ind}, 'male')
    Phi_usod_new = 1.2212; % Munger 1988, Karaaslan 2005
    Phi_u_new    = 0.0150; % Munger 1988, Layton 2016
    R_r_new      = 15    ; % Munger 1988
    W_b          = 238   ; % Munger 1988
elseif strcmp(sex{sex_ind}, 'female')
    Phi_usod_new = 1.2212; % Munger 1988, Karaaslan 2005
    Phi_u_new    = 0.0150; % Munger 1988, Layton 2016
    R_r_new      = 24    ; % Munger 1988
    W_b          = 194   ; % Munger 1988
end
V_b_new = 0.06 * W_b + 0.77; % Lee 1985

% Rat value = Human value x SF
% Note: This includes conversion of units.
if     strcmp(species{spe_ind}, 'human')
    SF_S =                    1; % sodium flow
    SF_U =                    1; % urine flow
    SF_R =                    1; % resistance
    SF_V =                    1; % volume
elseif strcmp(species{spe_ind}, 'rat')
    SF_S = Phi_usod_new / 0.126; % sodium flow
    SF_U = Phi_u_new    / 0.001; % urine flow
    SF_R = R_r_new      / 83.3 ; % resistance
    SF_V = V_b_new      / 5    ; % volume
end

%% Manual entry of values

rsna          = 1; 
% alpha_map     = 100; % fixed var
alpha_map     = 1; 
alpha_rap     = 1; 
if     strcmp(sex{sex_ind}, 'male')
R_r           = 15; 
elseif strcmp(sex{sex_ind}, 'female')
R_r           = 24; 
end
beta_rsna     = 1; 
if     strcmp(sex{sex_ind}, 'male')
Phi_rb        = 6.5; 
elseif strcmp(sex{sex_ind}, 'female')
Phi_rb        = 4.2; 
end
if     strcmp(sex{sex_ind}, 'male')
Phi_gfilt     = 1.22; 
elseif strcmp(sex{sex_ind}, 'female')
Phi_gfilt     = 0.84; 
end
P_f           = 16; 
P_gh          = 62; 
% Sigma_tgf     = 3.859 * SF_S; % fixed var
Sigma_tgf     = 1; 
Phi_filsod    = Phi_gfilt * 143; 
if     strcmp(sex{sex_ind}, 'male')
% eta_ptsodreab = 0.93; % layton 2016
% eta_dtsodreab = 0.77; 
% eta_cdsodreab = 0.15; 
eta_ptsodreab = 0.8; % karaaslan
eta_dtsodreab = 0.5; 
eta_cdsodreab = 0.93;
elseif strcmp(sex{sex_ind}, 'female')
% eta_ptsodreab = 0.90; % layton 2016
% eta_dtsodreab = 0.77; 
% eta_cdsodreab = 0.15; 
% eta_ptsodreab = 0.71; % karaaslan
% eta_dtsodreab = 0.5; 
% eta_cdsodreab = 0.93;
% eta_ptsodreab = 0.5; % anita suggested
% eta_dtsodreab = 0.5; 
% eta_cdsodreab = 0.96;
eta_ptsodreab = 0.5; % calibrated
eta_dtsodreab = 0.5; 
eta_cdsodreab = 0.96;
end
Phi_ptsodreab = Phi_filsod * eta_ptsodreab; 
% gamma_filsod  = 14 * SF_S; % fixed var
gamma_filsod  = 1; 
gamma_at      = 1; 
gamma_rsna    = 1; 
Phi_mdsod     = Phi_filsod - Phi_ptsodreab; 
Phi_dtsodreab = Phi_mdsod * eta_dtsodreab; 
psi_al        = 1; 
Phi_dtsod     = Phi_mdsod - Phi_dtsodreab; 
Phi_cdsodreab = Phi_dtsod * eta_cdsodreab; 
% lambda_dt     = 2.03 * SF_S; % fixed var
lambda_dt     = 1; 
lambda_anp    = 1;
lambda_al     = 1;
Phi_usod      = Phi_dtsod - Phi_cdsodreab; 
Phi_sodin     = Phi_usod;
V_b           = 0.06 * W_b + 0.77; % Lee 1985
V_ecf         = (1/ (-0.4744) * log((2.4312*SF_V / (V_b - 4.5479*SF_V)) - 1) + 18.1128)*SF_V; 
P_mf          = 7;
Phi_vr        = 5 / SF_R; 
Phi_co        = 5 / SF_R; 
P_ra          = 0; 
vas           = 1; 
vas_f         = 0.01; 
vas_d         = 0.01; 
R_a           = 16.6 * SF_R;
R_ba          = 16.6 * SF_R; 
R_vr          = 1.4  * SF_R; 
R_tp          = 20   * SF_R; 
P_ma          = 100; 
epsilon_aum   = 1; 
% a_auto        = 0.011; % fixed var
a_auto        = 1; 
a_chemo       = 0.24; 
a_baro        = 1; 
C_adh         = 4.2; 
N_adh         = 1; 
% N_adhs        = 141; % fixed var
N_adhs        = 1; 
delta_ra      = 0; 
if     strcmp(sex{sex_ind}, 'male')
eta_ptwreab   = 0.86; 
eta_dtwreab   = 0.60; 
eta_cdwreab   = 0.78; 
elseif strcmp(sex{sex_ind}, 'female')
% eta_ptwreab   = 0.80; 
% eta_dtwreab   = 0.60; 
% eta_cdwreab   = 0.78; 
    eta_ptwreab = 0.5; % calibrated
    eta_dtwreab = 0.6; 
    eta_cdwreab = 0.91;
end
Phi_ptwreab   = Phi_gfilt * eta_ptwreab; 
mu_ptsodreab  = 1; 
Phi_mdu       = Phi_gfilt - Phi_ptwreab; 
Phi_dtwreab   = Phi_mdu * eta_dtwreab; 
mu_dtsodreab  = 1; 
Phi_dtu       = Phi_mdu - Phi_dtwreab; 
Phi_cdwreab   = Phi_dtu * eta_cdwreab; 
mu_cdsodreab  = 1; 
mu_adh        = 1; 
Phi_u         = Phi_dtu - Phi_cdwreab;
Phi_win       = Phi_u;
M_sod         = 2160 * SF_V; 
C_sod         = 143; 
% nu_mdsod      = 1.731 * SF_S; % fixed var
nu_mdsod      = 1; 
nu_rsna       = 1; 
C_al          = 387; % Pendergrass 2008
N_al          = 1; 
N_als         = 1; 
% xi_ksod       = 23.6; % fixed var
xi_ksod       = 1; 
xi_map        = 1; 
xi_at         = 1; 
hatC_anp      = 1; 
AGT           = 576005; 
nu_AT1        = 1; 
R_sec         = 1; 
PRC           = 17; 
PRA           = 135; 
AngI          = 90; 
AngII         = 6; 
AT1R          = 20; 
AT2R          = 7; 
Ang17         = 50; 
AngIV         = 1.3; 
if     strcmp(sex{sex_ind},  'male')
R_aa          = 6.0; 
R_ea          = 9.8;
elseif strcmp(sex{sex_ind},'female')
R_aa          = 9.4; 
R_ea          = 15.;
end
% Sigma_myo     = 62; % fixed var
Sigma_myo     = 1; 
Psi_AT1RAA    = 1; 
Psi_AT1REA    = 1; 
Psi_AT2RAA    = 1; 
Psi_AT2REA    = 1; 
%%

% Order
x  = [rsna; alpha_map; alpha_rap; R_r; beta_rsna; Phi_rb; Phi_gfilt; ...
      P_f; P_gh; Sigma_tgf; Phi_filsod; Phi_ptsodreab; eta_ptsodreab; ...
      gamma_filsod; gamma_at; gamma_rsna; Phi_mdsod; Phi_dtsodreab; ...
      eta_dtsodreab; psi_al; Phi_dtsod; Phi_cdsodreab; eta_cdsodreab; ...
      lambda_dt; lambda_anp; lambda_al; Phi_usod; Phi_sodin; V_ecf; ...
      V_b; P_mf; Phi_vr; Phi_co; P_ra; vas; vas_f; vas_d; R_a; R_ba; ...
      R_vr; R_tp; P_ma; epsilon_aum; a_auto; a_chemo; a_baro; C_adh; ...
      N_adh; N_adhs; delta_ra; ...
      M_sod; C_sod; nu_mdsod; nu_rsna; C_al; N_al; N_als; xi_ksod; ...
      xi_map; xi_at; hatC_anp; AGT; nu_AT1; R_sec; PRC; PRA; AngI; ...
      AngII; AT1R; AT2R; Ang17; AngIV; R_aa; R_ea; Sigma_myo; ...
      Psi_AT1RAA; Psi_AT1REA; Psi_AT2RAA; Psi_AT2REA; ...
      Phi_ptwreab; eta_ptwreab; mu_ptsodreab; ...
      Phi_mdu; Phi_dtwreab; eta_dtwreab; mu_dtsodreab; Phi_dtu; ...
      Phi_cdwreab; eta_cdwreab; mu_cdsodreab; mu_adh; Phi_u; Phi_win];

SSdataIG = x;

% Save values.
save_data_name = sprintf('%s_%s_ss_data_IG.mat', species{spe_ind},sex{sex_ind});
save_data_name = strcat('Data/', save_data_name);
save(save_data_name, 'SSdataIG')

end % sex

end

























