% Baseline steady state data for the variables is given by Jessica's code
% in folder 01_ang_inf. Copy those data files to the folder Data in this 
% directory. Then run this script.
% 
% This script reodrders the variables in the data to fit the order of the
% variables in the code in this directory. It then saves it with a suitable
% name.

gender = {'male', 'female'};

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

for gg = 1:2 % gender

% Scaling factor
% Rat flow = Human flow x SF
if     strcmp(gender{gg}, 'male')
    SF = 4.5*10^(-3)*10^(3);
elseif strcmp(gender{gg}, 'female')
    SF = 2/3 * 4.5*10^(-3)*10^(3);
end

% Load data.
% if     strcmp(gender{i}, 'male')
%     SSdata = csvread(  'male_RAS1_ALD1_AT2R0_normo_0_0_0_0_0.txt');
% elseif strcmp(gender{i}, 'female')
%     SSdata = csvread('female_RAS0_ALD0_AT2R1_normo_0_0_0_0_0.txt');
% end
if     strcmp(gender{gg}, 'male')
    SSdata = csvread(  'male_RAS1_ALD1_normo_0_0_0_0_0_IG.txt');
elseif strcmp(gender{gg}, 'female')
    SSdata = csvread('female_RAS0_ALD0_normo_0_0_0_0_0_IG.txt');
end

% Rows 1-7 are some other quantities.
% SSdata(1:7) = '';
SSdata(1:6) = '';

% Jessica's order
alpha_map     = SSdata(1 );
alpha_rap     = SSdata(2 );
rsna          = SSdata(3 );
beta_rsna     = SSdata(4 );
R_aa          = SSdata(5 ) / SF;
R_ea          = SSdata(6 ) / SF;
R_r           = SSdata(7 ) / SF;
Phi_rb        = SSdata(8 ) * SF;
P_gh          = SSdata(9 );
P_f           = SSdata(10);
Phi_gfilt     = SSdata(11) * SF;
Sigma_tgf     = SSdata(12);
Phi_filsod    = SSdata(13) * SF;
gamma_filsod  = SSdata(14);
gamma_at      = SSdata(15);
gamma_rsna    = SSdata(16);
eta_ptsodreab = SSdata(17);
Phi_ptsodreab = SSdata(18) * SF;
Phi_mdsod     = SSdata(19) * SF;
psi_al        = SSdata(20);
eta_dtsodreab = SSdata(21);
Phi_dtsodreab = SSdata(22) * SF;
Phi_dtsod     = SSdata(23) * SF;
lambda_dt     = SSdata(24);
lambda_anp    = SSdata(25);
eta_cdsodreab = SSdata(26);
Phi_cdsodreab = SSdata(27) * SF;
Phi_usod      = SSdata(28) * SF;
Phi_win       = SSdata(29) * SF;
V_b           = SSdata(30) * SF;
P_mf          = SSdata(31);
Phi_vr        = SSdata(32) * SF;
Phi_co        = SSdata(33) * SF;
P_ra          = SSdata(34);
vas_f         = SSdata(35);
vas_d         = SSdata(36);
R_ba          = SSdata(37) / SF;
R_a           = SSdata(38) / SF;
R_vr          = SSdata(39) / SF;
R_tp          = SSdata(40) / SF;
P_ma          = SSdata(41);
a_auto        = SSdata(42);
a_chemo       = SSdata(43);
epsilon_aum   = SSdata(44);
N_adhs        = SSdata(45);
C_adh         = SSdata(46);
mu_al         = SSdata(47);
mu_adh        = SSdata(48);
Phi_twreab    = SSdata(49) * SF;
Phi_u         = SSdata(50) * SF;
C_sod         = SSdata(51);
nu_mdsod      = SSdata(52);
nu_rsna       = SSdata(53);
xi_ksod       = SSdata(54);
xi_map        = SSdata(55);
xi_at         = SSdata(56);
N_als         = SSdata(57);
C_al          = SSdata(58);
hatC_anp      = SSdata(59);
V_ecf         = SSdata(60) * SF;
vas           = SSdata(61);
a_baro        = SSdata(62);
N_adh         = SSdata(63);
delta_ra      = SSdata(64);
M_sod         = SSdata(65) * SF;
N_al          = SSdata(66);
R_sec         = SSdata(67);
nu_AT1        = SSdata(68);
PRC           = SSdata(69);
PRA           = SSdata(70);
AGT           = SSdata(71);
AngI          = SSdata(72);
AngII         = SSdata(73);
Ang17         = SSdata(74);
AngIV         = SSdata(75);
AT1R          = SSdata(76);
AT2R          = SSdata(77);
Psi_AT1RAA    = SSdata(78);
Psi_AT1REA    = SSdata(79);
Sigma_myo     = SSdata(80);
Psi_AT2RAA    = SSdata(81);
Psi_AT2REA    = SSdata(82);

% My order
x  = [rsna; alpha_map; alpha_rap; R_r; beta_rsna; Phi_rb; Phi_gfilt; ...
      P_f; P_gh; Sigma_tgf; Phi_filsod; Phi_ptsodreab; eta_ptsodreab; ...
      gamma_filsod; gamma_at; gamma_rsna; Phi_mdsod; Phi_dtsodreab; ...
      eta_dtsodreab; psi_al; Phi_dtsod; Phi_cdsodreab; eta_cdsodreab; ...
      lambda_dt; lambda_anp; Phi_usod; Phi_win; V_ecf; V_b; P_mf; ...
      Phi_vr; Phi_co; P_ra; vas; vas_f; vas_d; R_a; R_ba; R_vr; R_tp; ...
      P_ma; epsilon_aum; a_auto; a_chemo; a_baro; C_adh; N_adh; ...
      N_adhs; delta_ra; Phi_twreab; mu_al; mu_adh; Phi_u; M_sod; ...
      C_sod; nu_mdsod; nu_rsna; C_al; N_al; N_als; xi_ksod; xi_map; ...
      xi_at; hatC_anp; AGT; nu_AT1; R_sec; PRC; PRA; AngI; AngII; ...
      AT1R; AT2R; Ang17; AngIV; R_aa; R_ea; Sigma_myo; Psi_AT1RAA; ...
      Psi_AT1REA; Psi_AT2RAA; Psi_AT2REA];

SSdataIG = x;

save_data_name = sprintf('%s_ss_data_IG.mat', gender{gg});
save_data_name = strcat('Data/', save_data_name);
save(save_data_name, 'SSdataIG')

end % gender

clear
































