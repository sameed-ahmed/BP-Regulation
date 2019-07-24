% This script calculates the steady state solution to the system using
% fsolve. Some previous values are used as an initial guess. These are
% taken from Jessica, which are taken in part from some paper (Karaaslan
% 2005?).

function [exitflag,imag] = solve_ss_numerical(ACEI,diuretic,NSAID,IG)
human = 1;% 1 for human, 0 for rat
species = {'rat','human'};
gender     = {'male', 'female'};
AA = 1;

%SS_data_IG = zeros(82,2);
% X          = zeros(82,2);
% RESIDUAL   = zeros(82,2);
% EXITFLAG   = zeros(1 ,2);
% OUTPUT     = cell (1 ,2);

for gg = 1
%    for mm=0:0.25:1
%d_inhibit = mm
%% Parameters

if     strcmp(gender{gg}, 'male')
    gen = 1;
elseif strcmp(gender{gg}, 'female')
    gen = 0;
end
%disp([species{human+1},' ',gender{gg},' ',num2str(AA)])

% Scaling factor
% Rat flow = Human flow x SF
if human
    SF = 1.0;
else
    SF = 4.5*10^(-3);
end

N_rsna      = 1*AA;
R_aass      = 31.67 / SF;   % mmHg min / l
R_eass      = 51.66 / SF;   % mmHg min / l
P_B         = 18;           % mmHg
P_go        = 28;           % mmHg
C_gcf       = 0.00781 * SF;
eta_etapt   = 0.8; 
eta_epsdt   = 0.5; 
eta_etacd   = 0.93; 
K_vd        = 0.00001;
K_bar       = 16.6 / SF;    % mmHg min / l
R_bv        = 3.4 / SF;     % mmHg min / l
T_adh       = 6;            % min
Phi_sodin   = 0.126 * SF;   % mEq / min
%if     strcmp(gender{gg}, 'male')
%   C_K      = 6;            % mEq / l 
%elseif strcmp(gender{gg}, 'female')
   C_K      = 5;            % mEq / l 
%end
T_al        = 30;           % min LISTED AS 30 IN TABLE %listed as 60 in text will only change dN_al
N_rs        = 1;            % ng / ml / min

% RAS
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

if human
    X_PRCPRA = 61/60.0; %fmol/min/pg
    if strcmp(gender{gg}, 'male')     
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
        ALD_eq = 85;
        amount = 0.0831;
        slope = 8.5210;
        baseline = 0.018;
    elseif strcmp(gender{gg}, 'female')
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
        ALD_eq = 69.1775;
        amount = 0.0836;
        slope = 8.4953;
        baseline = 0.02;
    end
    pt_sod_reab_EQ = 13.909;%15.285;%;14.4;%
    dt_sod_reab_EQ = 1.5859;%2.0714;%1.8;
    cd_sod_reab_EQ = 1.6909;%1.6934;%1.674;

 
else %rat
    % Male and female different parameters for RAS
    if     strcmp(gender{gg}, 'male')
        C_K      = 6;
        X_PRCPRA = 135.59/17.312;
        k_AGT   = 801.02;
        c_ACE   = 0.096833;
        c_Chym  = 0.010833;
        c_NEP   = 0.012667;
        c_ACE2  = 0.0026667;
        c_IIIV  = 0.29800;
        c_AT1R  = 0.19700;
        c_AT2R  = 0.065667;
        AT1R_eq = 20.46;
        AT2R_eq = 6.82;
        ALD_eq = 395.3;

    elseif strcmp(gender{gg}, 'female')
        X_PRCPRA = 114.22/17.312;
        k_AGT   = 779.63;
        c_ACE   = 0.11600;
        c_Chym  = 0.012833;
        c_NEP   = 0.0076667;
        c_ACE2  = 0.00043333;
        c_IIIV  = 0.29800;
        c_AT1R  = 0.19700;
        c_AT2R  = 0.065667;
        AT1R_eq = 20.46;
        AT2R_eq = 6.82;
        ALD_eq = 379.4;

    end
    pt_sod_reab_EQ = 0.068294039337572;
    dt_sod_reab_EQ = 0.008094477703862;
    cd_sod_reab_EQ = 0.007616239742696;
end
amount = 0.08;
slope = 8.5;
%baseline = 0.02;
if     strcmp(gender{gg}, 'male')
    if     AA>1
        baseline = 0.0183;
    else
        baseline = 0.0199;
    end
else
    if     AA>1
        baseline = 0.01843;
    else
        baseline = 0.0222;
    end
end


pars = [N_rsna; R_aass; R_eass; P_B; P_go; C_gcf; eta_etapt; eta_epsdt; ...
        eta_etacd; K_vd; K_bar; R_bv; T_adh; Phi_sodin; C_K; T_al; ...
        N_rs; X_PRCPRA; h_renin; h_AGT; h_AngI; h_AngII; h_Ang17; ...
        h_AngIV; h_AT1R; h_AT2R; c_GPautoreg; P_ghnom; k_AGT; c_ACE; ...
        c_Chym; c_NEP; c_ACE2; c_IIIV; c_AT1R; c_AT2R; AT1R_eq; ...
        AT2R_eq; ALD_eq; gen; ...
        pt_sod_reab_EQ; dt_sod_reab_EQ; cd_sod_reab_EQ;amount;slope; baseline; human; SF];
%% Drug treatments
kappa_ACEI = ACEI;
kappa_d = 0;
kappa_d_tgf = 0;
kappa_d_renin = 0;


if ACEI ~= 0
    kappa_ACEI = 0.77;
end
if diuretic

           kappa_d = 0.15;
    kappa_d_tgf = 0.4;
    kappa_d_renin = 0.4;
end 


drugs = [kappa_ACEI,kappa_d,kappa_d_tgf,kappa_d_renin,NSAID];
%% Variables initial guess

% Add directory containing data.
% mypath = pwd;
% mypath = strcat(mypath, '/Data');
% addpath(genpath(mypath))

% Load data for steady state initial value. 
% Need to first run transform_data.m on Jessica's data files.
%if     strcmp(gender{gg}, 'male')
%    SS_data_IG(:,gg) = csvread(  '../00_Rat Stuff for Jessica/Rat Code for Jessica/Data/male_ss_data_IG.txt');
%elseif strcmp(gender{gg}, 'female')
%    SS_data_IG(:,gg) = csvread('../00_Rat Stuff for Jessica/Rat Code for Jessica/Data/female_ss_data_IG.txt');
%end


%SS_data_struct = load(sprintf('%s_ss_data_%s_%s_%s.mat', gender{gg},num2str(starts(1)*ACEI),num2str(starts(2)*diuretic),num2str(starts(3)*NSAID)),'SSdata');


SS_data_struct = load(IG,'SSdata');
SS_data_IG = SS_data_struct.SSdata;
%disp([IG,' ',num2str(length(SS_data_IG))])
 if length(SS_data_IG) == 82
     disp('took correct path')
     SS_data_IG_new(1:52) = SS_data_IG(1:52);
     SS_data_IG_new(53) = 1;
     SS_data_IG_new(54:83) = SS_data_IG(53:82);
 elseif length(SS_data_IG) ==83
     SS_data_IG_new = SS_data_IG;
 else
     disp('wrong IG size')
     return
 end
 SS_data_IG(84) = 0.126;

% 
%  time_data_struct = load(IG,'x');
%  x_data = time_data_struct.x;
%  %size(x_data)
%  SS_data_IG = transpose(x_data(117,:));%
%  size(SS_data_IG)
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
x0 = SS_data_IG; x_p0 = zeros(84,1); t = 0;

%% Find steady state solution
tchange=0;

%options = optimset(); options = optimset('Display','off');
options = optimset('Display','off','MaxFunEvals',8200+10000);
[SSdata, residual, ...
 exitflag, output] = fsolve(@(x) blood_press_reg_sim(t,x,x_p0,pars,tchange,drugs), ...
                            x0, options);

% Check for solver convergence.
if exitflag == 0
    %disp('Solver did not converge.')
    %disp(output)
end

% Check for imaginary solution.
if not (isreal(SSdata))
    %disp('Imaginary number returned.')
    imag = 1;
else
    imag = 0;
end

% Set any values that are within machine precision of 0 equal to 0.
for i = 1:length(SSdata)
    if abs(SSdata(i)) < eps*100
        SSdata(i) = 0;
    end
end
    
save_data_name = sprintf('%s_%s_ss_%s_%s_%s_rsna%s_fixedparams.mat', species{human+1},gender{gg},num2str(ACEI),num2str(diuretic),num2str(NSAID),num2str(AA));
save_data_name = strcat('Data/', save_data_name);
save(save_data_name, 'SSdata', 'residual', 'exitflag', 'output')

% X(:,g) = x;
% RESIDUAL(:,g) = residual;
% EXITFLAG(g) = exitflag;
% OUTPUT{g} = output;

end
end
%end






























