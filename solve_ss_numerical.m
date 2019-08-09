% This script calculates the steady state solution to the system using
% fsolve. Some previous values are used as an initial guess. These are
% taken from Jessica, which are taken in part from some paper (Karaaslan
% 2005?).

function [exitflag,imag] = solve_ss_numerical(ACEI,diuretic,NSAID,IG)
human = 1;% 1 for human, 0 for rat
species = {'rat','human'};
gender     = {'male', 'female'};
AA = 2.5;

%SS_data_IG = zeros(82,2);
% X          = zeros(82,2);
% RESIDUAL   = zeros(82,2);
% EXITFLAG   = zeros(1 ,2);
% OUTPUT     = cell (1 ,2);

for gg = 1
%    for mm=0:0.25:1
%d_inhibit = mm
%% Parameters
pars = get_params(species{human+1},gender{gg},AA);

%% Drug Treatments
kappa_ACEI = ACEI;
kappa_d = 0;
kappa_d_tgf = 0;
kappa_d_renin = 0;


if ACEI ~= 0
    kappa_ACEI = 0.76;
end
if diuretic

    kappa_d = 0.15;
    kappa_d_tgf = 0.4;
    kappa_d_renin = 0.4;
end 


drugs = [kappa_ACEI,kappa_d,kappa_d_tgf,kappa_d_renin,NSAID];
%% Variables initial guess


SS_data_struct = load(IG,'SSdata');
SS_data_IG = SS_data_struct.SSdata;

  if length(SS_data_IG) ==83
       SS_data_IG(84) = 0.126;
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
x0 = SS_data_IG; x_p0 = zeros(84,1); t = 2000;

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
    
save_data_name = sprintf('%s_%s_ss_%s_%s_%s_rsna%s_fixedparams_highsodin_highphiwin_impairedmyo2.mat', species{human+1},gender{gg},num2str(ACEI),num2str(diuretic),num2str(NSAID),num2str(AA));
save_data_name = strcat('Data/', save_data_name);
save(save_data_name, 'SSdata', 'residual', 'exitflag', 'output')

% X(:,g) = x;
% RESIDUAL(:,g) = residual;
% EXITFLAG(g) = exitflag;
% OUTPUT{g} = output;

end
end
%end






























