% 

function pars = get_pars(gender, scenario)

%% Parameters

% Gender boolean
% Pars which determine scaling factors: 
%    urine sodium flow
%    renal vascular resistance
%    body weight (~blood volume)
if     strcmp(gender, 'male')
    gen = 1;
    Phi_usod_new = 1.2212;
    R_r_new = 28.6;
    W_b = 238;
elseif strcmp(gender, 'female')
    gen = 0;
    Phi_usod_new = 1.2212;
    R_r_new = 44.8;
    W_b = 194;
end
% Blood volume as a function of body weight from Lee 1985
V_b_new = 0.06 * W_b + 0.77;

% Scaling factors
% Rat value = Human value x SF
% Note: This includes conversion of units.
SF_S = Phi_usod_new / 0.126; % sodium flow
SF_R = R_r_new      / 83.3 ; % resistance
SF_V = V_b_new      / 5    ; % volume

N_rsna     = 1.00;
% R_aass    = 31.67 / SF;   % mmHg min / ml
% R_eass    = 51.66 / SF;   % mmHg min / ml
if     strcmp(gender, 'male')
    R_aass = 10.87;    % mmHg min / ml
    R_eass = 17.74;    % mmHg min / ml
elseif strcmp(gender, 'female')
    R_aass = 17.02;    % mmHg min / ml
    R_eass = 27.76;    % mmHg min / ml
end
P_B        = 18;           % mmHg
P_go       = 28;           % mmHg
% C_gcf     = 0.00781 * SF;
if     strcmp(gender, 'male')
    C_gcf  = 0.068;    % ml / min / mmHg
elseif strcmp(gender, 'female')
    C_gcf  = 0.047;    % ml / min / mmHg
end

% Male and female different parameters for fractional reabsorption
if     strcmp(gender, 'male')
%     eta_ptsodreab_eq = 0.93;  % layton 2016
%     eta_dtsodreab_eq = 0.77; 
%     eta_cdsodreab_eq = 0.15;
    eta_ptsodreab_eq = 0.80; % karaaslan
    eta_dtsodreab_eq = 0.50; 
    eta_cdsodreab_eq = 0.93;
elseif strcmp(gender, 'female')
    if     strcmp(scenario, 'm_Reab'         ) || ...
           strcmp(scenario, 'm_RAS_&_m_Reab' ) || ...
           strcmp(scenario, 'm_RSNA_&_m_Reab')
    eta_ptsodreab_eq = 0.71; % male
    eta_dtsodreab_eq = 0.50; 
    eta_cdsodreab_eq = 0.93;
    else
    eta_ptsodreab_eq = 0.50; % calibrated
    eta_dtsodreab_eq = 0.50; 
    eta_cdsodreab_eq = 0.96;
    end
%     eta_ptsodreab_eq = 0.90; % layton 2016
%     eta_dtsodreab_eq = 0.77; 
%     eta_cdsodreab_eq = 0.15;
%     eta_ptsodreab_eq = 0.50; % anita suggested
%     eta_dtsodreab_eq = 0.50; 
%     eta_cdsodreab_eq = 0.96;
end
if     strcmp(gender, 'male')
    eta_ptwreab_eq = 0.86; 
    eta_dtwreab_eq = 0.60; 
    eta_cdwreab_eq = 0.78;
elseif strcmp(gender, 'female')
    if     strcmp(scenario, 'm_Reab'         ) || ...
           strcmp(scenario, 'm_RAS_&_m_Reab' ) || ...
           strcmp(scenario, 'm_RSNA_&_m_Reab')
    eta_ptwreab_eq = 0.80; % male 
    eta_dtwreab_eq = 0.60; 
    eta_cdwreab_eq = 0.78;
    else
    eta_ptwreab_eq = 0.50; % calibrated
    eta_dtwreab_eq = 0.60; 
    eta_cdwreab_eq = 0.91;
    end
end

% K_vd      = 0.00001;
K_vd      = 0.01;
% K_bar     = 16.6 / SF;  % mmHg min / ml
K_bar     = 16.6 * SF_R;  % mmHg min / ml
% R_bv      = 3.4 / SF;   % mmHg min / ml
R_bv      = 3.4 * SF_R;   % mmHg min / ml
N_adhs_eq = 1;
T_adh     = 6;            % min
% Phi_sodin = 1.2278;       % microEq / min % old
% Phi_sodin = 2.3875;       % microEq / min % layton 2016
Phi_sodin = 1.2212;       % microEq / min % karaaslan
N_als_eq  = 1;
C_K       = 5;            % microEq / ml 
T_al      = 30;           % min LISTED AS 30 IN TABLE %listed as 60 in text will only change dN_al
N_rs      = 1;            % ng / ml / min

% RAS
h_renin   = 12;      % min
h_AGT     = 10*60;   % min
h_AngI    = 0.5;     % min
h_AngII   = 0.66;    % min
h_Ang17   = 30;      % min
h_AngIV   = 0.5;     % min
h_AT1R    = 12;      % min
h_AT2R    = 12;      % min

% Male and female different parameters for RAS
if     strcmp(gender, 'male')
    X_PRCPRA = 135.59/17.312;
    k_AGT    = 801.02;
    c_ACE    = 0.096833;
    c_Chym   = 0.010833;
    c_NEP    = 0.012667;
    c_ACE2   = 0.0026667;
    c_IIIV   = 0.29800;
    c_AT1R   = 0.19700;
    c_AT2R   = 0.065667;
    AT1R_eq  = 20.4807902818665;
    AT2R_eq  = 6.82696474842298;
elseif strcmp(gender, 'female')
    if     strcmp(scenario, 'm_RAS'         ) || ...
           strcmp(scenario, 'm_RAS_&_m_Reab')
    X_PRCPRA = 135.59/17.312; % male
    k_AGT    = 801.02;
    c_ACE    = 0.096833;
    c_Chym   = 0.010833;
    c_NEP    = 0.012667;
    c_ACE2   = 0.0026667;
    c_IIIV   = 0.29800;
    c_AT1R   = 0.19700;
    c_AT2R   = 0.065667;
    AT1R_eq  = 20.4807902818665;
    AT2R_eq  = 6.82696474842298;
    else
    X_PRCPRA = 114.22/17.312;
    k_AGT    = 779.63;
    c_ACE    = 0.11600;
    c_Chym   = 0.012833;
    c_NEP    = 0.0076667;
    c_ACE2   = 0.00043333;
    c_IIIV   = 0.29800;
    c_AT1R   = 0.19700;
    c_AT2R   = 0.065667;
    AT1R_eq  = 20.4538920068419;
    AT2R_eq  = 6.81799861123497;
    end
end

Psi_AT2RAA_eq = 1;
Psi_AT2REA_eq = 1;

% Parameter input.
pars = [N_rsna; R_aass; R_eass; P_B; P_go; C_gcf; eta_ptsodreab_eq; ...
        eta_dtsodreab_eq; eta_cdsodreab_eq; eta_ptwreab_eq; ...
        eta_dtwreab_eq; eta_cdwreab_eq; K_vd; K_bar; R_bv; N_adhs_eq; ...
        T_adh; Phi_sodin; N_als_eq; C_K; T_al; N_rs; X_PRCPRA; h_renin; ...
        h_AGT; h_AngI; h_AngII; h_Ang17; h_AngIV; h_AT1R; h_AT2R; ...
        k_AGT; c_ACE; c_Chym; c_NEP; c_ACE2; c_IIIV; c_AT1R; c_AT2R; ...
        AT1R_eq; AT2R_eq; Psi_AT2RAA_eq; Psi_AT2REA_eq; ...
        gen; SF_S; SF_R; SF_V];
    
end































