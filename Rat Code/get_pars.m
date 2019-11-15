% 

function pars = get_pars(species, gender, scenario)

%% Parameters

% Species boolean
if     strcmp(species, 'human')
    spe = 1;
elseif strcmp(species, 'rat')
    spe = 0;
end

% Gender boolean
% Sodium intake
% Renal vascular resistance
% Body weight (~blood volume)
% Original values are added as separate parameters because these may be
% modified by another script.
if     strcmp(gender, 'male')
    gen = 1;
    Phi_sodin_orig = 1.2212;
    Phi_u_orig = 0.0150;
    R_r_orig = 5.981 + 9.756;
    W_b = 238;
elseif strcmp(gender, 'female')
    gen = 0;
    Phi_sodin_orig = 1.2212;
    Phi_u_orig = 0.0150;
    R_r_orig = 9.361 + 15.27;
    W_b = 194;
end

N_rsna    = 1.00  ; % -
P_B       = 18    ; % mmHg
P_go      = 28    ; % mmHg
% K_vd      = 0.00001;
K_vd      = 0.01  ; % -
K_bar     = 16.6  ; % mmHg min / ml
R_bv      = 3.4   ; % mmHg min / ml
N_adhs_eq = 1     ; % -
T_adh     = 6     ; % min
Phi_sodin = 1.2212; % microEq / min
N_als_eq  = 1     ; % -
C_K       = 5     ; % microEq / ml 
T_al      = 60    ; % min
N_rs      = 1     ; % ng / ml / min

% RAS
h_renin       = 12  ; % min
h_AGT         = 600 ; % min
h_AngI        = 0.5 ; % min
h_AngII       = 0.66; % min
h_Ang17       = 30  ; % min
h_AngIV       = 0.5 ; % min
h_AT1R        = 12  ; % min
h_AT2R        = 12  ; % min
Psi_AT2RAA_eq = 1   ; % -
Psi_AT2REA_eq = 1   ; % -

% Species specific parameters
if     strcmp(species, 'human')
    
elseif strcmp(species, 'rat')
    if     strcmp(gender, 'male')
        R_aass = 5.981; % mmHg min / ml
        R_eass = 9.756; % mmHg min / ml
        C_gcf  = 0.068; % ml / min / mmHg
        
        % Transport parameters
        eta_ptsodreab_eq = 0.80; % karaaslan
        eta_dtsodreab_eq = 0.50; 
        eta_cdsodreab_eq = 0.93;
        eta_ptwreab_eq = 0.86; 
        eta_dtwreab_eq = 0.60; 
        eta_cdwreab_eq = 0.78;
        
        % RAS
        X_PRCPRA = 135.59/17.312   ; % 
        k_AGT    = 801.02          ; % 
        c_ACE    = 0.096833        ; % 
        c_Chym   = 0.010833        ; % 
        c_NEP    = 0.012667        ; % 
        c_ACE2   = 0.0026667       ; % 
        c_IIIV   = 0.29800         ; % 
        c_AT1R   = 0.19700         ; % 
        c_AT2R   = 0.065667        ; % 
        AT1R_eq  = 20.4807902818665; % 
        AT2R_eq  = 6.82696474842298; % 
    elseif strcmp(gender, 'female')
        R_aass = 9.361; % mmHg min / ml
        R_eass = 15.27; % mmHg min / ml
        C_gcf  = 0.047; % ml / min / mmHg
        
        % Transport parameters
        if   strcmp(scenario, 'm_Reab'         ) || ...
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
        
        if   strcmp(scenario, 'm_Reab'         ) || ...
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
        
        % RAS
        if   strcmp(scenario, 'm_RAS'         ) || ...
             strcmp(scenario, 'm_RAS_&_m_Reab')
        X_PRCPRA = 135.59/17.312   ; %  % male
        k_AGT    = 801.02          ; % 
        c_ACE    = 0.096833        ; % 
        c_Chym   = 0.010833        ; % 
        c_NEP    = 0.012667        ; % 
        c_ACE2   = 0.0026667       ; % 
        c_IIIV   = 0.29800         ; % 
        c_AT1R   = 0.19700         ; % 
        c_AT2R   = 0.065667        ; % 
        AT1R_eq  = 20.4807902818665; % 
        AT2R_eq  = 6.82696474842298; % 
        else
        X_PRCPRA = 114.22/17.312;
        k_AGT    = 779.63          ; % 
        c_ACE    = 0.11600         ; % 
        c_Chym   = 0.012833        ; % 
        c_NEP    = 0.0076667       ; % 
        c_ACE2   = 0.00043333      ; % 
        c_IIIV   = 0.29800         ; % 
        c_AT1R   = 0.19700         ; % 
        c_AT2R   = 0.065667        ; % 
        AT1R_eq  = 20.4538920068419; % 
        AT2R_eq  = 6.81799861123497; % 
        end
    end
end

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

% Retrieve and replace parameters in fixed variable equations.
% These are the shift parameters which ensure that effect variables are 1.
% Set name for data file to be loaded based upon gender.    
if   strcmp(scenario, 'm_Reab'         ) || ...
     strcmp(scenario, 'm_RSNA_&_m_Reab')
    load_data_name = sprintf('%s_fixed_var_pars_scenario_m_Reab.mat', gender);
else
    load_data_name = sprintf('%s_fixed_var_pars_scenario_Normal.mat', gender);
end
load(load_data_name, 'fixed_var_pars');

% Parameter input
pars = [N_rsna; R_aass; R_eass; P_B; P_go; C_gcf; eta_ptsodreab_eq; ...
        eta_dtsodreab_eq; eta_cdsodreab_eq; eta_ptwreab_eq; ...
        eta_dtwreab_eq; eta_cdwreab_eq; K_vd; K_bar; R_bv; ...
        N_adhs_eq; T_adh; Phi_sodin; N_als_eq; C_K; T_al; N_rs; ...
        X_PRCPRA; h_renin; h_AGT; h_AngI; h_AngII; h_Ang17; h_AngIV; ...
        h_AT1R; h_AT2R; k_AGT; c_ACE; c_Chym; c_NEP; c_ACE2; c_IIIV; ...
        c_AT1R; c_AT2R; AT1R_eq; AT2R_eq; Psi_AT2RAA_eq; Psi_AT2REA_eq; ...
        spe; gen; Phi_sodin_orig; Phi_u_orig; R_r_orig; W_b; ...
        fixed_var_pars];
    
end































