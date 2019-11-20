% 

function pars = get_pars(species, sex, scenario,hyp)

%% Parameters

% Species boolean
if     strcmp(species, 'human')
    spe = 1;
elseif strcmp(species, 'rat')
    spe = 0;
end
% Gender boolean
if     strcmp(sex, 'male')
    gen = 1;
elseif strcmp(sex, 'female')
    gen = 0;
end


%parameters in common between human and rat
N_rsna    = 1.00  ; % -
if hyp
    N_rsna = 2.5;
end
P_B       = 18    ; % mmHg
P_go      = 28    ; % mmHg
K_bar     = 16.6  ; % mmHg min / ml
R_bv      = 3.4   ; % mmHg min / ml
N_adhs_eq = 1     ; % -
T_adh     = 6     ; % min
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
    R_aass      = 31.67;   % mmHg min / l
    R_eass      = 51.66;   % mmHg min / l
    C_gcf       = 0.00781;
    K_vd      = 0.00001;
    Phi_sodin = 0.126;
    %Sodium transport
    eta_ptsodreab_eq = 0.80; % karaaslan
    eta_dtsodreab_eq = 0.50; 
    eta_cdsodreab_eq = 0.93;
    %RAAS
    X_PRCPRA = 61/60.0; %fmol/min/pg
    if  strcmp(sex, 'male')     
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
        if  hyp > 1 %hypertensive
            A_twreab = 0.0182;
        else %normotensive
            A_twreab = 0.0193;
        end
    elseif strcmp(sex, 'female')
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
        if  hyp > 1 %hypertensive
            A_twreab = 0.0181;
        else %normotensive
            A_twreab = 0.0199;
        end

    end
    pt_sod_reab_EQ = 13.909;
    dt_sod_reab_EQ = 1.5859;
    cd_sod_reab_EQ = 1.6909;

elseif strcmp(species, 'rat')
    if     strcmp(sex, 'male')
        Phi_sodin_orig = 1.2212;
        Phi_u_orig = 0.0150;
        R_r_orig = 5.981 + 9.756;
        W_b = 238;
        
        R_aass = 5.981; % mmHg min / ml
        R_eass = 9.756; % mmHg min / ml
        C_gcf  = 0.068; % ml / min / mmHg
        K_vd      = 0.01  ; % -
        
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
    elseif strcmp(sex, 'female')
        Phi_sodin_orig = 1.2212;
        Phi_u_orig = 0.0150;
        R_r_orig = 9.361 + 15.27;
        W_b = 194;
        
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
    % Add directory containing data.
    mypath = pwd;
    mypath = strcat(mypath, '/Data');
    addpath(genpath(mypath))

    % Retrieve and replace parameters in fixed variable equations.
    % These are the shift parameters which ensure that effect variables are 1.
    % Set name for data file to be loaded based upon sex.    
    if   strcmp(scenario, 'm_Reab'         ) || ...
         strcmp(scenario, 'm_RSNA_&_m_Reab')
        load_data_name = sprintf('%s_fixed_var_pars_scenario_m_Reab.mat', sex);
    else
        load_data_name = sprintf('%s_fixed_var_pars_scenario_Normal.mat', sex);
    end
    load(load_data_name, 'fixed_var_pars');
end

% Parameter input
pars = [spe; gen;...
        N_rsna; R_aass; R_eass; P_B; P_go; C_gcf; ...
        eta_ptsodreab_eq; eta_dtsodreab_eq; eta_cdsodreab_eq; ...     
        K_vd; K_bar; R_bv; ...
        N_adhs_eq; ... %rat
        T_adh; Phi_sodin; ...
        N_als_eq; ... %rat
        C_K; T_al; N_rs; ...
        X_PRCPRA; h_renin; h_AGT; h_AngI; h_AngII; h_Ang17; h_AngIV; ...
        h_AT1R; h_AT2R; k_AGT; c_ACE; c_Chym; c_NEP; c_ACE2; c_IIIV; ...
        c_AT1R; c_AT2R; AT1R_eq; AT2R_eq; ...
        Psi_AT2RAA_eq; Psi_AT2REA_eq; ...
         ];%...
if spe == 1 %human
    pars = [pars; ...
            ALD_eq;...
            pt_sod_reab_EQ; dt_sod_reab_EQ; cd_sod_reab_EQ; A_twreab]; %human water reabsorption
else %rat
    pars = [pars;...
            eta_ptwreab_eq; eta_dtwreab_eq; eta_cdwreab_eq; ...
            Phi_sodin_orig; Phi_u_orig; R_r_orig; W_b; ...
            fixed_var_pars];
end






























