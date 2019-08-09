function pars = get_params(species,sex,AA)
if     strcmp(sex, 'male')
    gen = 1;
elseif strcmp(sex, 'female')
    gen = 0;
end
%disp([species{human+1},' ',gender{gg},' ',num2str(AA)])

% Scaling factor
% Rat flow = Human flow x SF
if strcmp(species, 'human') 
    SF = 1.0;
    human = 1;
elseif strcmp(species, 'rat') 
    SF = 4.5*10^(-3);
    human = 0;
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

if strcmp(species, 'human') 
    X_PRCPRA = 61/60.0; %fmol/min/pg
    if strcmp(sex, 'male')     
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

    end
    pt_sod_reab_EQ = 13.909;%15.285;%;14.4;%
    dt_sod_reab_EQ = 1.5859;%2.0714;%1.8;
    cd_sod_reab_EQ = 1.6909;%1.6934;%1.674;

 
elseif strcmp(species, 'rat') 
    % Male and female different parameters for RAS
    if     strcmp(sex, 'male')
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

    elseif strcmp(sex, 'female')
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
if     strcmp(sex, 'male')
    if     AA>1
        baseline = 0.0182;
    else
        baseline = 0.0193;
    end
else
    if     AA>1
        baseline = 0.0181;
    else
        baseline = 0.0199;
    end
end


pars = [N_rsna; R_aass; R_eass; P_B; P_go; C_gcf; eta_etapt; eta_epsdt; ...
        eta_etacd; K_vd; K_bar; R_bv; T_adh; Phi_sodin; C_K; T_al; ...
        N_rs; X_PRCPRA; h_renin; h_AGT; h_AngI; h_AngII; h_Ang17; ...
        h_AngIV; h_AT1R; h_AT2R; c_GPautoreg; P_ghnom; k_AGT; c_ACE; ...
        c_Chym; c_NEP; c_ACE2; c_IIIV; c_AT1R; c_AT2R; AT1R_eq; ...
        AT2R_eq; ALD_eq; gen; ...
        pt_sod_reab_EQ; dt_sod_reab_EQ; cd_sod_reab_EQ;amount;slope; baseline; human; SF];
end
