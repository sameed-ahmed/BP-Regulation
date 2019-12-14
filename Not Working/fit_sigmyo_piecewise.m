% This runs run_sim_P_ma.m in such a way that the parameters in Sigma_myo
% are fit using fmincon to the GFR vs MAP graph.

function [pars_min, ess_min, exitflag, output] = fit_sig_myo_piecewise

%% User input 

% These are roughly taken from some source. Adjust and cite later.

% Autoregulatory range
ARR = [90,160];
% Slope of GFR vs MAP relative to baseline in autoregulatory range, i.e.,
% (deltaGFR / deltaMAP) / GFRbaseline.
% Units are ([l/min] / [mm Hg]) / [l/min].
slope = 0.0069/1000;
% Factor to increase/decrease GFR by when 30 mm Hg outside of 
% autoregulatory range.
fact = 1/2;

%% Parameters initial guess 

% if     P_gh <= x1
%     Sigma_myo = y1;
% elseif x1 < P_gh < x2
%     Sigma_myo = (y2 - y1)/(x2 - x1)*P_gh + (x2*y1 - x1*y2)/(x2 - x1);
% elseif x2 <= P_gh
%     Sigma_myo = y2;
% end
x1 = 61.93  ; % left position
x2 = 62.64  ; % right position
y1 =  0.2818; % bottom position
y2 =  7.615 ; % top position

pars_est = [x1; x2; y1; y2];

gender = {'male', 'female'};

%%

for gg = 1:1 % gender

%% Parameters 

if     strcmp(gender{gg}, 'male')
    gen = 1;
elseif strcmp(gender{gg}, 'female')
    gen = 0;
end

% Scaling factor
% Rat flow = Human flow x SF
if     strcmp(gender{gg}, 'male')
    SF = 4.5*10^(-3);
elseif strcmp(gender{gg}, 'female')
    SF = 2/3 * 4.5*10^(-3);
end

N_rsna       = 1;
R_aass       = 31.67 / SF;   % mmHg min / l
R_eass       = 51.66 / SF;   % mmHg min / l
P_B          = 18;           % mmHg
P_go         = 28;           % mmHg
C_gcf        = 0.00781 * SF;
if     strcmp(gender{gg}, 'male')
   eta_etapt = 0.8; 
elseif strcmp(gender{gg}, 'female')
   eta_etapt = 0.5; 
end
eta_epsdt    = 0.5; 
if     strcmp(gender{gg}, 'male')
   eta_etacd = 0.93; 
elseif strcmp(gender{gg}, 'female')
   eta_etacd = 0.972; 
end
K_vd         = 0.00001;
K_bar        = 16.6 / SF;    % mmHg min / l
R_bv         = 3.4 / SF;     % mmHg min / l
T_adh        = 6;            % min
Phi_sodin    = 0.126 * SF;   % mEq / min
C_K          = 5;            % mEq / l 
T_al         = 30;           % min LISTED AS 30 IN TABLE %listed as 60 in text will only change dN_al
N_rs         = 1;            % ng / ml / min

% RAS
h_renin      = 12;      % min
h_AGT        = 10*60;   % min
h_AngI       = 0.5;     % min
h_AngII      = 0.66;    % min
h_Ang17      = 30;      % min
h_AngIV      = 0.5;     % min
h_AT1R       = 12;      % min
h_AT2R       = 12;      % min

% Male and female different parameters for RAS
if     strcmp(gender{gg}, 'male')
    X_PRCPRA = 135.59/17.312;
    k_AGT    = 801.02;
    c_ACE    = 0.096833;
    c_Chym   = 0.010833;
    c_NEP    = 0.012667;
    c_ACE2   = 0.0026667;
    c_IIIV   = 0.29800;
    c_AT1R   = 0.19700;
    c_AT2R   = 0.065667;
    AT1R_eq  = 20.46;
    AT2R_eq  = 6.82;
elseif strcmp(gender{gg}, 'female')
    X_PRCPRA = 114.22/17.312;
    k_AGT    = 779.63;
    c_ACE    = 0.11600;
    c_Chym   = 0.012833;
    c_NEP    = 0.0076667;
    c_ACE2   = 0.00043333;
    c_IIIV   = 0.29800;
    c_AT1R   = 0.19700;
    c_AT2R   = 0.065667;
    AT1R_eq  = 20.46;
    AT2R_eq  = 6.82;
end

pars = [N_rsna; R_aass; R_eass; P_B; P_go; C_gcf; eta_etapt; eta_epsdt; ...
        eta_etacd; K_vd; K_bar; R_bv; T_adh; Phi_sodin; C_K; T_al; ...
        N_rs; X_PRCPRA; h_renin; h_AGT; h_AngI; h_AngII; h_Ang17; ...
        h_AngIV; h_AT1R; h_AT2R; k_AGT; c_ACE; c_Chym; c_NEP; c_ACE2; ...
        c_IIIV; c_AT1R; c_AT2R; AT1R_eq; AT2R_eq; SF; gen];

%% Load steady state data initial guess 

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

% Load data for steady state initial value. 
if     strcmp(gender{gg},   'male')
    SSdata_IG = csvread(  'male_ss_data_IG.txt');
elseif strcmp(gender{gg}, 'female')
    SSdata_IG = csvread('female_ss_data_IG.txt');
end

%% Optimization 

% Place holders for fmincon.
Amat = []; bvec = []; Aeq = []; beq = []; nonlcon = @base_value; 

% Lower and upper bounds for parameters in fmincon.
% pars_est = [x1; x2; y1; y2];
lb         = [60; 62;  0;  0];
ub         = [62; 66;  1; 20];

% Edit options for optimizer.
options = optimset('Display','iter');
% options = optimset('Display','iter', 'MaxFunEvals',100, ...
%                    'MaxIter',100, 'FinDiffRelStep',sqrt(eps));

% Find minimizing parameters.
tic
[pars_min, ess_min, exitflag, output] = ...
    fmincon(@(pars_est) cost_fun(pars,pars_est,SSdata_IG,ARR,slope,fact),pars_est, ...
            Amat,bvec,Aeq,beq,lb,ub,nonlcon, options);
CPUtime = toc

%% Save result

save_data_name = sprintf('%s_pars_est_sigmamyo_piecewise.mat', gender{gg});
save(save_data_name, 'pars_min', 'ess_min', 'exitflag', 'output', 'CPUtime')

end % gender

end

% =========================================================================
% Cost function to calculate error between actual and calculated GFR.
% =========================================================================

function err = cost_fun(pars,pars_est,SSdata_IG,ARR,slope,fact)

%% GFR value and GFR/MAP slope in autoregulatory range.

gen = pars(end);
if     gen == 1
    gender = 'male';
elseif gen == 0
    gender = 'female';
end
if     strcmp(gender,'male')
    GFR_baseline = 5.6576e-04;
elseif strcmp(gender,'female')
    GFR_baseline = 3.7751e-04;
end
slope = slope * GFR_baseline;
deltaMAP = 5;
deltaGFR = slope * deltaMAP;

%% Form MAP

% MAP autoregulatory range
ARR = (ARR(1):deltaMAP:ARR(2))'; 
middle_num = length(ARR); middle_ind = ARR - (ARR(1)-30-1);
% MAP whole range
MAP_lower1 = ARR(1  )-30; MAP_lower2 = ARR(1  )-15; 
MAP_upper1 = ARR(end)+15; MAP_upper2 = ARR(end)+30;
MAP_range = [MAP_lower1; MAP_lower2; ARR; MAP_upper1; MAP_upper2];

%% Form GFR actual

% Increase/decrease linearly by factor when outside of autoregulatory range. 
GFR_act_lower1 = (1-fact  ) * GFR_baseline; GFR_act_lower2 = (1-fact/2) * GFR_baseline;
GFR_act_upper1 = (1+fact/2) * GFR_baseline; GFR_act_upper2 = (1+fact  ) * GFR_baseline;

GFR_act = zeros(middle_num,1) + GFR_baseline;
for i = 2:middle_num
    GFR_act(i) = GFR_act(i) + (i-1)*deltaGFR;
end
GFR_act = [GFR_act_lower1; GFR_act_lower2; GFR_act; GFR_act_upper1; GFR_act_upper2];

%% Modify MAP and GFR calculated

% Retrieve GFR calculated and steady state.
[MAP, GFR] = run_sim_P_ma(pars,pars_est,SSdata_IG);
% [MAP, GFR, SSdata] = run_sim_P_ma; % quickly test this function only

% Round values. Delete overlapping values.
MAP = round(MAP);
for i = length(MAP):-1:2
    if MAP(i) == MAP(i-1)
        MAP(i) = '';
        GFR(i) = '';
    end
end

% Trim excess.
for i = length(MAP):-1:1
    if MAP(i) > MAP_upper2 || MAP(i) < MAP_lower1
        MAP(i) = '';
        GFR(i) = '';
    end
end
if min(MAP) > MAP_lower1 || max(MAP) < MAP_upper2
    err = 1;
    return;
end

%% Form GFR calculated

GFR_cal_lower1 = GFR(1     ); GFR_cal_lower2 = GFR(1+15); 
GFR_cal_upper1 = GFR(end-15); GFR_cal_upper2 = GFR(end );
GFR_cal = GFR(middle_ind);
GFR_cal = [GFR_cal_lower1; GFR_cal_lower2; GFR_cal; GFR_cal_upper1; GFR_cal_upper2];

%% Compute error

err = norm((GFR_cal - GFR_act)./GFR_act);

end

% =========================================================================
% Run simulation of all pars vs MAP.
% =========================================================================

function [MAP, GFR] = run_sim_P_ma(pars,pars_est,SSdata_IG)

lower_time = 130;
upper_time = 160;

% X = (iteration, variable)
X = zeros(lower_time+1+upper_time,82);
% Initial condition for the variables and their derivatives. 
% System is initially at steady state, so the derivative is 0.
SSdata = solve_ss_numerical(pars,pars_est,SSdata_IG);
x0 = SSdata; x_p0 = zeros(82,1);
X(lower_time+1,:) = x0;

change = {'decrease', 'increase'};
for cc = 1:2 % change

if     strcmp(change{cc}, 'decrease')
    ch = 1;
elseif strcmp(change{cc}, 'increase')
    ch = 0;
end

pars = [pars; ch];

%% Solve DAE

% Time at which to keep steady state.
if     strcmp(change{cc}, 'decrease')
    tss = linspace(1,lower_time,lower_time);
elseif strcmp(change{cc}, 'increase')
    tss = linspace(1,upper_time,upper_time);
end
iter_range = length(tss);

for iter = 1:iter_range

% Run simulation to get steady state pressure value. ----------------------

% Initial time (min); Final time (min); Time vector;
t0 = 0; tf = tss(iter); tspan = [t0, tf];
% Options
options = odeset('MaxStep',50); % default is 0.1*abs(t0-tf)

% Solve dae
[~,x] = ode15i(@(t,x,x_p) ...
               bp_reg_sim_P_ma_1st(t,x,x_p, pars,pars_est), ...
               tspan, x0, x_p0, options);
P_ma_ss = x(end,41);

% Run simulation to get steady state values of everything else. -----------

% Initial time (min); Final time (min); Time vector;
t0 = 0; tf = tss(iter) + 50; tspan = [t0, tf];

% Solve dae
[~,x] = ode15i(@(t,x,x_p) ...
               bp_reg_sim_P_ma_2nd(t,x,x_p,pars,tss(iter),P_ma_ss,pars_est), ...
               tspan, x0, x_p0);
if     strcmp(change{cc}, 'decrease')
    X(lower_time+1-iter,:) = x(end,:);
elseif strcmp(change{cc}, 'increase')
    X(lower_time+1+iter,:) = x(end,:);
end

end % iteration

pars(end) = '';

end % change

%% Retrieve quantities of interest.

% X = (iteration, variable)
MAP = X(:,41);
GFR = X(:,7 );

end

% =========================================================================
% Solve for steady state baseline value.
% =========================================================================

function SSdata = solve_ss_numerical(pars,pars_est,SSdata_IG)

%% Variables initial guess

% Initial guess for the variables.
% Find the steady state solution, so the derivative is 0.
% Arbitrary value for time to input.
x0 = SSdata_IG; x_p0 = zeros(82,1); t = 0;

%% Find steady state solution

options = optimset('Display','off'); %options = optimset('MaxFunEvals',8200+10000);
[SSdata, ~, ...
 exitflag, output] = fsolve(@(x) bp_reg_solve(t,x,x_p0,pars,pars_est), ...
                            x0, options);

% Check for solver convergence.
if exitflag == 0
    disp('Solver did not converge.')
    disp(output)
    return
end

% Check for imaginary solution.
if not (isreal(SSdata))
    disp('Imaginary number returned.')
    return
end

% Set any values that are within machine precision of 0 equal to 0.
for i = 1:length(SSdata)
    if abs(SSdata(i)) < eps*100
        SSdata(i) = 0;
    end
end

end

% =========================================================================
% Model function for steady state soltuion.
% =========================================================================

function f = bp_reg_solve(t,x,x_p,pars,pars_est)

% pars_est = [x1; x2; y1; y2];
x1 = pars_est(1); x2 = pars_est(2); y1 = pars_est(3); y2 = pars_est(4); 

%% Retrieve parameters by name.

% pars = [N_rsna; R_aass; R_eass; P_B; P_go; C_gcf; eta_etapt; eta_epsdt; ...
%         eta_etacd; K_vd; K_bar; R_bv; T_adh; Phi_sodin; C_K; T_al; ...
%         N_rs; X_PRCPRA; h_renin; h_AGT; h_AngI; h_AngII; h_Ang17; ...
%         h_AngIV; h_AT1R; h_AT2R; c_GPautoreg; P_ghnom; k_AGT; c_ACE; ...
%         c_Chym; c_NEP; c_ACE2; c_IIIV; c_AT1R; c_AT2R; AT1R_eq; ...
%         AT2R_eq; SF; gen];

% Scaling factor
% Rat flow = Human flow x SF
SF = pars(end-1);
pars(end-1) = '';

N_rsna      = pars(1 );
R_aass      = pars(2 );
R_eass      = pars(3 );
P_B         = pars(4 );
P_go        = pars(5 );
C_gcf       = pars(6 );
eta_etapt   = pars(7 );
eta_epsdt   = pars(8 );
eta_etacd   = pars(9 );
K_vd        = pars(10);
K_bar       = pars(11);
R_bv        = pars(12);
T_adh       = pars(13);
Phi_sodin   = pars(14);
C_K         = pars(15);
T_al        = pars(16);
N_rs        = pars(17);
X_PRCPRA    = pars(18);
h_renin     = pars(19);
h_AGT       = pars(20);
h_AngI      = pars(21);
h_AngII     = pars(22);
h_Ang17     = pars(23);
h_AngIV     = pars(24);
h_AT1R      = pars(25);
h_AT2R      = pars(26);
k_AGT       = pars(27);
c_ACE       = pars(28);
c_Chym      = pars(29);
c_NEP       = pars(30);
c_ACE2      = pars(31);
c_IIIV      = pars(32);
c_AT1R      = pars(33);
c_AT2R      = pars(34);
AT1R_eq     = pars(35);
AT2R_eq     = pars(36);
gen         = pars(37);
if     gen == 1
    gender = 'male';
elseif gen == 0
    gender = 'female';
end

%% Retrieve variables by name.

rsna          = x(1 ); rsna_p          = x_p(1 ); 
alpha_map     = x(2 ); alpha_map_p     = x_p(2 ); 
alpha_rap     = x(3 ); alpha_rap_p     = x_p(3 ); 
R_r           = x(4 ); R_r_p           = x_p(4 ); 
beta_rsna     = x(5 ); beta_rsna_p     = x_p(5 ); 
Phi_rb        = x(6 ); Phi_rb_p        = x_p(6 ); 
Phi_gfilt     = x(7 ); Phi_gfilt_p     = x_p(7 ); 
P_f           = x(8 ); P_f_p           = x_p(8 ); 
P_gh          = x(9 ); P_gh_p          = x_p(9 ); 
Sigma_tgf     = x(10); Sigma_tgf_p     = x_p(10); 
Phi_filsod    = x(11); Phi_filsod_p    = x_p(11); 
Phi_ptsodreab = x(12); Phi_ptsodreab_p = x_p(12); 
eta_ptsodreab = x(13); eta_ptsodreab_p = x_p(13); 
gamma_filsod  = x(14); gamma_filsod_p  = x_p(14); 
gamma_at      = x(15); gamma_at_p      = x_p(15); 
gamma_rsna    = x(16); gamma_rsna_p    = x_p(16); 
Phi_mdsod     = x(17); Phi_mdsod_p     = x_p(17); 
Phi_dtsodreab = x(18); Phi_dtsodreab_p = x_p(18); 
eta_dtsodreab = x(19); eta_dtsodreab_p = x_p(19); 
psi_al        = x(20); psi_al_p        = x_p(20); 
Phi_dtsod     = x(21); Phi_dtsod_p     = x_p(21); 
Phi_cdsodreab = x(22); Phi_cdsodreab_p = x_p(22); 
eta_cdsodreab = x(23); eta_cdsodreab_p = x_p(23); 
lambda_dt     = x(24); lambda_dt_p     = x_p(24); 
lambda_anp    = x(25); lambda_anp_p    = x_p(25); 
Phi_usod      = x(26); Phi_usod_p      = x_p(26); 
Phi_win       = x(27); Phi_win_p       = x_p(27); 
V_ecf         = x(28); V_ecf_p         = x_p(28); 
V_b           = x(29); V_b_p           = x_p(29); 
P_mf          = x(30); P_mf_p          = x_p(30); 
Phi_vr        = x(31); Phi_vr_p        = x_p(31); 
Phi_co        = x(32); Phi_co_p        = x_p(32); 
P_ra          = x(33); P_ra_p          = x_p(33); 
vas           = x(34); vas_p           = x_p(34); 
vas_f         = x(35); vas_f_p         = x_p(35); 
vas_d         = x(36); vas_d_p         = x_p(36); 
R_a           = x(37); R_a_p           = x_p(37); 
R_ba          = x(38); R_ba_p          = x_p(38); 
R_vr          = x(39); R_vr_p          = x_p(39); 
R_tp          = x(40); R_tp_p          = x_p(40); 
P_ma          = x(41); P_ma_p          = x_p(41); 
epsilon_aum   = x(42); epsilon_aum_p   = x_p(42); 
a_auto        = x(43); a_auto_p        = x_p(43); 
a_chemo       = x(44); a_chemo_p       = x_p(44); 
a_baro        = x(45); a_baro_p        = x_p(45); 
C_adh         = x(46); C_adh_p         = x_p(46); 
N_adh         = x(47); N_adh_p         = x_p(47); 
N_adhs        = x(48); N_adhs_p        = x_p(48); 
delta_ra      = x(49); delta_ra_p      = x_p(49); 
Phi_twreab    = x(50); Phi_twreab_p    = x_p(50); 
mu_al         = x(51); mu_al_p         = x_p(51); 
mu_adh        = x(52); mu_adh_p        = x_p(52); 
Phi_u         = x(53); Phi_u_p         = x_p(53); 
M_sod         = x(54); M_sod_p         = x_p(54); 
C_sod         = x(55); C_sod_p         = x_p(55); 
nu_mdsod      = x(56); nu_mdsod_p      = x_p(56); 
nu_rsna       = x(57); nu_rsna_p       = x_p(57); 
C_al          = x(58); C_al_p          = x_p(58); 
N_al          = x(59); N_al_p          = x_p(59); 
N_als         = x(60); N_als_p         = x_p(60); 
xi_ksod       = x(61); xi_ksod_p       = x_p(61); 
xi_map        = x(62); xi_map_p        = x_p(62); 
xi_at         = x(63); xi_at_p         = x_p(63); 
hatC_anp      = x(64); hatC_anp_p      = x_p(64); 
AGT           = x(65); AGT_p           = x_p(65); 
nu_AT1        = x(66); nu_AT1_p        = x_p(66); 
R_sec         = x(67); R_sec_p         = x_p(67); 
PRC           = x(68); PRC_p           = x_p(68); 
PRA           = x(69); PRA_p           = x_p(69); 
AngI          = x(70); AngI_p          = x_p(70); 
AngII         = x(71); AngII_p         = x_p(71); 
AT1R          = x(72); AT1R_p          = x_p(72); 
AT2R          = x(73); AT2R_p          = x_p(73); 
Ang17         = x(74); Ang17_p         = x_p(74); 
AngIV         = x(75); AngIV_p         = x_p(75); 
R_aa          = x(76); R_aa_p          = x_p(76); 
R_ea          = x(77); R_ea_p          = x_p(77); 
Sigma_myo     = x(78); Sigma_myo_p     = x_p(78); 
Psi_AT1RAA    = x(79); Psi_AT1RAA_p    = x_p(79); 
Psi_AT1REA    = x(80); Psi_AT1REA_p    = x_p(80); 
Psi_AT2RAA    = x(81); Psi_AT2RAA_p    = x_p(81); 
Psi_AT2REA    = x(82); Psi_AT2REA_p    = x_p(82); 

%% Differential algebraic equation system f(t,x,x') = 0.

f = zeros(length(x),1);

% rsna
rsna0 = N_rsna * alpha_map * alpha_rap;
if     strcmp(gender,'male')
    f(1 ) = rsna - rsna0;
elseif strcmp(gender,'female')
    f(1 ) = rsna - rsna0^(1/rsna0);
%     f(1 ) = rsna - rsna0;
end
% alpha_map
f(2 ) = alpha_map - ( 0.5 + 1 / (1 + exp((P_ma - 100) / 15)) );
% alpha_rap
f(3 ) = alpha_rap - ( 1 - 0.008 * P_ra );
% R_r
f(4 ) = R_r - ( R_aa + R_ea );
% beta_rsna
f(5 ) = beta_rsna - ( 2 / (1 + exp(-3.16 * (rsna - 1))) );
% f(5 ) = beta_rsna - ( 1.5 * (rsna - 1) + 1 );
% Phi_rb
f(6 ) = Phi_rb - ( P_ma / R_r );
% Phi_gfilt
f(7 ) = Phi_gfilt - ( P_f * C_gcf );
% P_f
f(8 ) = P_f - ( P_gh - (P_B + P_go) );
% P_gh
f(9 ) = P_gh - ( P_ma - Phi_rb * R_aa );
% Sigma_tgf - rat - female reabsorption
% f(10) = Sigma_tgf - ( 0.3408 + 3.449 / (3.88 + exp((Phi_mdsod - 3.859) / (-0.9617))) );
if     strcmp(gender,'male')
    f(10) = Sigma_tgf - ( 0.3408 + 3.449 / (3.88 + exp((Phi_mdsod - 3.859 * SF) / (-0.9617 * SF) )) );
elseif strcmp(gender,'female')
    f(10) = Sigma_tgf - ( 0.3408 + 3.449 / (3.88 + exp((Phi_mdsod - 3.859 * SF * 2.500) / (-0.9617 * SF * 2.500) )) );
end
% Phi_filsod
f(11) = Phi_filsod - ( Phi_gfilt * C_sod );
% Phi_ptsodreab
f(12) = Phi_ptsodreab - ( Phi_filsod * eta_ptsodreab );
% eta_ptsodreab
f(13) = eta_ptsodreab - ( eta_etapt * gamma_filsod * gamma_at * gamma_rsna );
% gamma_filsod - rat
% f(14) = gamma_filsod - ( 0.85 + 0.3 / (1 + exp((Phi_filsod - 18)/138)) );
f(14) = gamma_filsod - ( 0.85 + 0.3 / (1 + exp((Phi_filsod - 18 * SF)/(138 * SF) )) );
% gamma_at
f(15) = gamma_at - ( 0.95 + 0.12 / (1 + exp(2.6 - 2.342 * (AT1R/AT1R_eq))) );
% gamma_rsna
f(16) = gamma_rsna - ( 0.72 + 0.56 / (1 + exp((1 - rsna) / 2.18)) );
% Phi_mdsod
f(17) = Phi_mdsod - ( Phi_filsod - Phi_ptsodreab );
% Phi_dtsodreab
f(18) = Phi_dtsodreab - ( Phi_mdsod * eta_dtsodreab );
% eta_dtsodreab
f(19) = eta_dtsodreab - ( eta_epsdt * psi_al );
% psi_al - rat
% f(20) = psi_al - ( 0.17 + 0.94 / (1 + exp((0.48 - 1.2 * log10(C_al)) / 0.88)) );
f(20) = psi_al - ( 0.17 + 0.94 / (1 + exp((0.48 - 0.87 * log10(C_al)) / 0.88)) );
% Phi_dtsod
f(21) = Phi_dtsod - ( Phi_mdsod - Phi_dtsodreab );
% Phi_cdsodreab
f(22) = Phi_cdsodreab - ( Phi_dtsod * eta_cdsodreab );
% eta_cdsodreab
f(23) = eta_cdsodreab - ( eta_etacd * lambda_dt * lambda_anp );
% lambda_dt - rat - female reabsorption
% f(24) = lambda_dt - ( 0.82 + 0.39 / (1 + exp((Phi_dtsod - 1.7625) / 0.375)) );
if     strcmp(gender,'male')
    f(24) = lambda_dt - ( 0.82 + 0.39 / (1 + exp((Phi_dtsod - 1.7625 * SF) / (0.375 * SF) )) );
elseif strcmp(gender,'female')
    f(24) = lambda_dt - ( 0.82 + 0.39 / (1 + exp((Phi_dtsod - 1.7625 * SF * 2.504) / (0.375 * SF * 2.504) )) );
end
% lambda_anp
f(25) = lambda_anp - ( -0.1 * hatC_anp + 1.1 );
% Phi_usod
f(26) = Phi_usod - ( Phi_dtsod - Phi_cdsodreab );
% Phi_win - rat
% f(27) = Phi_win - ( 0.003 / (1 + exp(-2.25 * (C_adh - 3.87))) );
f(27) = Phi_win - ( 0.003 * SF / (1 + exp(-2.25 * (C_adh - 3.87))) );
% V_ecf
f(28) = V_ecf_p - ( Phi_win - Phi_u );
% V_b - rat
% f(29) = V_b - ( 4.5479 + 2.4312 / (1 + exp(-(V_ecf - 18.1128) * 0.4744)) );
f(29) = V_b - ( 4.5479 * SF + 2.4312 * SF / (1 + exp(-(V_ecf - 18.1128 * SF) * (0.4744 / SF) )) );
% P_mf - rat
% f(30) = P_mf - ( (7.436 * V_b - 30.18) * epsilon_aum );
f(30) = P_mf - ( ( (7.436 / SF) * V_b - 30.18) * epsilon_aum );
% Phi_vr
f(31) = Phi_vr - ( (P_mf - P_ra) / R_vr );
% Phi_co
f(32) = Phi_co - ( Phi_vr );
% P_ra - rat
% f(33) = P_ra - ( max( 0, 0.2787 * exp(Phi_co * 0.2281) - 0.8256 ) );
f(33) = P_ra - ( max( 0, 0.2787 * exp(Phi_co * 0.2281 / SF) - 0.8256 ) );
% vas
f(34) = vas_p - ( vas_f - vas_d );
% vas_f - rat
% f(35) = vas_f - ( (11.312 * exp(-Phi_co * 0.4799)) / 100000 );
f(35) = vas_f - ( (11.312 * exp(-Phi_co * 0.4799 / SF)) / 100000 );
% vas_d
f(36) = vas_d - ( vas * K_vd );
% R_a
f(37) = R_a - ( R_ba * epsilon_aum );
% R_ba
f(38) = R_ba - ( K_bar / vas );
% R_vr
f(39) = R_vr - ( (8 * R_bv + R_a) / 31 );
% R_tp
f(40) = R_tp - ( R_a + R_bv );
% P_ma
f(41) = P_ma - ( Phi_co * R_tp );
% epsilon_aum
f(42) = epsilon_aum - ( 4/5 * (a_chemo + a_baro) );
% a_auto
f(43) = a_auto - ( 3.0042 * exp(-0.011 * P_ma) );
% a_chemo
f(44) = a_chemo - ( 1/4 * a_auto );
% a_baro
f(45) = a_baro_p - ( 3/4 * (a_auto_p - 0.0000667 * (a_baro - 1)) );
% C_adh
f(46) = C_adh - ( 4 * N_adh );
% N_adh
f(47) = N_adh_p - ( 1/T_adh * (N_adhs - N_adh) );
% N_adhs
f(48) = N_adhs - ( (C_sod - 141 + max( 0, epsilon_aum - 1 ) - delta_ra) / 3 );
% delta_ra
f(49) = delta_ra_p - ( 0.2 * P_ra_p - 0.0007 * delta_ra );
% Phi_twreab - rat
% f(50) = Phi_twreab - ( 0.025 - 0.001 / (mu_al * mu_adh) + 0.8 * Phi_gfilt );
f(50) = Phi_twreab - ( 0.025 * SF - 0.001 * SF / (mu_al * mu_adh) + 0.8 * Phi_gfilt );
% mu_al - rat
% f(51) = mu_al - ( 0.17 + 0.94 / (1 + exp((0.48 - 1.2 * log10(C_al)) / 0.88)) );
f(51) = mu_al - ( 0.17 + 0.94 / (1 + exp((0.48 - 0.87 * log10(C_al)) / 0.88)) );
% mu_adh
f(52) = mu_adh - ( 0.3313 + 0.8 / (1 + exp(0.6 - 3.7 * log10(C_adh))) ); 
% Phi_u - rat
% f(53) = Phi_u - ( max( 0.0003, Phi_gfilt - Phi_twreab ) );
f(53) = Phi_u - ( max( 0.0003 * SF, Phi_gfilt - Phi_twreab ) );
% M_sod
f(54) = M_sod_p - ( Phi_sodin - Phi_usod );
% C_sod
f(55) = C_sod - ( M_sod / V_ecf );
% nu_mdsod - rat - female reabsorption
% if     strcmp(gender,'male')
%     f(56) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.731) / 0.6056)) );
% elseif strcmp(gender,'female')
%     f(56) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.637) / 0.6056)) );
% end
% % f(56) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.667) / 0.6056)) );
if     strcmp(gender,'male')
    f(56) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.731 * SF) / (0.6056 * SF) )) );
elseif strcmp(gender,'female')
    f(56) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.637 * SF * 2.500) / (0.6056 * SF * 2.500) )) );
end
% nu_rsna
f(57) = nu_rsna - ( 1.822 - 2.056 / (1.358 + exp(rsna - 0.8667)) );
% C_al - rat
% if     strcmp(gender,'male')
%     f(58) = C_al - ( max( 1, N_al * 85      ) );
% elseif strcmp(gender,'female')
%     f(58) = C_al - ( max( 1, N_al * 69.1775 ) );
% end
if     strcmp(gender,'male')
    f(58) = C_al - ( max( 1, N_al * 395.3 ) );
elseif strcmp(gender,'female')
    f(58) = C_al - ( max( 1, N_al * 379.4 ) );
end
% N_al
f(59) = N_al_p - ( 1/T_al * (N_als - N_al) );
% N_als
f(60) = N_als - ( xi_ksod * xi_map * xi_at );
% xi_ksod
f(61) = xi_ksod - ( 5 / ( 1 + exp(0.265 * (C_sod/C_K - 23.6)) ) ); 
% xi_map
if P_ma <= 100
    f(62) = xi_map - ( 70.1054 * exp(-0.0425 * P_ma) );
else
    f(62) = xi_map - ( 1 );
end
% xi_at
f(63) = xi_at - ( 0.47 + 2.4 / (1 + exp((2.82 - 1.952 * (AT1R/AT1R_eq)) / 0.8)) );
% hatC_anp
f(64) = hatC_anp - ( 7.4052 - 6.554 / (1 + exp(P_ra - 3.762)) ); 
% AGT
f(65) = AGT_p - ( k_AGT - PRA - log(2)/h_AGT * AGT );
% nu_AT1
f(66) = nu_AT1 - ( 10^(0.0102 - 0.95 * log10(AT1R / AT1R_eq)) );
% R_sec
f(67) = R_sec - ( N_rs * nu_mdsod * nu_rsna * nu_AT1 );
% PRC
f(68) = PRC_p - ( R_sec - log(2)/h_renin * PRC );
% PRA
f(69) = PRA - ( PRC * X_PRCPRA );
% AngI
f(70) = AngI_p - ( PRA - (c_ACE + c_Chym + c_NEP) * AngI - log(2)/h_AngI * AngI );
% AngII
f(71) = AngII_p - ( (c_ACE + c_Chym) * AngI - (c_ACE2 + c_IIIV + c_AT1R + c_AT2R) * AngII - log(2)/h_AngII * AngII );
% AT1R
f(72) = AT1R_p - ( c_AT1R * AngII - log(2)/h_AT1R * AT1R );
% AT2R
f(73) = AT2R_p - ( c_AT2R * AngII - log(2)/h_AT2R * AT2R );
% Ang17
f(74) = Ang17_p - ( c_NEP * AngI + c_ACE2 * AngII - log(2)/h_Ang17 * Ang17 );
% AngIV
f(75) = AngIV_p - ( c_IIIV * AngII - log(2)/h_AngIV * AngIV );
% R_aa
f(76) = R_aa - ( R_aass * beta_rsna * Sigma_tgf * Sigma_myo * Psi_AT1RAA * Psi_AT2RAA);
% R_ea
f(77) = R_ea - ( R_eass * Psi_AT1REA * Psi_AT2REA );
% ESTIMATE ----------------------------------------------------------------
% Sigma_myo
if     P_gh <= x1
    f(78) = Sigma_myo - ( y1 );
elseif P_gh > x1 && P_gh < x2
    f(78) = Sigma_myo - ( (y2 - y1)/(x2 - x1)*P_gh + (x2*y1 - x1*y2)/(x2 - x1) );
elseif P_gh >= x2
    f(78) = Sigma_myo - ( y2 );
end
% ESTIMATE ----------------------------------------------------------------
% Psi_AT1RAA
f(79) = Psi_AT1RAA - ( 0.8   + 0.2092 * (AT1R / AT1R_eq) - 0.00925 / (AT1R / AT1R_eq) );
% Psi_AT1REA
f(80) = Psi_AT1REA - ( 0.925 + 0.0835 * (AT1R / AT1R_eq) - 0.0085  / (AT1R / AT1R_eq) );
% Psi_AT2RAA
if     strcmp(gender,'male')
    f(81) = Psi_AT2RAA - ( 1 );
elseif strcmp(gender,'female')
%     f(81) = Psi_AT2RAA - ( 0.025 * (AT2R_eq - AT2R) + 1 );
    f(81) = Psi_AT2RAA - ( 0.9 + 0.1 * exp(-(AT2R/AT2R_eq - 1)) );
%     f(81) = Psi_AT2RAA - ( 1 );
end
% Psi_AT2REA
if     strcmp(gender,'male')
    f(82) = Psi_AT2REA - ( 1 );
elseif strcmp(gender,'female')
%     f(82) = Psi_AT2REA - ( 0.01  * (AT2R_eq - AT2R) + 1 );
    f(82) = Psi_AT2REA - ( 0.9 + 0.1 * exp(-(AT2R/AT2R_eq - 1)) );
%     f(82) = Psi_AT2REA - ( 1 );
end

end

% =========================================================================
% Model function for decreasing MAP.
% =========================================================================

function f = bp_reg_sim_P_ma_1st(t,x,x_p,pars,pars_est)

% pars_est = [x1; x2; y1; y2];
x1 = pars_est(1); x2 = pars_est(2); y1 = pars_est(3); y2 = pars_est(4); 

%% Retrieve parameters by name.

% pars = [N_rsna; R_aass; R_eass; P_B; P_go; C_gcf; eta_etapt; eta_epsdt; ...
%         eta_etacd; K_vd; K_bar; R_bv; T_adh; Phi_sodin; C_K; T_al; ...
%         N_rs; X_PRCPRA; h_renin; h_AGT; h_AngI; h_AngII; h_Ang17; ...
%         h_AngIV; h_AT1R; h_AT2R; c_GPautoreg; P_ghnom; k_AGT; c_ACE; ...
%         c_Chym; c_NEP; c_ACE2; c_IIIV; c_AT1R; c_AT2R; AT1R_eq; ...
%         AT2R_eq; SF; gen; ch];

% Scaling factor
% Rat flow = Human flow x SF
SF = pars(end-2);
pars(end-2) = '';

N_rsna      = pars(1 );
R_aass      = pars(2 );
R_eass      = pars(3 );
P_B         = pars(4 );
P_go        = pars(5 );
C_gcf       = pars(6 );
eta_etapt   = pars(7 );
eta_epsdt   = pars(8 );
eta_etacd   = pars(9 );
K_vd        = pars(10);
K_bar       = pars(11);
R_bv        = pars(12);
T_adh       = pars(13);
Phi_sodin   = pars(14);
C_K         = pars(15);
T_al        = pars(16);
N_rs        = pars(17);
X_PRCPRA    = pars(18);
h_renin     = pars(19);
h_AGT       = pars(20);
h_AngI      = pars(21);
h_AngII     = pars(22);
h_Ang17     = pars(23);
h_AngIV     = pars(24);
h_AT1R      = pars(25);
h_AT2R      = pars(26);
k_AGT       = pars(27);
c_ACE       = pars(28);
c_Chym      = pars(29);
c_NEP       = pars(30);
c_ACE2      = pars(31);
c_IIIV      = pars(32);
c_AT1R      = pars(33);
c_AT2R      = pars(34);
AT1R_eq     = pars(35);
AT2R_eq     = pars(36);
gen         = pars(37);
if     gen == 1
    gender = 'male';
elseif gen == 0
    gender = 'female';
end
ch          = pars(38);
if     ch == 1
    change = 'decrease';
elseif ch == 0
    change = 'increase';
end

%% Retrieve variables by name.

rsna          = x(1 ); rsna_p          = x_p(1 ); 
alpha_map     = x(2 ); alpha_map_p     = x_p(2 ); 
alpha_rap     = x(3 ); alpha_rap_p     = x_p(3 ); 
R_r           = x(4 ); R_r_p           = x_p(4 ); 
beta_rsna     = x(5 ); beta_rsna_p     = x_p(5 ); 
Phi_rb        = x(6 ); Phi_rb_p        = x_p(6 ); 
Phi_gfilt     = x(7 ); Phi_gfilt_p     = x_p(7 ); 
P_f           = x(8 ); P_f_p           = x_p(8 ); 
P_gh          = x(9 ); P_gh_p          = x_p(9 ); 
Sigma_tgf     = x(10); Sigma_tgf_p     = x_p(10); 
Phi_filsod    = x(11); Phi_filsod_p    = x_p(11); 
Phi_ptsodreab = x(12); Phi_ptsodreab_p = x_p(12); 
eta_ptsodreab = x(13); eta_ptsodreab_p = x_p(13); 
gamma_filsod  = x(14); gamma_filsod_p  = x_p(14); 
gamma_at      = x(15); gamma_at_p      = x_p(15); 
gamma_rsna    = x(16); gamma_rsna_p    = x_p(16); 
Phi_mdsod     = x(17); Phi_mdsod_p     = x_p(17); 
Phi_dtsodreab = x(18); Phi_dtsodreab_p = x_p(18); 
eta_dtsodreab = x(19); eta_dtsodreab_p = x_p(19); 
psi_al        = x(20); psi_al_p        = x_p(20); 
Phi_dtsod     = x(21); Phi_dtsod_p     = x_p(21); 
Phi_cdsodreab = x(22); Phi_cdsodreab_p = x_p(22); 
eta_cdsodreab = x(23); eta_cdsodreab_p = x_p(23); 
lambda_dt     = x(24); lambda_dt_p     = x_p(24); 
lambda_anp    = x(25); lambda_anp_p    = x_p(25); 
Phi_usod      = x(26); Phi_usod_p      = x_p(26); 
Phi_win       = x(27); Phi_win_p       = x_p(27); 
V_ecf         = x(28); V_ecf_p         = x_p(28); 
V_b           = x(29); V_b_p           = x_p(29); 
P_mf          = x(30); P_mf_p          = x_p(30); 
Phi_vr        = x(31); Phi_vr_p        = x_p(31); 
Phi_co        = x(32); Phi_co_p        = x_p(32); 
P_ra          = x(33); P_ra_p          = x_p(33); 
vas           = x(34); vas_p           = x_p(34); 
vas_f         = x(35); vas_f_p         = x_p(35); 
vas_d         = x(36); vas_d_p         = x_p(36); 
R_a           = x(37); R_a_p           = x_p(37); 
R_ba          = x(38); R_ba_p          = x_p(38); 
R_vr          = x(39); R_vr_p          = x_p(39); 
R_tp          = x(40); R_tp_p          = x_p(40); 
P_ma          = x(41); P_ma_p          = x_p(41); 
epsilon_aum   = x(42); epsilon_aum_p   = x_p(42); 
a_auto        = x(43); a_auto_p        = x_p(43); 
a_chemo       = x(44); a_chemo_p       = x_p(44); 
a_baro        = x(45); a_baro_p        = x_p(45); 
C_adh         = x(46); C_adh_p         = x_p(46); 
N_adh         = x(47); N_adh_p         = x_p(47); 
N_adhs        = x(48); N_adhs_p        = x_p(48); 
delta_ra      = x(49); delta_ra_p      = x_p(49); 
Phi_twreab    = x(50); Phi_twreab_p    = x_p(50); 
mu_al         = x(51); mu_al_p         = x_p(51); 
mu_adh        = x(52); mu_adh_p        = x_p(52); 
Phi_u         = x(53); Phi_u_p         = x_p(53); 
M_sod         = x(54); M_sod_p         = x_p(54); 
C_sod         = x(55); C_sod_p         = x_p(55); 
nu_mdsod      = x(56); nu_mdsod_p      = x_p(56); 
nu_rsna       = x(57); nu_rsna_p       = x_p(57); 
C_al          = x(58); C_al_p          = x_p(58); 
N_al          = x(59); N_al_p          = x_p(59); 
N_als         = x(60); N_als_p         = x_p(60); 
xi_ksod       = x(61); xi_ksod_p       = x_p(61); 
xi_map        = x(62); xi_map_p        = x_p(62); 
xi_at         = x(63); xi_at_p         = x_p(63); 
hatC_anp      = x(64); hatC_anp_p      = x_p(64); 
AGT           = x(65); AGT_p           = x_p(65); 
nu_AT1        = x(66); nu_AT1_p        = x_p(66); 
R_sec         = x(67); R_sec_p         = x_p(67); 
PRC           = x(68); PRC_p           = x_p(68); 
PRA           = x(69); PRA_p           = x_p(69); 
AngI          = x(70); AngI_p          = x_p(70); 
AngII         = x(71); AngII_p         = x_p(71); 
AT1R          = x(72); AT1R_p          = x_p(72); 
AT2R          = x(73); AT2R_p          = x_p(73); 
Ang17         = x(74); Ang17_p         = x_p(74); 
AngIV         = x(75); AngIV_p         = x_p(75); 
R_aa          = x(76); R_aa_p          = x_p(76); 
R_ea          = x(77); R_ea_p          = x_p(77); 
Sigma_myo     = x(78); Sigma_myo_p     = x_p(78); 
Psi_AT1RAA    = x(79); Psi_AT1RAA_p    = x_p(79); 
Psi_AT1REA    = x(80); Psi_AT1REA_p    = x_p(80); 
Psi_AT2RAA    = x(81); Psi_AT2RAA_p    = x_p(81); 
Psi_AT2REA    = x(82); Psi_AT2REA_p    = x_p(82); 

%% Differential algebraic equation system f(t,x,x') = 0.

f = zeros(length(x),1);

% rsna
rsna0 = N_rsna * alpha_map * alpha_rap;
if     strcmp(gender,'male')
    f(1 ) = rsna - rsna0;
elseif strcmp(gender,'female')
    f(1 ) = rsna - rsna0^(1/rsna0);
%     f(1 ) = rsna - rsna0;
end
% alpha_map
f(2 ) = alpha_map - ( 0.5 + 1 / (1 + exp((P_ma - 100) / 15)) );
% alpha_rap
f(3 ) = alpha_rap - ( 1 - 0.008 * P_ra );
% R_r
f(4 ) = R_r - ( R_aa + R_ea );
% beta_rsna
f(5 ) = beta_rsna - ( 2 / (1 + exp(-3.16 * (rsna - 1))) );
% f(5 ) = beta_rsna - ( 1.5 * (rsna - 1) + 1 );
% Phi_rb
f(6 ) = Phi_rb - ( P_ma / R_r );
% Phi_gfilt
f(7 ) = Phi_gfilt - ( P_f * C_gcf );
% P_f
f(8 ) = P_f - ( P_gh - (P_B + P_go) );
% P_gh
f(9 ) = P_gh - ( P_ma - Phi_rb * R_aa );
% Sigma_tgf - rat - female reabsorption
% f(10) = Sigma_tgf - ( 0.3408 + 3.449 / (3.88 + exp((Phi_mdsod - 3.859) / (-0.9617))) );
if     strcmp(gender,'male')
    f(10) = Sigma_tgf - ( 0.3408 + 3.449 / (3.88 + exp((Phi_mdsod - 3.859 * SF) / (-0.9617 * SF) )) );
elseif strcmp(gender,'female')
    f(10) = Sigma_tgf - ( 0.3408 + 3.449 / (3.88 + exp((Phi_mdsod - 3.859 * SF * 2.500) / (-0.9617 * SF * 2.500) )) );
end
% Phi_filsod
f(11) = Phi_filsod - ( Phi_gfilt * C_sod );
% Phi_ptsodreab
f(12) = Phi_ptsodreab - ( Phi_filsod * eta_ptsodreab );
% eta_ptsodreab
f(13) = eta_ptsodreab - ( eta_etapt * gamma_filsod * gamma_at * gamma_rsna );
% gamma_filsod - rat
% f(14) = gamma_filsod - ( 0.85 + 0.3 / (1 + exp((Phi_filsod - 18)/138)) );
f(14) = gamma_filsod - ( 0.85 + 0.3 / (1 + exp((Phi_filsod - 18 * SF)/(138 * SF) )) );
% gamma_at
f(15) = gamma_at - ( 0.95 + 0.12 / (1 + exp(2.6 - 2.342 * (AT1R/AT1R_eq))) );
% gamma_rsna
f(16) = gamma_rsna - ( 0.72 + 0.56 / (1 + exp((1 - rsna) / 2.18)) );
% Phi_mdsod
f(17) = Phi_mdsod - ( Phi_filsod - Phi_ptsodreab );
% Phi_dtsodreab
f(18) = Phi_dtsodreab - ( Phi_mdsod * eta_dtsodreab );
% eta_dtsodreab
f(19) = eta_dtsodreab - ( eta_epsdt * psi_al );
% psi_al - rat
% f(20) = psi_al - ( 0.17 + 0.94 / (1 + exp((0.48 - 1.2 * log10(C_al)) / 0.88)) );
f(20) = psi_al - ( 0.17 + 0.94 / (1 + exp((0.48 - 0.87 * log10(C_al)) / 0.88)) );
% Phi_dtsod
f(21) = Phi_dtsod - ( Phi_mdsod - Phi_dtsodreab );
% Phi_cdsodreab
f(22) = Phi_cdsodreab - ( Phi_dtsod * eta_cdsodreab );
% eta_cdsodreab
f(23) = eta_cdsodreab - ( eta_etacd * lambda_dt * lambda_anp );
% lambda_dt - rat - female reabsorption
% f(24) = lambda_dt - ( 0.82 + 0.39 / (1 + exp((Phi_dtsod - 1.7625) / 0.375)) );
if     strcmp(gender,'male')
    f(24) = lambda_dt - ( 0.82 + 0.39 / (1 + exp((Phi_dtsod - 1.7625 * SF) / (0.375 * SF) )) );
elseif strcmp(gender,'female')
    f(24) = lambda_dt - ( 0.82 + 0.39 / (1 + exp((Phi_dtsod - 1.7625 * SF * 2.504) / (0.375 * SF * 2.504) )) );
end
% lambda_anp
f(25) = lambda_anp - ( -0.1 * hatC_anp + 1.1 );
% Phi_usod
f(26) = Phi_usod - ( Phi_dtsod - Phi_cdsodreab );
% Phi_win - rat
% f(27) = Phi_win - ( 0.003 / (1 + exp(-2.25 * (C_adh - 3.87))) );
f(27) = Phi_win - ( 0.003 * SF / (1 + exp(-2.25 * (C_adh - 3.87))) );
% V_ecf
f(28) = V_ecf_p - ( Phi_win - Phi_u );
% V_b - rat
% f(29) = V_b - ( 4.5479 + 2.4312 / (1 + exp(-(V_ecf - 18.1128) * 0.4744)) );
f(29) = V_b - ( 4.5479 * SF + 2.4312 * SF / (1 + exp(-(V_ecf - 18.1128 * SF) * (0.4744 / SF) )) );
% P_mf - rat
% f(30) = P_mf - ( (7.436 * V_b - 30.18) * epsilon_aum );
f(30) = P_mf - ( ( (7.436 / SF) * V_b - 30.18) * epsilon_aum );
% Phi_vr
f(31) = Phi_vr - ( (P_mf - P_ra) / R_vr );
% Phi_co
f(32) = Phi_co - ( Phi_vr );
% P_ra - rat
% f(33) = P_ra - ( max( 0, 0.2787 * exp(Phi_co * 0.2281) - 0.8256 ) );
f(33) = P_ra - ( max( 0, 0.2787 * exp(Phi_co * 0.2281 / SF) - 0.8256 ) );
% vas
f(34) = vas_p - ( vas_f - vas_d );
% vas_f - rat
% f(35) = vas_f - ( (11.312 * exp(-Phi_co * 0.4799)) / 100000 );
f(35) = vas_f - ( (11.312 * exp(-Phi_co * 0.4799 / SF)) / 100000 );
% vas_d
f(36) = vas_d - ( vas * K_vd );
% R_a
f(37) = R_a - ( R_ba * epsilon_aum );
% R_ba
f(38) = R_ba - ( K_bar / vas );
% R_vr
f(39) = R_vr - ( (8 * R_bv + R_a) / 31 );
% R_tp
f(40) = R_tp - ( R_a + R_bv );
% VARY --------------------------------------------------------------------
% P_ma
if     strcmp(change, 'decrease')
    f(41) = P_ma - ( Phi_co * R_tp - t );
elseif strcmp(change, 'increase')
    f(41) = P_ma - ( Phi_co * R_tp + t );
end
% VARY --------------------------------------------------------------------
% epsilon_aum
f(42) = epsilon_aum - ( 4/5 * (a_chemo + a_baro) ); 
% a_auto
f(43) = a_auto - ( 3.0042 * exp(-0.011 * P_ma) );
% a_chemo
f(44) = a_chemo - ( 1/4 * a_auto );
% a_baro
f(45) = a_baro_p - ( 3/4 * (a_auto_p - 0.0000667 * (a_baro - 1)) );
% C_adh
f(46) = C_adh - ( 4 * N_adh );
% N_adh
f(47) = N_adh_p - ( 1/T_adh * (N_adhs - N_adh) );
% N_adhs
f(48) = N_adhs - ( (C_sod - 141 + max( 0, epsilon_aum - 1 ) - delta_ra) / 3 );
% delta_ra
f(49) = delta_ra_p - ( 0.2 * P_ra_p - 0.0007 * delta_ra );
% Phi_twreab - rat
% f(50) = Phi_twreab - ( 0.025 - 0.001 / (mu_al * mu_adh) + 0.8 * Phi_gfilt );
f(50) = Phi_twreab - ( 0.025 * SF - 0.001 * SF / (mu_al * mu_adh) + 0.8 * Phi_gfilt );
% mu_al - rat
% f(51) = mu_al - ( 0.17 + 0.94 / (1 + exp((0.48 - 1.2 * log10(C_al)) / 0.88)) );
f(51) = mu_al - ( 0.17 + 0.94 / (1 + exp((0.48 - 0.87 * log10(C_al)) / 0.88)) );
% mu_adh
f(52) = mu_adh - ( 0.3313 + 0.8 / (1 + exp(0.6 - 3.7 * log10(C_adh))) ); 
% Phi_u - rat
% f(53) = Phi_u - ( max( 0.0003, Phi_gfilt - Phi_twreab ) );
f(53) = Phi_u - ( max( 0.0003 * SF, Phi_gfilt - Phi_twreab ) );
% M_sod
f(54) = M_sod_p - ( Phi_sodin - Phi_usod );
% C_sod
f(55) = C_sod - ( M_sod / V_ecf );
% nu_mdsod - rat - female reabsorption
% if     strcmp(gender,'male')
%     f(56) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.731) / 0.6056)) );
% elseif strcmp(gender,'female')
%     f(56) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.637) / 0.6056)) );
% end
% % f(56) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.667) / 0.6056)) );
if     strcmp(gender,'male')
    f(56) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.731 * SF) / (0.6056 * SF) )) );
elseif strcmp(gender,'female')
    f(56) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.637 * SF * 2.500) / (0.6056 * SF * 2.500) )) );
end
% nu_rsna
f(57) = nu_rsna - ( 1.822 - 2.056 / (1.358 + exp(rsna - 0.8667)) );
% C_al - rat
% if     strcmp(gender,'male')
%     f(58) = C_al - ( max( 1, N_al * 85      ) );
% elseif strcmp(gender,'female')
%     f(58) = C_al - ( max( 1, N_al * 69.1775 ) );
% end
if     strcmp(gender,'male')
    f(58) = C_al - ( max( 1, N_al * 395.3 ) );
elseif strcmp(gender,'female')
    f(58) = C_al - ( max( 1, N_al * 379.4 ) );
end
% N_al
f(59) = N_al_p - ( 1/T_al * (N_als - N_al) );
% N_als
f(60) = N_als - ( xi_ksod * xi_map * xi_at );
% xi_ksod
f(61) = xi_ksod - ( 5 / ( 1 + exp(0.265 * (C_sod/C_K - 23.6)) ) ); 
% xi_map
if P_ma <= 100
    f(62) = xi_map - ( 70.1054 * exp(-0.0425 * P_ma) );
else
    f(62) = xi_map - ( 1 );
end
% xi_at
f(63) = xi_at - ( 0.47 + 2.4 / (1 + exp((2.82 - 1.952 * (AT1R/AT1R_eq)) / 0.8)) );
% hatC_anp
f(64) = hatC_anp - ( 7.4052 - 6.554 / (1 + exp(P_ra - 3.762)) ); 
% AGT
f(65) = AGT_p - ( k_AGT - PRA - log(2)/h_AGT * AGT );
% nu_AT1
f(66) = nu_AT1 - ( 10^(0.0102 - 0.95 * log10(AT1R / AT1R_eq)) );
% R_sec
f(67) = R_sec - ( N_rs * nu_mdsod * nu_rsna * nu_AT1 );
% PRC
f(68) = PRC_p - ( R_sec - log(2)/h_renin * PRC );
% PRA
f(69) = PRA - ( PRC * X_PRCPRA );
% AngI
f(70) = AngI_p - ( PRA - (c_ACE + c_Chym + c_NEP) * AngI - log(2)/h_AngI * AngI );
% AngII
f(71) = AngII_p - ( (c_ACE + c_Chym) * AngI - (c_ACE2 + c_IIIV + c_AT1R + c_AT2R) * AngII - log(2)/h_AngII * AngII );
% AT1R
f(72) = AT1R_p - ( c_AT1R * AngII - log(2)/h_AT1R * AT1R );
% AT2R
f(73) = AT2R_p - ( c_AT2R * AngII - log(2)/h_AT2R * AT2R );
% Ang17
f(74) = Ang17_p - ( c_NEP * AngI + c_ACE2 * AngII - log(2)/h_Ang17 * Ang17 );
% AngIV
f(75) = AngIV_p - ( c_IIIV * AngII - log(2)/h_AngIV * AngIV );
% R_aa
f(76) = R_aa - ( R_aass * beta_rsna * Sigma_tgf * Sigma_myo * Psi_AT1RAA * Psi_AT2RAA);
% R_ea
f(77) = R_ea - ( R_eass * Psi_AT1REA * Psi_AT2REA );
% ESTIMATE ----------------------------------------------------------------
% Sigma_myo
if     P_gh <= x1
    f(78) = Sigma_myo - ( y1 );
elseif P_gh > x1 && P_gh < x2
    f(78) = Sigma_myo - ( (y2 - y1)/(x2 - x1)*P_gh + (x2*y1 - x1*y2)/(x2 - x1) );
elseif P_gh >= x2
    f(78) = Sigma_myo - ( y2 );
end
% ESTIMATE ----------------------------------------------------------------
% Psi_AT1RAA
f(79) = Psi_AT1RAA - ( 0.8   + 0.2092 * (AT1R / AT1R_eq) - 0.00925 / (AT1R / AT1R_eq) );
% Psi_AT1REA
f(80) = Psi_AT1REA - ( 0.925 + 0.0835 * (AT1R / AT1R_eq) - 0.0085  / (AT1R / AT1R_eq) );
% Psi_AT2RAA
if     strcmp(gender,'male')
    f(81) = Psi_AT2RAA - ( 1 );
elseif strcmp(gender,'female')
%     f(81) = Psi_AT2RAA - ( 0.025 * (AT2R_eq - AT2R) + 1 );
    f(81) = Psi_AT2RAA - ( 0.9 + 0.1 * exp(-(AT2R/AT2R_eq - 1)) );
%     f(81) = Psi_AT2RAA - ( 1 );
end
% Psi_AT2REA
if     strcmp(gender,'male')
    f(82) = Psi_AT2REA - ( 1 );
elseif strcmp(gender,'female')
%     f(82) = Psi_AT2REA - ( 0.01  * (AT2R_eq - AT2R) + 1 );
    f(82) = Psi_AT2REA - ( 0.9 + 0.1 * exp(-(AT2R/AT2R_eq - 1)) );
%     f(82) = Psi_AT2REA - ( 1 );
end

end

% =========================================================================
% Model function for increasing MAP.
% =========================================================================

function f = bp_reg_sim_P_ma_2nd(t,x,x_p,pars,tss,P_ma_ss,pars_est)

% pars_est = [x1; x2; y1; y2];
x1 = pars_est(1); x2 = pars_est(2); y1 = pars_est(3); y2 = pars_est(4); 

%% Retrieve parameters by name.

% pars = [N_rsna; R_aass; R_eass; P_B; P_go; C_gcf; eta_etapt; eta_epsdt; ...
%         eta_etacd; K_vd; K_bar; R_bv; T_adh; Phi_sodin; C_K; T_al; ...
%         N_rs; X_PRCPRA; h_renin; h_AGT; h_AngI; h_AngII; h_Ang17; ...
%         h_AngIV; h_AT1R; h_AT2R; c_GPautoreg; P_ghnom; k_AGT; c_ACE; ...
%         c_Chym; c_NEP; c_ACE2; c_IIIV; c_AT1R; c_AT2R; AT1R_eq; ...
%         AT2R_eq; SF; gen; ch];

% Scaling factor
% Rat flow = Human flow x SF
SF = pars(end-2);
pars(end-2) = '';

N_rsna      = pars(1 );
R_aass      = pars(2 );
R_eass      = pars(3 );
P_B         = pars(4 );
P_go        = pars(5 );
C_gcf       = pars(6 );
eta_etapt   = pars(7 );
eta_epsdt   = pars(8 );
eta_etacd   = pars(9 );
K_vd        = pars(10);
K_bar       = pars(11);
R_bv        = pars(12);
T_adh       = pars(13);
Phi_sodin   = pars(14);
C_K         = pars(15);
T_al        = pars(16);
N_rs        = pars(17);
X_PRCPRA    = pars(18);
h_renin     = pars(19);
h_AGT       = pars(20);
h_AngI      = pars(21);
h_AngII     = pars(22);
h_Ang17     = pars(23);
h_AngIV     = pars(24);
h_AT1R      = pars(25);
h_AT2R      = pars(26);
k_AGT       = pars(27);
c_ACE       = pars(28);
c_Chym      = pars(29);
c_NEP       = pars(30);
c_ACE2      = pars(31);
c_IIIV      = pars(32);
c_AT1R      = pars(33);
c_AT2R      = pars(34);
AT1R_eq     = pars(35);
AT2R_eq     = pars(36);
gen         = pars(37);
if     gen == 1
    gender = 'male';
elseif gen == 0
    gender = 'female';
end
ch          = pars(38);
if     ch == 1
    change = 'decrease';
elseif ch == 0
    change = 'increase';
end

%% Retrieve variables by name.

rsna          = x(1 ); rsna_p          = x_p(1 ); 
alpha_map     = x(2 ); alpha_map_p     = x_p(2 ); 
alpha_rap     = x(3 ); alpha_rap_p     = x_p(3 ); 
R_r           = x(4 ); R_r_p           = x_p(4 ); 
beta_rsna     = x(5 ); beta_rsna_p     = x_p(5 ); 
Phi_rb        = x(6 ); Phi_rb_p        = x_p(6 ); 
Phi_gfilt     = x(7 ); Phi_gfilt_p     = x_p(7 ); 
P_f           = x(8 ); P_f_p           = x_p(8 ); 
P_gh          = x(9 ); P_gh_p          = x_p(9 ); 
Sigma_tgf     = x(10); Sigma_tgf_p     = x_p(10); 
Phi_filsod    = x(11); Phi_filsod_p    = x_p(11); 
Phi_ptsodreab = x(12); Phi_ptsodreab_p = x_p(12); 
eta_ptsodreab = x(13); eta_ptsodreab_p = x_p(13); 
gamma_filsod  = x(14); gamma_filsod_p  = x_p(14); 
gamma_at      = x(15); gamma_at_p      = x_p(15); 
gamma_rsna    = x(16); gamma_rsna_p    = x_p(16); 
Phi_mdsod     = x(17); Phi_mdsod_p     = x_p(17); 
Phi_dtsodreab = x(18); Phi_dtsodreab_p = x_p(18); 
eta_dtsodreab = x(19); eta_dtsodreab_p = x_p(19); 
psi_al        = x(20); psi_al_p        = x_p(20); 
Phi_dtsod     = x(21); Phi_dtsod_p     = x_p(21); 
Phi_cdsodreab = x(22); Phi_cdsodreab_p = x_p(22); 
eta_cdsodreab = x(23); eta_cdsodreab_p = x_p(23); 
lambda_dt     = x(24); lambda_dt_p     = x_p(24); 
lambda_anp    = x(25); lambda_anp_p    = x_p(25); 
Phi_usod      = x(26); Phi_usod_p      = x_p(26); 
Phi_win       = x(27); Phi_win_p       = x_p(27); 
V_ecf         = x(28); V_ecf_p         = x_p(28); 
V_b           = x(29); V_b_p           = x_p(29); 
P_mf          = x(30); P_mf_p          = x_p(30); 
Phi_vr        = x(31); Phi_vr_p        = x_p(31); 
Phi_co        = x(32); Phi_co_p        = x_p(32); 
P_ra          = x(33); P_ra_p          = x_p(33); 
vas           = x(34); vas_p           = x_p(34); 
vas_f         = x(35); vas_f_p         = x_p(35); 
vas_d         = x(36); vas_d_p         = x_p(36); 
R_a           = x(37); R_a_p           = x_p(37); 
R_ba          = x(38); R_ba_p          = x_p(38); 
R_vr          = x(39); R_vr_p          = x_p(39); 
R_tp          = x(40); R_tp_p          = x_p(40); 
P_ma          = x(41); P_ma_p          = x_p(41); 
epsilon_aum   = x(42); epsilon_aum_p   = x_p(42); 
a_auto        = x(43); a_auto_p        = x_p(43); 
a_chemo       = x(44); a_chemo_p       = x_p(44); 
a_baro        = x(45); a_baro_p        = x_p(45); 
C_adh         = x(46); C_adh_p         = x_p(46); 
N_adh         = x(47); N_adh_p         = x_p(47); 
N_adhs        = x(48); N_adhs_p        = x_p(48); 
delta_ra      = x(49); delta_ra_p      = x_p(49); 
Phi_twreab    = x(50); Phi_twreab_p    = x_p(50); 
mu_al         = x(51); mu_al_p         = x_p(51); 
mu_adh        = x(52); mu_adh_p        = x_p(52); 
Phi_u         = x(53); Phi_u_p         = x_p(53); 
M_sod         = x(54); M_sod_p         = x_p(54); 
C_sod         = x(55); C_sod_p         = x_p(55); 
nu_mdsod      = x(56); nu_mdsod_p      = x_p(56); 
nu_rsna       = x(57); nu_rsna_p       = x_p(57); 
C_al          = x(58); C_al_p          = x_p(58); 
N_al          = x(59); N_al_p          = x_p(59); 
N_als         = x(60); N_als_p         = x_p(60); 
xi_ksod       = x(61); xi_ksod_p       = x_p(61); 
xi_map        = x(62); xi_map_p        = x_p(62); 
xi_at         = x(63); xi_at_p         = x_p(63); 
hatC_anp      = x(64); hatC_anp_p      = x_p(64); 
AGT           = x(65); AGT_p           = x_p(65); 
nu_AT1        = x(66); nu_AT1_p        = x_p(66); 
R_sec         = x(67); R_sec_p         = x_p(67); 
PRC           = x(68); PRC_p           = x_p(68); 
PRA           = x(69); PRA_p           = x_p(69); 
AngI          = x(70); AngI_p          = x_p(70); 
AngII         = x(71); AngII_p         = x_p(71); 
AT1R          = x(72); AT1R_p          = x_p(72); 
AT2R          = x(73); AT2R_p          = x_p(73); 
Ang17         = x(74); Ang17_p         = x_p(74); 
AngIV         = x(75); AngIV_p         = x_p(75); 
R_aa          = x(76); R_aa_p          = x_p(76); 
R_ea          = x(77); R_ea_p          = x_p(77); 
Sigma_myo     = x(78); Sigma_myo_p     = x_p(78); 
Psi_AT1RAA    = x(79); Psi_AT1RAA_p    = x_p(79); 
Psi_AT1REA    = x(80); Psi_AT1REA_p    = x_p(80); 
Psi_AT2RAA    = x(81); Psi_AT2RAA_p    = x_p(81); 
Psi_AT2REA    = x(82); Psi_AT2REA_p    = x_p(82); 

%% Differential algebraic equation system f(t,x,x') = 0.

f = zeros(length(x),1);

% rsna
rsna0 = N_rsna * alpha_map * alpha_rap;
if     strcmp(gender,'male')
    f(1 ) = rsna - rsna0;
elseif strcmp(gender,'female')
    f(1 ) = rsna - rsna0^(1/rsna0);
%     f(1 ) = rsna - rsna0;
end
% alpha_map
f(2 ) = alpha_map - ( 0.5 + 1 / (1 + exp((P_ma - 100) / 15)) );
% alpha_rap
f(3 ) = alpha_rap - ( 1 - 0.008 * P_ra );
% R_r
f(4 ) = R_r - ( R_aa + R_ea );
% beta_rsna
f(5 ) = beta_rsna - ( 2 / (1 + exp(-3.16 * (rsna - 1))) );
% f(5 ) = beta_rsna - ( 1.5 * (rsna - 1) + 1 );
% Phi_rb
f(6 ) = Phi_rb - ( P_ma / R_r );
% Phi_gfilt
f(7 ) = Phi_gfilt - ( P_f * C_gcf );
% P_f
f(8 ) = P_f - ( P_gh - (P_B + P_go) );
% P_gh
f(9 ) = P_gh - ( P_ma - Phi_rb * R_aa );
% Sigma_tgf - rat - female reabsorption
% f(10) = Sigma_tgf - ( 0.3408 + 3.449 / (3.88 + exp((Phi_mdsod - 3.859) / (-0.9617))) );
if     strcmp(gender,'male')
    f(10) = Sigma_tgf - ( 0.3408 + 3.449 / (3.88 + exp((Phi_mdsod - 3.859 * SF) / (-0.9617 * SF) )) );
elseif strcmp(gender,'female')
    f(10) = Sigma_tgf - ( 0.3408 + 3.449 / (3.88 + exp((Phi_mdsod - 3.859 * SF * 2.500) / (-0.9617 * SF * 2.500) )) );
end
% Phi_filsod
f(11) = Phi_filsod - ( Phi_gfilt * C_sod );
% Phi_ptsodreab
f(12) = Phi_ptsodreab - ( Phi_filsod * eta_ptsodreab );
% eta_ptsodreab
f(13) = eta_ptsodreab - ( eta_etapt * gamma_filsod * gamma_at * gamma_rsna );
% gamma_filsod - rat
% f(14) = gamma_filsod - ( 0.85 + 0.3 / (1 + exp((Phi_filsod - 18)/138)) );
f(14) = gamma_filsod - ( 0.85 + 0.3 / (1 + exp((Phi_filsod - 18 * SF)/(138 * SF) )) );
% gamma_at
f(15) = gamma_at - ( 0.95 + 0.12 / (1 + exp(2.6 - 2.342 * (AT1R/AT1R_eq))) );
% gamma_rsna
f(16) = gamma_rsna - ( 0.72 + 0.56 / (1 + exp((1 - rsna) / 2.18)) );
% Phi_mdsod
f(17) = Phi_mdsod - ( Phi_filsod - Phi_ptsodreab );
% Phi_dtsodreab
f(18) = Phi_dtsodreab - ( Phi_mdsod * eta_dtsodreab );
% eta_dtsodreab
f(19) = eta_dtsodreab - ( eta_epsdt * psi_al );
% psi_al - rat
% f(20) = psi_al - ( 0.17 + 0.94 / (1 + exp((0.48 - 1.2 * log10(C_al)) / 0.88)) );
f(20) = psi_al - ( 0.17 + 0.94 / (1 + exp((0.48 - 0.87 * log10(C_al)) / 0.88)) );
% Phi_dtsod
f(21) = Phi_dtsod - ( Phi_mdsod - Phi_dtsodreab );
% Phi_cdsodreab
f(22) = Phi_cdsodreab - ( Phi_dtsod * eta_cdsodreab );
% eta_cdsodreab
f(23) = eta_cdsodreab - ( eta_etacd * lambda_dt * lambda_anp );
% lambda_dt - rat - female reabsorption
% f(24) = lambda_dt - ( 0.82 + 0.39 / (1 + exp((Phi_dtsod - 1.7625) / 0.375)) );
if     strcmp(gender,'male')
    f(24) = lambda_dt - ( 0.82 + 0.39 / (1 + exp((Phi_dtsod - 1.7625 * SF) / (0.375 * SF) )) );
elseif strcmp(gender,'female')
    f(24) = lambda_dt - ( 0.82 + 0.39 / (1 + exp((Phi_dtsod - 1.7625 * SF * 2.504) / (0.375 * SF * 2.504) )) );
end
% lambda_anp
f(25) = lambda_anp - ( -0.1 * hatC_anp + 1.1 );
% Phi_usod
f(26) = Phi_usod - ( Phi_dtsod - Phi_cdsodreab );
% Phi_win - rat
% f(27) = Phi_win - ( 0.003 / (1 + exp(-2.25 * (C_adh - 3.87))) );
f(27) = Phi_win - ( 0.003 * SF / (1 + exp(-2.25 * (C_adh - 3.87))) );
% V_ecf
f(28) = V_ecf_p - ( Phi_win - Phi_u );
% V_b - rat
% f(29) = V_b - ( 4.5479 + 2.4312 / (1 + exp(-(V_ecf - 18.1128) * 0.4744)) );
f(29) = V_b - ( 4.5479 * SF + 2.4312 * SF / (1 + exp(-(V_ecf - 18.1128 * SF) * (0.4744 / SF) )) );
% P_mf - rat
% f(30) = P_mf - ( (7.436 * V_b - 30.18) * epsilon_aum );
f(30) = P_mf - ( ( (7.436 / SF) * V_b - 30.18) * epsilon_aum );
% Phi_vr
f(31) = Phi_vr - ( (P_mf - P_ra) / R_vr );
% Phi_co
f(32) = Phi_co - ( Phi_vr );
% P_ra - rat
% f(33) = P_ra - ( max( 0, 0.2787 * exp(Phi_co * 0.2281) - 0.8256 ) );
f(33) = P_ra - ( max( 0, 0.2787 * exp(Phi_co * 0.2281 / SF) - 0.8256 ) );
% vas
f(34) = vas_p - ( vas_f - vas_d );
% vas_f - rat
% f(35) = vas_f - ( (11.312 * exp(-Phi_co * 0.4799)) / 100000 );
f(35) = vas_f - ( (11.312 * exp(-Phi_co * 0.4799 / SF)) / 100000 );
% vas_d
f(36) = vas_d - ( vas * K_vd );
% R_a
f(37) = R_a - ( R_ba * epsilon_aum );
% R_ba
f(38) = R_ba - ( K_bar / vas );
% R_vr
f(39) = R_vr - ( (8 * R_bv + R_a) / 31 );
% R_tp
f(40) = R_tp - ( R_a + R_bv );
% VARY --------------------------------------------------------------------
% P_ma
% f(41) = P_ma - ( Phi_co * R_tp );
% Test varying blood pressure.
if     t <= tss
    if     strcmp(change, 'decrease')
        f(41) = P_ma - ( Phi_co * R_tp - t );
    elseif strcmp(change, 'increase')
        f(41) = P_ma - ( Phi_co * R_tp + t );
    end
elseif t > tss 
    f(41) = P_ma - ( P_ma_ss );
end
% VARY --------------------------------------------------------------------
% epsilon_aum
f(42) = epsilon_aum - ( 4/5 * (a_chemo + a_baro) );
% a_auto
f(43) = a_auto - ( 3.0042 * exp(-0.011 * P_ma) );
% a_chemo
f(44) = a_chemo - ( 1/4 * a_auto );
% a_baro
f(45) = a_baro_p - ( 3/4 * (a_auto_p - 0.0000667 * (a_baro - 1)) );
% C_adh
f(46) = C_adh - ( 4 * N_adh );
% N_adh
f(47) = N_adh_p - ( 1/T_adh * (N_adhs - N_adh) );
% N_adhs
f(48) = N_adhs - ( (C_sod - 141 + max( 0, epsilon_aum - 1 ) - delta_ra) / 3 );
% delta_ra
f(49) = delta_ra_p - ( 0.2 * P_ra_p - 0.0007 * delta_ra );
% Phi_twreab - rat
% f(50) = Phi_twreab - ( 0.025 - 0.001 / (mu_al * mu_adh) + 0.8 * Phi_gfilt );
f(50) = Phi_twreab - ( 0.025 * SF - 0.001 * SF / (mu_al * mu_adh) + 0.8 * Phi_gfilt );
% mu_al - rat
% f(51) = mu_al - ( 0.17 + 0.94 / (1 + exp((0.48 - 1.2 * log10(C_al)) / 0.88)) );
f(51) = mu_al - ( 0.17 + 0.94 / (1 + exp((0.48 - 0.87 * log10(C_al)) / 0.88)) );
% mu_adh
f(52) = mu_adh - ( 0.3313 + 0.8 / (1 + exp(0.6 - 3.7 * log10(C_adh))) ); 
% Phi_u - rat
% f(53) = Phi_u - ( max( 0.0003, Phi_gfilt - Phi_twreab ) );
f(53) = Phi_u - ( max( 0.0003 * SF, Phi_gfilt - Phi_twreab ) );
% M_sod
f(54) = M_sod_p - ( Phi_sodin - Phi_usod );
% C_sod
f(55) = C_sod - ( M_sod / V_ecf );
% nu_mdsod - rat - female reabsorption
% if     strcmp(gender,'male')
%     f(56) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.731) / 0.6056)) );
% elseif strcmp(gender,'female')
%     f(56) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.637) / 0.6056)) );
% end
% % f(56) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.667) / 0.6056)) );
if     strcmp(gender,'male')
    f(56) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.731 * SF) / (0.6056 * SF) )) );
elseif strcmp(gender,'female')
    f(56) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.637 * SF * 2.500) / (0.6056 * SF * 2.500) )) );
end
% nu_rsna
f(57) = nu_rsna - ( 1.822 - 2.056 / (1.358 + exp(rsna - 0.8667)) );
% C_al - rat
% if     strcmp(gender,'male')
%     f(58) = C_al - ( max( 1, N_al * 85      ) );
% elseif strcmp(gender,'female')
%     f(58) = C_al - ( max( 1, N_al * 69.1775 ) );
% end
if     strcmp(gender,'male')
    f(58) = C_al - ( max( 1, N_al * 395.3 ) );
elseif strcmp(gender,'female')
    f(58) = C_al - ( max( 1, N_al * 379.4 ) );
end
% N_al
f(59) = N_al_p - ( 1/T_al * (N_als - N_al) );
% N_als
f(60) = N_als - ( xi_ksod * xi_map * xi_at );
% xi_ksod
f(61) = xi_ksod - ( 5 / ( 1 + exp(0.265 * (C_sod/C_K - 23.6)) ) ); 
% xi_map
if P_ma <= 100
    f(62) = xi_map - ( 70.1054 * exp(-0.0425 * P_ma) );
else
    f(62) = xi_map - ( 1 );
end
% xi_at
f(63) = xi_at - ( 0.47 + 2.4 / (1 + exp((2.82 - 1.952 * (AT1R/AT1R_eq)) / 0.8)) );
% hatC_anp
f(64) = hatC_anp - ( 7.4052 - 6.554 / (1 + exp(P_ra - 3.762)) ); 
% AGT
f(65) = AGT_p - ( k_AGT - PRA - log(2)/h_AGT * AGT );
% nu_AT1
f(66) = nu_AT1 - ( 10^(0.0102 - 0.95 * log10(AT1R / AT1R_eq)) );
% R_sec
f(67) = R_sec - ( N_rs * nu_mdsod * nu_rsna * nu_AT1 );
% PRC
f(68) = PRC_p - ( R_sec - log(2)/h_renin * PRC );
% PRA
f(69) = PRA - ( PRC * X_PRCPRA );
% AngI
f(70) = AngI_p - ( PRA - (c_ACE + c_Chym + c_NEP) * AngI - log(2)/h_AngI * AngI );
% AngII
f(71) = AngII_p - ( (c_ACE + c_Chym) * AngI - (c_ACE2 + c_IIIV + c_AT1R + c_AT2R) * AngII - log(2)/h_AngII * AngII );
% AT1R
f(72) = AT1R_p - ( c_AT1R * AngII - log(2)/h_AT1R * AT1R );
% AT2R
f(73) = AT2R_p - ( c_AT2R * AngII - log(2)/h_AT2R * AT2R );
% Ang17
f(74) = Ang17_p - ( c_NEP * AngI + c_ACE2 * AngII - log(2)/h_Ang17 * Ang17 );
% AngIV
f(75) = AngIV_p - ( c_IIIV * AngII - log(2)/h_AngIV * AngIV );
% R_aa
f(76) = R_aa - ( R_aass * beta_rsna * Sigma_tgf * Sigma_myo * Psi_AT1RAA * Psi_AT2RAA);
% R_ea
f(77) = R_ea - ( R_eass * Psi_AT1REA * Psi_AT2REA );
% ESTIMATE ----------------------------------------------------------------
% Sigma_myo
if     P_gh <= x1
    f(78) = Sigma_myo - ( y1 );
elseif P_gh > x1 && P_gh < x2
    f(78) = Sigma_myo - ( (y2 - y1)/(x2 - x1)*P_gh + (x2*y1 - x1*y2)/(x2 - x1) );
elseif P_gh >= x2
    f(78) = Sigma_myo - ( y2 );
end
% ESTIMATE ----------------------------------------------------------------
% Psi_AT1RAA
f(79) = Psi_AT1RAA - ( 0.8   + 0.2092 * (AT1R / AT1R_eq) - 0.00925 / (AT1R / AT1R_eq) );
% Psi_AT1REA
f(80) = Psi_AT1REA - ( 0.925 + 0.0835 * (AT1R / AT1R_eq) - 0.0085  / (AT1R / AT1R_eq) );
% Psi_AT2RAA
if     strcmp(gender,'male')
    f(81) = Psi_AT2RAA - ( 1 );
elseif strcmp(gender,'female')
%     f(81) = Psi_AT2RAA - ( 0.025 * (AT2R_eq - AT2R) + 1 );
    f(81) = Psi_AT2RAA - ( 0.9 + 0.1 * exp(-(AT2R/AT2R_eq - 1)) );
%     f(81) = Psi_AT2RAA - ( 1 );
end
% Psi_AT2REA
if     strcmp(gender,'male')
    f(82) = Psi_AT2REA - ( 1 );
elseif strcmp(gender,'female')
%     f(82) = Psi_AT2REA - ( 0.01  * (AT2R_eq - AT2R) + 1 );
    f(82) = Psi_AT2REA - ( 0.9 + 0.1 * exp(-(AT2R/AT2R_eq - 1)) );
%     f(82) = Psi_AT2REA - ( 1 );
end

end

% =========================================================================
% Constraint for fmincon so that sigma_myo = 1 at baseline.
% =========================================================================

function [cc,ceq] = base_value(pars_est)

% Sigma_myo = (y2 - y1)/(x2 - x1)*P_gh + (x2*y1 - x1*y2)/(x2 - x1)
% pars_est = [x1; x2; y1; y2];
x1 = pars_est(1); x2 = pars_est(2); y1 = pars_est(3); y2 = pars_est(4); 

% c  (x) <= 0
% ceq(x)  = 0
cc = [];
ceq = (y2 - y1)/(x2 - x1)*62 + (x2*y1 - x1*y2)/(x2 - x1) - 1;

end
































