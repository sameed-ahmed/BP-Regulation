% This simulates the model for varying blood pressure. The simulation will 
% run for a physiologically reasonable range of blood pressure, (60,190),
% and then plot all variables versus blood pressure.

% Certain variables do not reach steady state for perturbed blood pressure.
% As a result, this simulation is run by manually changing blood pressure,
% keeping that constant blood pressure for some time in order for most
% variables to reach steady state, simulating the model for this change in
% blood pressure followed by constant blood pressure, recording all 
% variables at that end time, and then plotting all variables versus 
% blood pressure.

% function [MAP, GFR] = run_sim_P_ma
function run_sim_Pma_manual_iter

close all

lower_time = 130;
upper_time = 170;
time_to_ss = 100;

gender = {'male', 'female'};
change = {'decrease', 'increase'};
% X = (variable, iteration, gender)
% T      = cell(2,200);
X      = zeros(82,lower_time+1+upper_time,2);
% X      = zeros(82,120+1,2);

for gg = 1:2 % gender
for cc = 1:2 % change

%% Parameters

if     strcmp(gender{gg}, 'male')
    gen = 1;
elseif strcmp(gender{gg}, 'female')
    gen = 0;
end
if     strcmp(change{cc}, 'decrease')
    ch = 1;
elseif strcmp(change{cc}, 'increase')
    ch = 0;
end

% Scaling factor
% Rat flow = Human flow x SF
if     strcmp(gender{gg}, 'male')
    SF = 4.5*10^(-3)*10^(3);
elseif strcmp(gender{gg}, 'female')
    SF = 2/3 * 4.5*10^(-3)*10^(3);
end

N_rsna    = 1;
R_aass    = 31.67 / SF;   % mmHg min / l
R_eass    = 51.66 / SF;   % mmHg min / l
P_B       = 18;           % mmHg
P_go      = 28;           % mmHg
C_gcf     = 0.00781 * SF;
if     strcmp(gender{gg}, 'male')
   eta_etapt = 0.8; 
elseif strcmp(gender{gg}, 'female')
   eta_etapt = 0.5; 
end
eta_epsdt = 0.5; 
if     strcmp(gender{gg}, 'male')
   eta_etacd = 0.93; 
elseif strcmp(gender{gg}, 'female')
   eta_etacd = 0.972; 
end
K_vd      = 0.00001;
K_bar     = 16.6 / SF;    % mmHg min / l
R_bv      = 3.4 / SF;     % mmHg min / l
T_adh     = 6;            % min
Phi_sodin = 0.126 * SF;   % mEq / min
C_K       = 5;            % mEq / l 
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
        c_IIIV; c_AT1R; c_AT2R; AT1R_eq; AT2R_eq; gen; ch; SF];

%% Solve DAE

% Initial value
% This initial condition is the steady state data value taken from
% experiments (CITE). Therefore, the initial condition of the derivative is
% 0.

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

% Load data for steady state initial value. 
% Need to first run transform_data.m on Jessica's data files.
% if     strcmp(gender{gg}, 'male')
%     load(  'male_ss_data.mat', 'SSdata');
% elseif strcmp(gender{gg}, 'female')
%     load('female_ss_data.mat', 'SSdata');
% end

if     strcmp(gender{gg}, 'male')
    load(  'male_ss_data_new_sigmamyo.mat', 'SSdata');
elseif strcmp(gender{gg}, 'female')
    load('female_ss_data_new_sigmamyo.mat', 'SSdata');
end

% Load estimated pars for sigma_myo.
load('male_pars_est_sigmamyo0.0.mat', 'pars_min');

names  = {'$rsna$'; '$\alpha_{map}$'; '$\alpha_{rap}$'; '$R_{r}$'; ...
          '$\beta_{rsna}$'; '$\Phi_{rb}$'; '$\Phi_{gfilt}$'; '$P_{f}$'; ...
          '$P_{gh}$'; '$\Sigma_{tgf}$'; '$\Phi_{filsod}$'; ...
          '$\Phi_{pt-sodreab}$'; '$\eta_{pt-sodreab}$'; ...
          '$\gamma_{filsod}$'; '$\gamma_{at}$'; '$\gamma_{rsna}$'; ...
          '$\Phi_{md-sod}$'; '$\Phi_{dt-sodreab}$'; ...
          '$\eta_{dt-sodreab}$'; '$\psi_{al}$'; '$\Phi_{dt-sod}$'; ...
          '$\Phi_{cd-sodreab}$'; '$\eta_{cd-sodreab}$'; ...
          '$\lambda_{dt}$'; '$\lambda_{anp}$'; '$\Phi_{u-sod}$'; ...
          '$\Phi_{win}$'; '$V_{ecf}$'; '$V_{b}$'; '$P_{mf}$'; ...
          '$\Phi_{vr}$'; '$\Phi_{co}$'; '$P_{ra}$'; '$vas$'; ...
          '$vas_{f}$'; '$vas_{d}$'; '$R_{a}$'; '$R_{ba}$'; '$R_{vr}$'; ...
          '$R_{tp}$'; '$P_{ma}$'; '$\epsilon_{aum}$'; '$a_{auto}$'; ...
          '$a_{chemo}$'; '$a_{baro}$'; '$C_{adh}$'; '$N_{adh}$'; ...
          '$N_{adhs}$'; '$\delta_{ra}$'; '$\Phi_{t-wreab}$'; ...
          '$\mu_{al}$'; '$\mu_{adh}$'; '$\Phi_{u}$'; '$M_{sod}$'; ...
          '$C_{sod}$'; '$\nu_{md-sod}$'; '$\nu_{rsna}$'; '$C_{al}$'; ...
          '$N_{al}$'; '$N_{als}$'; '$\xi_{k/sod}$'; '$\xi_{map}$'; ...
          '$\xi_{at}$'; '$\hat{C}_{anp}$'; '$AGT$'; '$\nu_{AT1}$'; ...
          '$R_{sec}$'; '$PRC$'; '$PRA$'; '$Ang I$'; '$Ang II$'; ...
          '$Ang II_{AT1R-bound}$'; '$Ang II_{AT2R-bound}$'; ...
          '$Ang (1-7)$'; '$Ang IV$'; '$R_{aa}$'; '$R_{ea}$'; ...
          '$\Sigma_{myo}$'; '$\Psi_{AT1R-AA}$'; '$\Psi_{AT1R-EA}$'; ...
          '$\Psi_{AT2R-AA}$'; '$\Psi_{AT2R-EA}$'};

% Initial condition for the variables and their derivatives. 
% System is initially at steady state, so the derivative is 0.
x0 = SSdata; x_p0 = zeros(82,1);
X(:,lower_time+1,gg) = x0;

% Time at which to keep steady state.
if     strcmp(change{cc}, 'decrease')
    tss = linspace(1,lower_time,lower_time);
elseif strcmp(change{cc}, 'increase')
    tss = linspace(1,upper_time,upper_time);
end
iter_range = length(tss);

% ODE options
options = odeset('MaxStep',50); % default is 0.1*abs(t0-tf)

for iter = 1:iter_range

    % Initial time (min); Final time (min); Time vector;
    t0 = 0; tf = tss(iter) + time_to_ss; tspan = [t0, tf];

    % Solve dae
    [t,x] = ode15i(@(t,x,x_p) ...
                   bp_reg_sim_Pma_manual(t,x,x_p, pars,pars_min,tss(iter),SSdata), ...
                   tspan, x0, x_p0, options);
    if     strcmp(change{cc}, 'decrease')
        X(:,lower_time+1-iter,gg) = x(end,:);
    elseif strcmp(change{cc}, 'increase')
        X(:,lower_time+1+iter,gg) = x(end,:);
    end

end % iteration

end % change
end % gender

%% Retrieve data and visualize

% Retrieve male and female.
% X = (iteration, variable, gender)
X_m = X(:,:,1);     X_f = X(:,:,2);
P_ma_m = X_m(41,:); P_ma_f = X_f(41,:);
% MAP = X_m(41,:);
% GFR = X_m(7 ,:);

% Set limits for x axis.
min1 = min(P_ma_m); min2 = min(P_ma_f); xlower = min(min1,min2);
max1 = max(P_ma_m); max2 = max(P_ma_f); xupper = max(max1,max2);

% y-axis limits
ylower = zeros(length(X_m(:,1)),1); yupper = ylower; 
for i = 1:length(ylower)
    ylower(i) = 0.9*min(min(X_m(i,:)), min(X_f(i,:)));
    yupper(i) = 1.1*max(max(X_m(i,:)), max(X_f(i,:)));
    if ylower(i) == yupper(i)
        ylower(i) = ylower(i) - 10^(-5); yupper(i) = yupper(i) + 10^(-5);
    end
end

f = gobjects(6,1);
s = gobjects(6,15);
% Loop through each set of subplots.
for i = 1:6
    f(i) = figure('pos',[750 500 650 450]);
    % This is to avoid the empty plots in the last subplot set.
    if i == 6
        last_plot = 7;
    else
        last_plot = 15;
    end
    % Loop through each subplot within a set of subplots.
    for j = 1:last_plot
        s(i,j) = subplot(3,5,j);
        
        plot(s(i,j), P_ma_m,X_m((i-1)*15 + j,:), P_ma_f,X_f((i-1)*15 + j,:));
        xlim([xlower, xupper])
        ylim([ylower((i-1)*15 + j), yupper((i-1)*15 + j)])
        
        xlabel(names(41), 'Interpreter','latex', 'FontSize',15)
        legend('Male', 'Female')
        title(names((i-1)*15 + j), 'Interpreter','latex', 'FontSize',15)
    end
end

% savefig(f, 'all_vars_vs_Pma.fig')
% savefig(f, 'all_vars_vs_Pma_new_sigmamyo.fig')
% savefig(f, 'all_vars_vs_Pma_new_sigmamyo_no_betarsna.fig')

% Plot GFR vs BP

g = figure('pos',[100 100 675 450]);
plot(P_ma_m,X_m(7,:))
xlim([xlower, xupper])

set(gca,'FontSize',14)
xlabel('BP (mmHg)'  ,'FontSize',14, 'FontWeight','bold')
ylabel('GFR (l/min)','FontSize',14, 'FontWeight','bold')

% str = {'Normal', 'conditions', '\downarrow'};
% annotation('textbox',[.15 .895 0 0],'String',str,...
%            'FontSize',18,'FitBoxToText','on');

% savefig(g, 'GFRvsBP_cal.fig')

end






























