% This simulates the blood pressure regulation model bp_reg.m for Ang II infusion.
% 
% Steady state data is calculated by solve_ss_baseline.m or solve_ss_scenario.m.

% function [SSdata, f] = run_sim
function run_sim_AngII

close all

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Begin user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scenarios
% Normal - Normal conditions
% m_RSNA - male RSNA
% m_AT2R - male AT2R
% m_RAS  - male RAS pars
% m_Reab - male fractional sodium and water reabsorption
% m_RAS_&_m_Reab - male RAS pars & fractional sodium and water reabsorption
scenario = {'Normal', 'm_RSNA', 'm_AT2R', 'm_RAS', 'm_Reab', 'm_RSNA_&_m_Reab'};
num_scen = length(scenario);
fixed_ss = 1;

% Number of days to run simulation after change; Day at which to induce change;
days = 13; day_change = 1;
% Number of points for plotting resolution
N = ((days+1)*1440) / 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of variables
num_vars = 92;

% Initialize variables.
% X = (variables, points, gender, scenario)
X = zeros(num_vars,N,2,num_scen);
T = zeros(N,2,num_scen);

gender = {'male', 'female'};

for ss = 1:num_scen % scenario
for gg = 1:2        % gender

%% Parameters

if     strcmp(gender{gg}, 'male')
    gen = 1;
elseif strcmp(gender{gg}, 'female')
    gen = 0;
end

% Scaling factors
% Rat value = Human value x SF
% Note: This includes conversion of units.
if     strcmp(gender{gg}, 'male')
    SF_S = 9.69;  % sodium flow % karaaslan
    SF_R = 0.343; % resistance
    SF_V = 3;     % volume
elseif strcmp(gender{gg}, 'female')
    SF_S = 9.69;  % sodium flow % karaaslan
    SF_R = 0.537; % resistance
    SF_V = 2.4;   % volume
end

N_rsna      = 1;
if     strcmp(gender{gg}, 'male')
R_aass    = 10.87;   % mmHg min / ml
R_eass    = 17.74;   % mmHg min / ml
elseif strcmp(gender{gg}, 'female')
R_aass    = 17.02;   % mmHg min / ml
R_eass    = 27.76;   % mmHg min / ml
end
P_B         = 18;           % mmHg
P_go        = 28;           % mmHg
if     strcmp(gender{gg}, 'male')
    C_gcf     = 0.068;
elseif strcmp(gender{gg}, 'female')
    C_gcf     = 0.047;
end

% Male and female different parameters for fractional reabsorption
if     strcmp(gender{gg}, 'male')
    eta_ptsodreab_eq = 0.80; % karaaslan
    eta_dtsodreab_eq = 0.5; 
    eta_cdsodreab_eq = 0.93;
elseif strcmp(gender{gg}, 'female')
    if   strcmp(scenario{ss}, 'm_Reab'         ) || ...
         strcmp(scenario{ss}, 'm_RAS_&_m_Reab' ) || ...
         strcmp(scenario{ss}, 'm_RSNA_&_m_Reab')
    eta_ptsodreab_eq = 0.71; % male
    eta_dtsodreab_eq = 0.5; 
    eta_cdsodreab_eq = 0.93;
    else
    eta_ptsodreab_eq = 0.5; % calibrated
    eta_dtsodreab_eq = 0.5; 
    eta_cdsodreab_eq = 0.96;
    end
end
if     strcmp(gender{gg}, 'male')
    eta_ptwreab_eq = 0.86; 
    eta_dtwreab_eq = 0.60; 
    eta_cdwreab_eq = 0.78;
elseif strcmp(gender{gg}, 'female')
    if   strcmp(scenario{ss}, 'm_Reab'         ) || ...
         strcmp(scenario{ss}, 'm_RAS_&_m_Reab' ) || ...
         strcmp(scenario{ss}, 'm_RSNA_&_m_Reab')
    eta_ptwreab_eq = 0.80; % male 
    eta_dtwreab_eq = 0.60; 
    eta_cdwreab_eq = 0.78;
    else
    eta_ptwreab_eq = 0.5; % calibrated
    eta_dtwreab_eq = 0.6; 
    eta_cdwreab_eq = 0.91;
    end
end

K_vd      = 0.01;
K_bar     = 16.6 * SF_R;  % mmHg min / ml
R_bv      = 3.4 * SF_R;   % mmHg min / ml
T_adh     = 6;            % min
Phi_sodin = 1.2212;       % microEq / min
C_K       = 5;            % microEq / ml 
T_al      = 30;           % min LISTED AS 30 IN TABLE %listed as 60 in text will only change dN_al
N_rs      = 1;            % ng / ml / min

% RAS
h_renin     = 12;      % min
h_AGT       = 10*60;   % min
h_AngI      = 0.5;     % min
h_AngII     = 0.66;    % min
h_Ang17     = 30;      % min
h_AngIV     = 0.5;     % min
h_AT1R      = 12;      % min
h_AT2R      = 12;      % min

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
    AT1R_eq  = 20.4807902818665;
    AT2R_eq  = 6.82696474842298;
elseif strcmp(gender{gg}, 'female')
    if     strcmp(scenario{ss}, 'm_RAS') || strcmp(scenario{ss}, 'm_RAS_&_m_Reab')
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

% Parameter input.
pars = [N_rsna; R_aass; R_eass; P_B; P_go; C_gcf; eta_ptsodreab_eq; ...
        eta_dtsodreab_eq; eta_cdsodreab_eq; eta_ptwreab_eq; ...
        eta_dtwreab_eq; eta_cdwreab_eq; K_vd; K_bar; R_bv; T_adh; ...
        Phi_sodin; C_K; T_al; N_rs; X_PRCPRA; h_renin; h_AGT; h_AngI; ...
        h_AngII; h_Ang17; h_AngIV; h_AT1R; h_AT2R; k_AGT; c_ACE; ...
        c_Chym; c_NEP; c_ACE2; c_IIIV; c_AT1R; c_AT2R; AT1R_eq; ...
        AT2R_eq; gen; SF_S; SF_R; SF_V];

%% Drugs

% drugs = [Ang II inf rate fmol/(ml min), ACEi target level, ARB target level]

if     strcmp(gender{gg}, 'male')
    drugs = [2022, 0, 0]; % Sampson 2008 male + female; 13 days
%     drugs = [(3/3)*5492, 0, 0]; % Zimmerman 2015 male + female; 14 days
elseif strcmp(gender{gg}, 'female')
    drugs = [2060, 0, 0]; % Sampson 2008 male + female; 13 days
%     drugs = [(2/3)*5492, 0, 0]; % Zimmerman 2015 male + female; 14 days
end

%% Solve DAE

% Initial value
% This initial condition is the steady state data value taken from
% Karaaslan 2005 and Leete 2018. 

% Set name for data file to be loaded based upon gender and scenario.    
load_data_name = sprintf('%s_ss_data_scenario_%s.mat', gender{gg},scenario{ss});

% Load data for steady state initial value. 
load(load_data_name, 'SSdata');

% Retrieve and replace parameters in fixed variable equations.
% These are the shift parameters which ensure that effect variables are 1.
fixed_ind = [2, 10, 14, 24, 44, 49, 66, 71, 88];
fixed_var_pars = SSdata(fixed_ind);
SSdata(fixed_ind) = 1;

% Variable names for plotting.
names  = {'$rsna$'; '$\alpha_{map}$'; '$\alpha_{rap}$'; '$R_{r}$'; ...
          '$\beta_{rsna}$'; '$\Phi_{rb}$'; '$\Phi_{gfilt}$'; '$P_{f}$'; ...
          '$P_{gh}$'; '$\Sigma_{tgf}$'; '$\Phi_{filsod}$'; ...
          '$\Phi_{pt-sodreab}$'; '$\eta_{pt-sodreab}$'; ...
          '$\gamma_{filsod}$'; '$\gamma_{at}$'; '$\gamma_{rsna}$'; ...
          '$\Phi_{md-sod}$'; '$\Phi_{dt-sodreab}$'; ...
          '$\eta_{dt-sodreab}$'; '$\psi_{al}$'; '$\Phi_{dt-sod}$'; ...
          '$\Phi_{cd-sodreab}$'; '$\eta_{cd-sodreab}$'; ...
          '$\lambda_{dt}$'; '$\lambda_{anp}$'; '$\lambda_{al}$'; ...
          '$\Phi_{u-sod}$'; '$\Phi_{win}$'; '$V_{ecf}$'; '$V_{b}$'; ...
          '$P_{mf}$'; '$\Phi_{vr}$'; '$\Phi_{co}$'; '$P_{ra}$'; ...
          '$vas$'; '$vas_{f}$'; '$vas_{d}$'; '$R_{a}$'; '$R_{ba}$'; ...
          '$R_{vr}$'; '$R_{tp}$'; '$P_{ma}$'; '$\epsilon_{aum}$'; ...
          '$a_{auto}$'; '$a_{chemo}$'; '$a_{baro}$'; '$C_{adh}$'; ...
          '$N_{adh}$'; '$N_{adhs}$'; '$\delta_{ra}$'; ...
          '$\Phi_{pt-wreab}$'; '$\eta_{pt-wreab}$'; ...
          '$\mu_{pt-sodreab}$'; '$\Phi_{md-u}$'; '$\Phi_{dt-wreab}$'; ...
          '$\eta_{dt-wreab}$'; '$\mu_{dt-sodreab}$'; '$\Phi_{dt-u}$'; ...
          '$\Phi_{cd-wreab}$'; '$\eta_{cd-wreab}$'; ...
          '$\mu_{cd-sodreab}$'; '$\mu_{adh}$'; '$\Phi_{u}$'; ...
          '$M_{sod}$'; '$C_{sod}$'; '$\nu_{md-sod}$'; '$\nu_{rsna}$'; ...
          '$C_{al}$'; '$N_{al}$'; '$N_{als}$'; '$\xi_{k/sod}$'; ...
          '$\xi_{map}$'; '$\xi_{at}$'; '$\hat{C}_{anp}$'; '$AGT$'; ...
          '$\nu_{AT1}$'; '$R_{sec}$'; '$PRC$'; '$PRA$'; '$Ang I$'; ...
          '$Ang II$'; '$Ang II_{AT1R-bound}$'; '$Ang II_{AT2R-bound}$'; ...
          '$Ang (1-7)$'; '$Ang IV$'; '$R_{aa}$'; '$R_{ea}$'; ...
          '$\Sigma_{myo}$'; '$\Psi_{AT1R-AA}$'; '$\Psi_{AT1R-EA}$'; ...
          '$\Psi_{AT2R-AA}$'; '$\Psi_{AT2R-EA}$'};

% Initial condition for the variables and their derivatives. 
% System is initially at steady state, so the derivative is 0.
x0 = SSdata; x_p0 = zeros(num_vars,1);

% Factor by which to change something.
fact = 1;
fact_var = '';

% Time at which to keep steady state, change a parameter, etc.
tchange = day_change*1440;
% tchange = 10;

% Initial time (min); Final time (min);
t0 = 0*1440; tend = tchange + days*1440;
% t0 = 0; tend = tchange + 1000;

% Time vector
tspan = linspace(t0,tend,N);

% ode options
% options = odeset('RelTol',1e-1, 'AbsTol',1e-4); % default is -3, -6
options = odeset('MaxStep',1); % default is 0.1*abs(t0-tf)
% options = odeset('RelTol',1e-2, 'AbsTol',1e-4, 'MaxStep',1e-0);

% Solve dae
[t,x] = ode15i(@(t,x,x_p) ...
               bp_reg_sim(t,x,x_p,pars,fixed_var_pars,SSdata,drugs,...
                          tchange,fact,fact_var,scenario{ss})     ,...
               tspan, x0, x_p0, options);

% X = (variables, points, gender, scenario)
T(  :,gg,ss) = t';
X(:,:,gg,ss) = x';

end % gender
end % scenario

%% Plot

% Retrieve male and female.
% X_m/f = (variables, points, scenario)
t_m = T(  :,1,1); t_f = T(  :,2,1);
X_m = reshape(X(:,:,1,:), [num_vars,N,num_scen]);
X_f = reshape(X(:,:,2,:), [num_vars,N,num_scen]);

% x-axis limits
xlower = t0; xupper = tend; 

% Convert minutes to days for longer simulations.
t_m = t_m/1440; t_f = t_f/1440; tchange = tchange/1440; 
xlower = xlower/1440; xupper = xupper/1440; 

% y-axis limits
ylower = zeros(num_vars,1); yupper = ylower; 
for i = 1:length(ylower)
    ylower(i) = 0.95*min( min(X_m(i,:,fixed_ss)), min(X_f(i,:,fixed_ss)) );
    yupper(i) = 1.05*max( max(X_m(i,:,fixed_ss)), max(X_f(i,:,fixed_ss)) );
    if abs(yupper(i)) < eps*100
        ylower(i) = -10^(-5); yupper(i) = 10^(-5);
    end
end

f = gobjects(7,1);
s = gobjects(7,15);
% Loop through each set of subplots.
for i = 1:7
%     f(i) = figure; 
    f(i) = figure('pos',[750 500 650 450]);
    % This is to avoid the empty plots in the last subplot set.
    if i == 7
        last_plot = 2;
    else
        last_plot = 15;
    end
    % Loop through each subplot within a set of subplots.
    for j = 1:last_plot
        s(i,j) = subplot(3,5,j);
        s(i,j).Position = s(i,j).Position + [0 0 0.01 0];
        
        plot(s(i,j), t_m,X_m((i-1)*15 + j,:,fixed_ss),'b', ...
                     t_f,X_f((i-1)*15 + j,:,fixed_ss),'r');
        
        xlim([xlower, xupper])
        ylim([ylower((i-1)*15 + j), yupper((i-1)*15 + j)])
        
        ax.XTick = (tchange+0*(1) : 2 : tchange+days*(1));
        ax.XTickLabel = {'0' ,'2' ,'4' ,'6' ,'8' ,'10','12',...
                         '14','16','18','20','22','24','26'};

%         legend('Male', 'Female')
        title(names((i-1)*15 + j), 'Interpreter','latex', 'FontSize',15)
    end
end

% Plot Mean Arterial Pressure vs Time. ------------------------------------

% Data from Sampson 2008. MAP is in difference from baseline.
tdata     = [0+1  ,1+1  ,2+1  ,3+1  ,4+1  ,5+1  ,6+1  ,...
             7+1  ,8+1  ,9+1  ,10+1 ,11+1 ,12+1 ,13+1 ];
MAPdata_m = [0.035,7.218,18.33,19.48,17.76,14.59,19.58,...
             26.18,28.87,29.54,31.26,34.71,36.53,42.18];
MAPdata_f = [0.011,10.85,15.98,14.31,14.31,18.44,14.71,...
             13.91,17.31,17.04,18.37,19.63,23.23,24.42];
% Substract MAP by baseline.
MAP_m = reshape(X_m(42,:,:) - X_m(42,1,:), [N,num_scen]);
MAP_f = reshape(X_f(42,:,:) - X_f(42,1,:), [N,num_scen]);

g = figure('DefaultAxesFontSize',14);%, 'pos',[100 100 650 450]);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 2.5]);
plot(t_m,MAP_m(:,fixed_ss),'-', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3);
xlim([xlower, xupper]); ylim([0, 60]);
ax = gca;
ax.XTick = (tchange+0*(1) : 2 : tchange+days*(1));
% ax.XTickLabel = {'0','1','2' ,'3' ,'4' ,'5' ,'6' ,'7', ...
%                  '8','9','10','11','12','13','14'};
ax.XTickLabel = {'0','2','4','6','8','10','12','14'};
xlabel('Time (days)')
ylabel('\DeltaMAP (mmHg)')
hold on
plot(t_f,MAP_f(:,fixed_ss),'-', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3)
plot(tdata,MAPdata_m,'o', 'DisplayName',  'Male data', 'Color',[0.203, 0.592, 0.835], 'MarkerSize',6, 'LineWidth',2)
plot(tdata,MAPdata_f,'o', 'DisplayName','Female data', 'Color',[0.835, 0.203, 0.576], 'MarkerSize',6, 'LineWidth',2)
[~, hobj, ~, ~] = legend({'Male sim','Female sim','Male data','Female data'}, 'FontSize',7,'Location','Northwest');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);

% Plot male - female bar graph for each scenario --------------------------

% X_m/f = (variable, points, scenario)
deltaMAP_m = reshape(X_m(42,end,1:end) - X_m(42,1,1:end), [1,num_scen]);
deltaMAP_f = reshape(X_f(42,end,1:end) - X_f(42,1,1:end), [1,num_scen]);
MAP_comp = deltaMAP_m(1) - [deltaMAP_m(1), deltaMAP_f];
scen_comp = categorical({'M - M'              , 'M - F', ...
                         'M - F M RSNA', 'M - F M AT2R', ...
                         'M - F M RAS' , 'M - F M Reab', ...
                         'M - FM RSNA\newline& Reab'});
scen_comp = reordercats(scen_comp,{'M - M'              , 'M - F', ...
                                   'M - F M RSNA', 'M - F M AT2R', ...
                                   'M - F M RAS' , 'M - F M Reab', ...
                                   'M - FM RSNA\newline& Reab'});

h = figure('DefaultAxesFontSize',10);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.15, 3.5]);
s1(1) = subplot(1,2,1); 
s1(2) = subplot(1,2,2); 
% s1(2).Position = s1(2).Position + [0.0, 0.0, 0.0, -0.01];

plot(s1(1), t_m,MAP_m(:,fixed_ss),'-', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3);
xlim(s1(1), [xlower, xupper]); ylim(s1(1), [0, 60]);
xticks(s1(1), (tchange+0*(1) : 2 : tchange+days*(1))); xticklabels(s1(1), {'0','2','4','6','8','10','12','14'});
xlabel(s1(1), 'Time (days)', 'FontSize',14*1.1); ylabel(s1(1), '\DeltaMAP (mmHg)', 'FontSize',14*1.1);
hold(s1(1), 'on')
plot(s1(1), t_f,MAP_f(:,fixed_ss),'-', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3)
plot(s1(1), tdata,MAPdata_m,'o', 'DisplayName',  'Male data', 'Color',[0.203, 0.592, 0.835], 'MarkerSize',6, 'LineWidth',2)
plot(s1(1), tdata,MAPdata_f,'o', 'DisplayName','Female data', 'Color',[0.835, 0.203, 0.576], 'MarkerSize',6, 'LineWidth',2)
[~, hobj, ~, ~] = legend(s1(1), {'Male sim','Female sim','Male data','Female data'}, 'FontSize',7,'Location','Northwest');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
title(s1(1), 'A', 'FontSize',14)

bar(s1(2), scen_comp,MAP_comp,'k');
% set(gca,'xticklabel',scen_comp_text);
% xtickget = get(gca,'xticklabel');  
% set(gca,'xticklabel',xtickget,'fontsize',6)
% xtickangle(s1(2),90)
% xlim(s1(2), [1-1,6+1])
ylim(s1(2), [0,15])
xlabel(s1(2), 'Scenario', 'FontSize',14); ylabel(s1(2), '\DeltaMAP (mmHg)', 'FontSize',14);
% hAxes.XAxis.FontSize = 6;
title(s1(2), 'B', 'FontSize',14)

% Save figures. -----------------------------------------------------------

% save_data_name = sprintf('all_vars_AngII_inf.fig');
% save_data_name = strcat('Figures/', save_data_name);
% savefig([f;g;h], save_data_name)

end






























