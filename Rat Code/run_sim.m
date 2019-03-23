% This simulates the blood pressure regulation model bp_reg.m.
% 
% Parameters are given by:
% "Long-Term Mathematical Model Involving Renal Sympathetic Nerve Activity,
% Arterial Pressure, and Sodium Excretion" - 2005 - Karaaslan, et. al.
% "Sex-specific Long-term Blood Pressure Regulation: Modeling and Analysis"
% - 2018 - Leete, Layton.
% 
% Steady state data is calculated by solve_ss_numerical.m.

% function [SSdata, f] = run_sim
function run_sim

close all

num_vars = 82;

gender   = {'male', 'female'};
SSDATA   = zeros(num_vars,2);
residual = zeros(num_vars,2);
X        = cell(1,2);
T        = cell(1,2);

% % Jacobian sparsity pattern
% [dfdy_s,dfdy_p_s] = jac_spar;

for gg = 1:2 % gender

%% Parameters

if     strcmp(gender{gg}, 'male')
    gen = 1;
elseif strcmp(gender{gg}, 'female')
    gen = 0;
end

% Scaling factor
% Rat flow = Human flow x SF
if     strcmp(gender{gg}, 'male')
    SF = 4.5*10^(-3)*10^(3);
elseif strcmp(gender{gg}, 'female')
    SF = 2/3 * 4.5*10^(-3)*10^(3);
end

N_rsna      = 1;
R_aass      = 31.67 / SF;   % mmHg min / l
R_eass      = 51.66 / SF;   % mmHg min / l
P_B         = 18;           % mmHg
P_go        = 28;           % mmHg
C_gcf       = 0.00781 * SF;
if     strcmp(gender{gg}, 'male')
   eta_etapt   = 0.8; 
%    eta_etapt   = 0.5; % female
elseif strcmp(gender{gg}, 'female')
   eta_etapt   = 0.5; 
%    eta_etapt   = 0.8; % male
end
eta_epsdt   = 0.5; 
if     strcmp(gender{gg}, 'male')
   eta_etacd   = 0.93; 
%    eta_etacd   = 0.972; % female
elseif strcmp(gender{gg}, 'female')
   eta_etacd   = 0.972; 
%    eta_etacd   = 0.93; % male
end
K_vd        = 0.00001;
K_bar       = 16.6 / SF;    % mmHg min / l
R_bv        = 3.4 / SF;     % mmHg min / l
T_adh       = 6;            % min
Phi_sodin   = 0.126 * SF;   % mEq / min
C_K         = 5;            % mEq / l 
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
%     % male
%     X_PRCPRA = 135.59/17.312;
%     k_AGT    = 801.02;
%     c_ACE    = 0.096833;
%     c_Chym   = 0.010833;
%     c_NEP    = 0.012667;
%     c_ACE2   = 0.0026667;
%     c_IIIV   = 0.29800;
%     c_AT1R   = 0.19700;
%     c_AT2R   = 0.065667;
%     AT1R_eq  = 20.46;
%     AT2R_eq  = 6.82;
%     % male
end

pars = [N_rsna; R_aass; R_eass; P_B; P_go; C_gcf; eta_etapt; eta_epsdt; ...
        eta_etacd; K_vd; K_bar; R_bv; T_adh; Phi_sodin; C_K; T_al; ...
        N_rs; X_PRCPRA; h_renin; h_AGT; h_AngI; h_AngII; h_Ang17; ...
        h_AngIV; h_AT1R; h_AT2R; k_AGT; c_ACE; c_Chym; c_NEP; c_ACE2; ...
        c_IIIV; c_AT1R; c_AT2R; AT1R_eq; AT2R_eq; gen; SF];

%% Drugs

% drugs = [Ang II inf rate fmol/(ml min), ACEi target level, ARB target level]
% drugs = [Ang II infusion rate fmol/(ml min)]
% drugs = [13625, 0, 0]; % Rajagopalan 1996 male; 5 days
% drugs = [13625, 0, 0] * (10/7); % Mollnau 2002 male
% drugs = [13625, 0, 0] / 5; % Ran 2006 male - ?
% drugs = [13625, 0, 0] / 25; % Brown 1981 female; 7 days
% if     strcmp(gender{gg}, 'male')
%     drugs = [(3/3)*5492, 0, 0]; % Zimmerman 2015 male + female; 14 days
% elseif strcmp(gender{gg}, 'female')
%     drugs = [(2/3)*5492, 0, 0]; % Zimmerman 2015 male + female; 14 days
% end
% drugs = [0, 1]; % Total ACEi

% if     strcmp(gender{gg}, 'male')
%     drugs = [(3/3)*10984, 0, 0]; % Sampson 2008 male + female; 13 days
% elseif strcmp(gender{gg}, 'female')
%     drugs = [(2/3)*10984, 0, 0]; % Sampson 2008 male + female; 13 days
% end

% drugs = [0, 0.78, 0]; % Leete 2018 ACEi
% drugs = [0, 0, 0.67]; % Leete 2018 ARB

drugs = [0, 0, 0]; % No drug

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
% %     load('male_ss_data_female_sodreab.mat', 'SSdata'); % female
% elseif strcmp(gender{gg}, 'female')
%     load('female_ss_data.mat', 'SSdata');
% %     load('female_ss_data_male_sodreab.mat', 'SSdata'); % male
% %     load('female_ss_dtata_male_raas.mat', 'SSdata'); % male
% end

% if     strcmp(gender{gg}, 'male')
%     load(  'male_ss_data_new_sigmamyo.mat', 'SSdata');
% elseif strcmp(gender{gg}, 'female')
%     load('female_ss_data_new_sigmamyo.mat', 'SSdata');
% end

% if     strcmp(gender{gg}, 'male')
%     load(  'male_ss_data_new_Psi.mat', 'SSdata');
% elseif strcmp(gender{gg}, 'female')
%     load('female_ss_data_new_Psi.mat', 'SSdata');
% end

% if     strcmp(gender{gg}, 'male')
%     load(  'male_ss_data_scenario_AT2R-.mat', 'SSdata');
% elseif strcmp(gender{gg}, 'female')
%     load('female_ss_data_scenario_AT2R-.mat', 'SSdata');
% end

if     strcmp(gender{gg}, 'male')
    load(  'male_ss_data_scenario_Normal.mat', 'SSdata');
elseif strcmp(gender{gg}, 'female')
    load('female_ss_data_scenario_Normal.mat', 'SSdata');
end

% if     strcmp(gender{gg}, 'male')
%     load(  'male_ss_data_new_Phitwreab.mat', 'SSdata');
% elseif strcmp(gender{gg}, 'female')
%     load('female_ss_data_new_Phitwreab.mat', 'SSdata');
% end

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
x0 = SSdata; x_p0 = zeros(num_vars,1);

% Factor by which to change something.
fact = 3.5;

% Time at which to keep steady state, change a parameter, etc.
% tchange = 1440;
tchange = 10;
days = 13;

% Initial time (min); Final time (min);
% t0 = 0*1440; tend = tchange + days*1440;
t0 = 0; tend = tchange + 90;

% Time vector
tspan = [t0, tend]; %linspace(t0,tf,N);

% ode options
options = odeset();
% options = odeset('Jacobian',@(t,x,x_p)jac_anal(pars, t,x,x_p));
% options = odeset('JPattern',{ dfdy_s{gg},dfdy_p_s{gg} });
% options = odeset('RelTol',1e-1, 'AbsTol',1e-4); % default is -3, -6
options = odeset('MaxStep',0.1); % default is 0.1*abs(t0-tf)
% options = odeset('RelTol',1e-7, 'AbsTol',1e-10, 'MaxStep',1e-1);

% Solve dae
[t,x] = ode15i(@(t,x,x_p) ...
               bp_reg_sim(t,x,x_p,pars,drugs,tchange,fact), ...
               tspan, x0, x_p0, options);

T{gg} = t';
X{gg} = x';

% residual(:,g) = bp_reg(0,x0,x_p0,pars);

end % gender

%% Plot

% Retrieve male and female.
t_m = T{1}; t_f = T{2};
X_m = X{1}; X_f = X{2};

% x-axis limits
xlower = t0; xupper = tend; 

% % Convert minutes to days for longer simulations.
% t_m = t_m/1440; t_f = t_f/1440; tchange = tchange/1440; 
% xlower = xlower/1440; xupper = xupper/1440; 

% y-axis limits
ylower = zeros(length(X_m(:,1)),1); yupper = ylower; 
for i = 1:length(ylower)
    ylower(i) = 0.9*min(min(X_m(i,:)), min(X_f(i,:)));
    yupper(i) = 1.1*max(max(X_m(i,:)), max(X_f(i,:)));
    if abs(yupper(i)) < eps*100
        ylower(i) = -10^(-5); yupper(i) = 10^(-5);
    end
end

f = gobjects(6,1);
s = gobjects(6,15);
% Loop through each set of subplots.
for i = 1:6
    f(i) = figure; 
%     f(i) = figure('pos',[750 500 650 450]);
    % This is to avoid the empty plots in the last subplot set.
    if i == 6
        last_plot = 7;
    else
        last_plot = 15;
    end
    % Loop through each subplot within a set of subplots.
    for j = 1:last_plot
        s(i,j) = subplot(3,5,j);
        s(i,j).Position = s(i,j).Position + [0 0 0.01 0];
        
        plot(s(i,j), t_m,X_m((i-1)*15 + j,:),'b', t_f,X_f((i-1)*15 + j,:),'r');
        
        xlim([xlower, xupper])
        ylim([ylower((i-1)*15 + j), yupper((i-1)*15 + j)])
        
%         Minutes
        xlabel('Time (min)')
% %         Days
%         ax = gca;
% %         ax.XTick = (tchange+0*(1*1440) : 1440 : tchange+days*(1*1440));
%         ax.XTick = (tchange+0*(1) : 1 : tchange+days*(1));
%         ax.XTickLabel = {'0' ,'1' ,'2' ,'3' ,'4' ,'5' ,'6' ,...
%                          '7' ,'8' ,'9' ,'10','11','12','13',...
%                          '14','15','16','17','18','19','20',...
%                          '21','22','23','24','25','26'};
%         xlabel('Time (days)')
% %         Weeks
%         ax = gca;
%         ax.XTick = [tchange+0*(7*1440); tchange+1*(7*1440); ...
%                     tchange+2*(7*1440); tchange+3*(7*1440)];
%         ax.XTickLabel = {'0','1','2','3'};
%         xlabel('Time (weeks)')
        
%         legend('Male', 'Female')
        title(names((i-1)*15 + j), 'Interpreter','latex', 'FontSize',15)
    end
end

% Plot Mean Arterial Pressure vs Time

% % Data from Zimmerman 2015. MAP is in relative change.
% tdata     = [0    ,1    ,2    ,3    ,4    ,5    ,6    ,7    ,...
%              8    ,9    ,10   ,11   ,12   ,13   ,14    ];
% MAPdata_m = [1.000,1.018,1.061,1.114,1.123,1.167,1.167,1.237,...
%              1.316,1.333,1.333,1.368,1.404,1.430,1.430,];
% MAPdata_f = [1.000,1.010,1.069,1.176,1.225,1.275,1.304,1.304,...
%              1.333,1.392,1.402,1.422,1.441,1.451,1.441,];
% % Multiply MAP relative change by baseline.
% MAPdata_m = X_m(41,1) * MAPdata_m;
% MAPdata_f = X_f(41,1) * MAPdata_f;
% 
% g = figure('pos',[100 100 675 450]);
% plot(t_m,X_m(41,:),'b-', t_f,X_f(41,:),'r-', 'LineWidth',3)
% xlim([xlower, xupper])
% ylim([80, 160])
% set(gca,'FontSize',14)
% ax = gca;
% ax.XTick = (tchange+0*(1) : 1 : tchange+days*(1));
% ax.XTickLabel = {'0','1','2' ,'3' ,'4' ,'5' ,'6' ,'7', ...
%                  '8','9','10','11','12','13','14'};
% xlabel('t (days)', 'Interpreter','latex', 'FontSize',22, 'FontWeight','bold')
% ylabel(names(41) , 'Interpreter','latex', 'FontSize',22, 'FontWeight','bold')
% legend('Male sim','Female sim');
% hold on
% plot(tdata,MAPdata_m,'bx', 'DisplayName',  'Male data', 'MarkerSize',10, 'LineWidth',3)
% plot(tdata,MAPdata_f,'rx', 'DisplayName','Female data', 'MarkerSize',10, 'LineWidth',3)

% % Data from Sampson 2008. MAP is in difference from baseline.
% tdata     = [0+1  ,1+1  ,2+1  ,3+1  ,4+1  ,5+1  ,6+1  ,...
%              7+1  ,8+1  ,9+1  ,10+1 ,11+1 ,12+1 ,13+1 ];
% MAPdata_m = [0.035,7.218,18.33,19.48,17.76,14.59,19.58,...
%              26.18,28.87,29.54,31.26,34.71,36.53,42.18];
% MAPdata_f = [0.011,10.85,15.98,14.31,14.31,18.44,14.71,...
%              13.91,17.31,17.04,18.37,19.63,23.23,24.42];
% % Substract MAP by baseline.
% MAP_m = X_m(41,:) - X_m(41,1);
% MAP_f = X_f(41,:) - X_f(41,1);
% 
% g = figure('DefaultAxesFontSize',30, 'pos',[100 100 650 450]);
% plot(t_m,MAP_m,'-', 'Color',[0.203, 0.592, 0.835], 'LineWidth',5);
% 
% xlim([xlower, xupper])
% ylim([0, 45])
% ax = gca;
% ax.XTick = (tchange+0*(1) : 1 : tchange+days*(1));
% ax.XTickLabel = {'0','1','2' ,'3' ,'4' ,'5' ,'6' ,'7', ...
%                  '8','9','10','11','12','13','14'};
% xlabel('Time (days)')
% ylabel('Change in MAP (mmHg)')
% % 'FontSize',22, 'FontWeight','bold'
% hold on
% plot(t_f,MAP_f,'-', 'Color',[0.835, 0.203, 0.576], 'LineWidth',5)
% plot(tdata,MAPdata_m,'o', 'DisplayName',  'Male data', 'Color',[0.203, 0.592, 0.835], 'MarkerSize',12, 'LineWidth',5)
% plot(tdata,MAPdata_f,'o', 'DisplayName','Female data', 'Color',[0.835, 0.203, 0.576], 'MarkerSize',12, 'LineWidth',5)
% legend('Male sim','Female sim','Male data','Female data', 'Location','Northwest');
% hold off

% Save figures.

% savefig(f, 'all_vars.fig')

% savefig(f, 'all_vars_new_sigmamyo0.0.fig')

% savefig(f, 'all_vars_new_Phitwreab.fig')

% savefig(f, 'all_vars_stepwise_Phisodin.fig')
% savefig(f, 'all_vars_stepwise_Phisodin_female_sodreab.fig')

% savefig(f, 'all_vars_new_Phitwreab_new_Sigmamyo_AngII_inf.fig')
% savefig(g, 'Pma_vs_t_new_Phitwreab_new_Sigmamyo_AngII_inf.fig')

% savefig(f, 'all_vars_Phisodin_inc.fig')

% male_Pma   = X_m(41,end)
% female_Pma = X_f(41,end)
% 
% BP = round(X_m(41,end));
% save_fig_name = sprintf('BP=%s.fig', num2str(BP));
% savefig(f ,save_fig_name)

% savefig(f ,'0.06x_Phisodin_no_rsna.fig')


end





























