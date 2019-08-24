% This simulates the blood pressure regulation model bp_reg_RPP.m for
% perturbations in renal perfusion pressure over the entire physiological range.
% 
% Steady state data is calculated by solve_ss_baseline.m or solve_ss_scenario.m.

function run_sim_RPP_whole

close all

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Begin user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Renal perfusion pressure perturbation
% Enter postive for increase or negative for decrease.
lower_per = -40; upper_per = 100;
% Perturbation increment
inc_per   = 5;
RPP_per   = [lower_per : inc_per : upper_per]';
num_per   = length(RPP_per);
% Index for baseline perturbation
bl_per    = (0 - RPP_per(1)) / inc_per + 1;

% Scenarios
% Normal          - normal conditions
% Denerve         - cut off rsna from kidney
% Denerve & AT2R- - cut off rsna from kidney and block AT2R
% Denerve & No Myo      - cut off rsna from kidney and block myogenic response
% Denerve & No Myo      - cut off rsna from kidney and block tubuloglomerular feedback
% Denerve & No Myo, TGF - cut off rsna from kidney and block myogenic response and tubuloglomerular feedback
scenario = {'Normal', 'Denerve', 'Denerve & AT2R-', ...
            'Denerve & Linear Myo', 'Denerve & No Myo', 'Denerve & No TGF', 'Denerve & No Myo, TGF'};
num_scen = length(scenario);

% Number of points for plotting resolution
num_points = 121;

% Index of RPP to plot for all variables
exact_per = 3;

% Index of scenario to plot for all variables
% Scenario 'Denerve' is the one from Hilliard 2011.
exact_scen = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of variables
% 1 less due to fixed water intake.
num_vars = 92-1;

% Initialize variables.
% X = (variables, points, gender, perturbation, scenario)
X = zeros(num_vars+1,num_points,2,num_per,num_scen);
% Retrieve male/female. 
% X_m/f = (variables, points, perturbation, scenario)
X_m = zeros(num_vars+1,num_points,num_per,num_scen);
X_f = zeros(num_vars+1,num_points,num_per,num_scen);

% Need to store male and female RPP for plotting later.
% RPP = (gender, scenario)
RPP = zeros(2,num_scen);

gender = {'male', 'female'};

for pp = 1:num_per  % perturbation
for ss = 1:num_scen % scenario
for gg = 1:2        % gender

% Retrieve and replace parameters in fixed variable equations.
% Set name for data file to be loaded based upon gender.    
load_data_name = sprintf('%s_ss_data_scenario_Normal.mat', gender{gg});
load(load_data_name, 'SSdata');
fixed_ind = [2, 10, 14, 24, 44, 49, 66, 71, 88];
fixed_var_pars = SSdata(fixed_ind);
phicophico = SSdata(33); cadhcadh = SSdata(47);
fixed_var_pars = [fixed_var_pars; cadhcadh; phicophico];

%% Parameters

% Parameter input
pars = get_pars(gender{gg}, scenario{ss});

%% Drugs

% drugs = [Ang II inf rate fmol/(ml min), ACEi target level, ARB target level]
drugs = [0, 0, 0]; % No drug

%% Solve DAE

% Initial value
% This initial condition is the steady state data value taken from
% Karaaslan 2005 and Leete 2018. 

% Load data for steady state initial value. 
% Fixed variable parameters are only retrieved and replaced for scenarios
% which are solved for in the solve_ss_baseline because the parameter has
% changed for the fixed variable to remain 1.
% Otherwise scenarios which are solved for in solve_ss_scenario do not 
% require this because they load the fixed variable parameter from the
% baseline scenario, and they are perturbed scenarios in which the fixed
% variable is no longer 1.
if   strcmp(scenario{ss},'Denerve & AT2R-')
    load_data_name = sprintf('%s_ss_data_scenario_AT2R-.mat', gender{gg});
    load(load_data_name, 'SSdata');
else
    load_data_name = sprintf('%s_ss_data_scenario_Normal.mat', gender{gg});
    load(load_data_name, 'SSdata');
    SSdata(fixed_ind) = 1;
end

% Store water intake as an input and delete it as a variable.
Phi_win_input = SSdata(28);
SSdata(28) = '';

% Input Renal Perfusion Pressure.
RPP(gg,ss) = SSdata(42-1);

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

% Time at which to keep steady state, change a parameter, etc.
tchange = 10;

% Initial time (min); Final time (min); Points per minute;
t0 = 0; tend = tchange + 50; ppm = (num_points-1)/(tend-t0);

% Time vector
tspan = linspace(t0,tend,num_points);

% ode options
options = odeset();
% options = odeset('RelTol',1e-1, 'AbsTol',1e-4); % default is -3, -6
% options = odeset('MaxStep',1e-3); % default is 0.1*abs(tf-t0)
options = odeset('RelTol',1e-1, 'AbsTol',1e-2, 'MaxStep',1e-2);

% Solve dae
[t,x] = ode15i(@(t,x,x_p) ...
                bp_reg_sim_RPP(t,x,x_p,pars,fixed_var_pars,...
                               Phi_win_input,tchange,drugs,...
                               RPP(gg,ss),RPP_per(pp)     ,...
                               SSdata,scenario{ss})       ,...
                tspan, x0, x_p0, options);
t = t'; x = x';

% Add in Phi_win where it originally was.
Phi_win = Phi_win_input*ones(1,length(t));
x = [x(1:27,:); Phi_win; x(28:end,:)];
% Store solution.
% X = (variables, points, gender, perturbation, scenario)
X(:,:,gg,pp,ss) = x;

end % gender
end % scenario
end % perturbation

%% Retrieve data and visualize

% X_m/f = (variables, points, perturbation, scenario)
X_m(:,:,:,:) = X(:,:,1,:,:);
X_f(:,:,:,:) = X(:,:,2,:,:); % X_f = X_m;

% % x-axis limits
% xlower = t0; xupper = tend; 
% 
% % y-axis limits
% ylower = zeros(num_vars+1); yupper = ylower; 
% for i = 1:num_vars+1
%     ylower(i) = 0.9*min(min(X_m(i,:,exact_per,exact_scen)), min(X_f(i,:,exact_per,exact_scen)));
%     yupper(i) = 1.1*max(max(X_m(i,:,exact_per,exact_scen)), max(X_f(i,:,exact_per,exact_scen)));
%     if abs(yupper(i)) < eps*100
%         ylower(i) = -10^(-5); yupper(i) = 10^(-5);
%     end
% end
% 
% % Plot all variables vs time. ---------------------------------------------
% 
% f  = gobjects(7,1);
% s1 = gobjects(7,15);
% % Loop through each set of subplots.
% for i = 1:7
% %     f(i) = figure; 
%     f(i) = figure('pos',[750 500 650 450]);
%     % This is to avoid the empty plots in the last subplot set.
%     if i == 7
%         last_plot = 2;
%     else
%         last_plot = 15;
%     end
%     % Loop through each subplot within a set of subplots.
%     for j = 1:last_plot
%         s1(i,j) = subplot(3,5,j);
%         s1(i,j).Position = s1(i,j).Position + [0 0 0.01 0];
%         
%         plot(s1(i,j), t,X_m((i-1)*15+j,:,exact_per,exact_scen),'b' , ...
%                       t,X_f((i-1)*15+j,:,exact_per,exact_scen),'r');
%         
% %         xlim([xlower, xupper])
%         ylim([ylower((i-1)*15 + j), yupper((i-1)*15 + j)])
%         
% %         Minutes
% %         ax = gca;
% %         ax.XTick = (tchange : 10 : tend);
% %         ax.XTickLabel = {'0'  ,'20' ,'40' ,'60' ,'80' ,'100','120',...
% %                          '140','160','180','200','220','140','260',...
% %                          '280','300','320','340','360','380','400',...
% %                          '420','440','460','480','500','520'};
% %         xlabel('$t$ (min)', 'Interpreter','latex')
%         
% %         legend('Male', 'Female')
%         title(names((i-1)*15 + j), 'Interpreter','latex', 'FontSize',15)
%     end
% end
% 
% % Plot renal perfusion pressure input vs time. ----------------------------
% 
% tplot   = [t0:1:tend];
% RPPplot = zeros(1,length(tplot));
% RPPplot(1        :tchange) = RPP(1,2);
% RPPplot(tchange+1:tend+1 ) = RPP(1,2) + RPP_per(exact_per);
% g = figure('pos',[100 100 675 450]);
% plot(tplot,RPPplot, 'LineWidth',3)
% xlabel('$t$ (min)', 'Interpreter','latex')
% ylabel('$RPP$'    , 'Interpreter','latex')

% Plot data as in Hilliard 2011. ------------------------------------------

% Time average quantity from 10-30 minutes after perturbation in RPP.
% Phi_rb = var(6), Phi_gfilt = var(7), Phi_u = var(63), Phi_usod = var(27)
% R_aa = var(86), P_gh = var(9)

% X_m/f = (variables, points, perturbation, scenario)
time_int    = (tchange+10)*ppm+1:(tchange+30)*ppm+1;
time_points = length(time_int);
time_value = (tchange+150)*ppm+1;
RBF_m  = zeros(num_per,num_scen); RBF_f  = zeros(num_per,num_scen);  
GFR_m  = zeros(num_per,num_scen); GFR_f  = zeros(num_per,num_scen); 
UF_m   = zeros(num_per,num_scen); UF_f   = zeros(num_per,num_scen); 
USOD_m = zeros(num_per,num_scen); USOD_f = zeros(num_per,num_scen); 
RAA_m  = zeros(num_per,num_scen); RAA_f  = zeros(num_per,num_scen); 
PGH_m  = zeros(num_per,num_scen); PGH_f  = zeros(num_per,num_scen); 
for ss = 1:num_scen
    for pp = 1:num_per
        RBF_m (pp,ss) = (sum(X_m(6 , time_int, pp, ss)) / time_points) ...
                      / (sum(X_m(6 , time_int, bl_per , ss)) / time_points);
        GFR_m (pp,ss) = (sum(X_m(7 , time_int, pp, ss)) / time_points) ...
                      / (sum(X_m(7 , time_int, bl_per , ss)) / time_points);
        UF_m  (pp,ss) = (sum(X_m(63, time_int, pp, ss)) / time_points) ...
                      / (sum(X_m(63, time_int, bl_per , ss)) / time_points);
        USOD_m(pp,ss) = (sum(X_m(27, time_int, pp, ss)) / time_points) ...
                      / (sum(X_m(27, time_int, bl_per , ss)) / time_points);
        RAA_m (pp,ss) = (sum(X_m(86, time_int, pp, ss)) / time_points) ...
                      / (sum(X_m(86, time_int, bl_per , ss)) / time_points);
        PGH_m (pp,ss) = (sum(X_m(9 , time_int, pp, ss)) / time_points) ...
                      / (sum(X_m(9 , time_int, bl_per , ss)) / time_points);
        
        RBF_f (pp,ss) = (sum(X_f(6 , time_int, pp, ss)) / time_points) ...
                      / (sum(X_f(6 , time_int, bl_per , ss)) / time_points);
        GFR_f (pp,ss) = (sum(X_f(7 , time_int, pp, ss)) / time_points) ...
                      / (sum(X_f(7 , time_int, bl_per , ss)) / time_points);
        UF_f  (pp,ss) = (sum(X_f(63, time_int, pp, ss)) / time_points) ...
                      / (sum(X_f(63, time_int, bl_per , ss)) / time_points);
        USOD_f(pp,ss) = (sum(X_f(27, time_int, pp, ss)) / time_points) ...
                      / (sum(X_f(27, time_int, bl_per , ss)) / time_points);
        RAA_f (pp,ss) = (sum(X_m(86, time_int, pp, ss)) / time_points) ...
                      / (sum(X_m(86, time_int, bl_per , ss)) / time_points);
        PGH_f (pp,ss) = (sum(X_m(9 , time_int, pp, ss)) / time_points) ...
                      / (sum(X_m(9 , time_int, bl_per , ss)) / time_points);
    end
end

% RPP
RPP_m = RPP(1,2) + RPP_per; RPP_f = RPP(2,2) + RPP_per; 

% Autoregulatory range vertical lines
arr_lower = [93,93]; arr_upper = [173,173]; arr_line = [-1;14];

% Plots -------------------------------------------------------------------

g = gobjects(3,1);

g(1) = figure('DefaultAxesFontSize',14);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.15, 2.5]);
s1(1) = subplot(1,2,1); 
s1(2) = subplot(1,2,2); 

plot(s1(1), RPP_m,RBF_m (:,2),'-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3);
xlim(s1(1), [55 ,210]); xticks(s1(1), [60:30 :210]);
ylim(s1(1), [0  ,2.5]); yticks(s1(1), [0 :0.5:2.5]);
xlabel(s1(1), 'RPP (mmHg)'); ylabel(s1(1), 'RBF (relative)');
hold(s1(1), 'on')
plot(s1(1), RPP_f,RBF_f (:,2),'-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3);
legend(s1(1), 'Male','Female', 'Location','Northwest')
plot(s1(1), arr_lower,arr_line,'k--', 'LineWidth',1.5,'HandleVisibility','off'); 
plot(s1(1), arr_upper,arr_line,'k--', 'LineWidth',1.5,'HandleVisibility','off'); 
hold(s1(1), 'off')
title(s1(1), 'A')

plot(s1(2), RPP_m,GFR_m (:,2),'-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3);
xlim(s1(2), [55 ,210]); xticks(s1(2), [60:30 :210]);
ylim(s1(2), [0  ,2.5]); yticks(s1(2), [0 :0.5:2.5]);
xlabel(s1(2), 'RPP (mmHg)'); ylabel(s1(2), 'GFR (relative)');
hold(s1(2), 'on')
plot(s1(2), RPP_f,GFR_f (:,2),'-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3);
plot(s1(2), arr_lower,arr_line,'k--', 'LineWidth',1.5,'HandleVisibility','off'); 
plot(s1(2), arr_upper,arr_line,'k--', 'LineWidth',1.5,'HandleVisibility','off'); 
hold(s1(2), 'off')
title(s1(2), 'B')

% ---

g(2) = figure('DefaultAxesFontSize',14);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.15, 2.5]);
s2(1) = subplot(1,2,1); 
s2(2) = subplot(1,2,2); 

plot(s2(1), RPP_m,UF_m  (:,2),'-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3);
xlim(s2(1), [55 ,210]); xticks(s2(1), [60:30:210]);
ylim(s2(1), [0  ,13 ]); yticks(s2(1), [0 :3 :13 ]);
xlabel(s2(1), 'RPP (mmHg)'); ylabel(s2(1), 'UF (relative)');
hold(s2(1), 'on')
plot(s2(1), RPP_f,UF_f  (:,2),'-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3);
legend(s2(1), 'Male','Female', 'Location','Northwest')
plot(s2(1), arr_lower,arr_line,'k--', 'LineWidth',1.5,'HandleVisibility','off'); 
plot(s2(1), arr_upper,arr_line,'k--', 'LineWidth',1.5,'HandleVisibility','off'); 
hold(s2(1), 'off')
title(s2(1), 'A')

plot(s2(2), RPP_m,USOD_m(:,2),'-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3);
xlim(s2(2), [55 ,210]); xticks(s2(2), [60:30:210]);
ylim(s2(2), [0  ,13 ]); yticks(s2(2), [0 :3 :13 ]);
xlabel(s2(2), 'RPP (mmHg)'); ylabel(s2(2), 'USOD (relative)');
hold(s2(2), 'on')
plot(s2(2), RPP_f,USOD_f(:,2),'-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3);
plot(s2(2), arr_lower,arr_line,'k--', 'LineWidth',1.5,'HandleVisibility','off'); 
plot(s2(2), arr_upper,arr_line,'k--', 'LineWidth',1.5,'HandleVisibility','off'); 
hold(s2(2), 'off')
title(s2(2), 'B')

% ---

g(3) = figure('DefaultAxesFontSize',14);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.15, 2.5]);
s3(1) = subplot(1,2,1); 
s3(2) = subplot(1,2,2); 

plot(s3(1), RPP_m,RAA_m (:,2),'-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3);
xlim(s3(1), [55 ,210]); xticks(s3(1), [60:30 :210]);
ylim(s3(1), [0  ,2.5]); yticks(s3(1), [0 :0.5:2.5]);
xlabel(s3(1), 'RPP (mmHg)'); ylabel(s3(1), 'AAR (relative)');
hold(s3(1), 'on')
plot(s3(1), RPP_f,RAA_f (:,2),'-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3);
legend(s3(1), 'Male','Female', 'Location','Southeast')
plot(s3(1), arr_lower,arr_line,'k--', 'LineWidth',1.5,'HandleVisibility','off'); 
plot(s3(1), arr_upper,arr_line,'k--', 'LineWidth',1.5,'HandleVisibility','off'); 
hold(s3(1), 'off')
title(s3(1), 'A')

plot(s3(2), RPP_m,PGH_m (:,2),'-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3);
xlim(s3(2), [55 ,210]); xticks(s3(2), [60 :30 :210]);
ylim(s3(2), [0.7,1.4]); yticks(s3(2), [0.7:0.2:1.4]);
xlabel(s3(2), 'RPP (mmHg)'); ylabel(s3(2), 'GHP (relative)');
hold(s3(2), 'on')
plot(s3(2), RPP_f,PGH_f (:,2),'-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3);
plot(s3(2), arr_lower,arr_line,'k--', 'LineWidth',1.5,'HandleVisibility','off'); 
plot(s3(2), arr_upper,arr_line,'k--', 'LineWidth',1.5,'HandleVisibility','off'); 
hold(s3(2), 'off')
title(s3(2), 'B')

% ---

h = figure('DefaultAxesFontSize',14);%, 'pos',[100 100 675 450]);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.15, 5]);
ss1(1) = subplot(2,2,1); 
ss1(2) = subplot(2,2,2); 
ss1(3) = subplot(2,2,3);
ss1(4) = subplot(2,2,4); 

plot(ss1(1), RPP_m,RAA_m     (:,2) ,'-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim(ss1(1), [55,210]); xticks(ss1(1), [60 :30 :210]);
ylim(ss1(1), [0,2.5]); yticks(ss1(1), [0 :0.5:2.5]);
xlabel(ss1(1), 'RPP (mmHg)'); ylabel(ss1(1), 'AAR (relative)');
hold(ss1(1), 'on')
plot(ss1(1), RPP_f,RAA_f     (:,2) ,'-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(ss1(1), arr_lower,arr_line,'k--', 'LineWidth',1.5,'HandleVisibility','off'); 
plot(ss1(1), arr_upper,arr_line,'k--', 'LineWidth',1.5,'HandleVisibility','off'); 
hold(ss1(1), 'off')
legend(ss1(1), 'Male','Female', 'Location','Southeast');
title(ss1(1), 'A')

plot(ss1(2), RPP_m,PGH_m     (:,2) ,'-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim(ss1(2), [55,210]); xticks(ss1(2), [60 :30 :210]);
ylim(ss1(2), [0.7,1.4]); yticks(ss1(2), [0.7:0.2:1.4]);
xlabel(ss1(2), 'RPP (mmHg)'); ylabel(ss1(2), 'GHP (relative)');
hold(ss1(2), 'on')
plot(ss1(2), RPP_f,PGH_f     (:,2) ,'-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(ss1(2), arr_lower,arr_line,'k--', 'LineWidth',1.5,'HandleVisibility','off'); 
plot(ss1(2), arr_upper,arr_line,'k--', 'LineWidth',1.5,'HandleVisibility','off'); 
hold(ss1(2), 'off')
title(ss1(2), 'B')

plot(ss1(3), RPP_m,RBF_m      (:,2) ,'-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim(ss1(3), [55,210]); xticks(ss1(3), [60 :30 :210]);
ylim(ss1(3), [0.7,1.4]); yticks(ss1(3), [0.7:0.2:1.4]);
xlabel(ss1(3), 'RPP (mmHg)'); ylabel(ss1(3), 'RBF (relative)');
hold(ss1(3), 'on')
plot(ss1(3), RPP_f,RBF_f      (:,2) ,'-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(ss1(3), arr_lower,arr_line,'k--', 'LineWidth',1.5,'HandleVisibility','off'); 
plot(ss1(3), arr_upper,arr_line,'k--', 'LineWidth',1.5,'HandleVisibility','off'); 
hold(ss1(3), 'off')
title(ss1(3), 'C')

plot(ss1(4), RPP_m,GFR_m    (:,2) ,'-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim(ss1(4), [55,210]); xticks(ss1(4), [60 :30 :210]);
ylim(ss1(4), [0.0,2.5]); yticks(ss1(4), [0 :0.5:2.5]);
xlabel(ss1(4), 'RPP (mmHg)'); ylabel(ss1(4), 'GFR (relative)');
hold(ss1(4), 'on')
plot(ss1(4), RPP_f,GFR_f    (:,2) ,'-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(ss1(4), arr_lower,arr_line,'k--', 'LineWidth',1.5,'HandleVisibility','off'); 
plot(ss1(4), arr_upper,arr_line,'k--', 'LineWidth',1.5,'HandleVisibility','off'); 
hold(ss1(4), 'off')
title(ss1(4), 'D')

% ---

% Plot all scenarios
i = figure('DefaultAxesFontSize',14);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 2.5]);
plot(RPP_m,GFR_m(:,2) ,'-', 'LineWidth',3); % k-
xlim([55,210]); xticks([60:30:210]);
ylim([-1,5]); yticks([-1:2:5]);
xlabel('RPP (mmHg)'); ylabel('GFR (relative)');
hold on
plot(RPP_m,GFR_m(:,4) ,'-', 'LineWidth',3, 'MarkerIndices',1:4:length(RPP_m)); % lin myo k--o
plot(RPP_m,GFR_m(:,5) ,'-', 'LineWidth',3                                   ); % no  myo k--
plot(RPP_m,GFR_m(:,6) ,'-', 'LineWidth',3                                   ); % no  tgf k:
plot(RPP_m,GFR_m(:,7) ,'-', 'LineWidth',3, 'MarkerIndices',1:4:length(RPP_m)); % no  myo & tgf k--x
[~, hobj, ~, ~] = legend({'Full AR','Linear MR','No MR','No TGF','No MR and TGF'}, 'FontSize',7,'Location','Northwest');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
plot(arr_lower,arr_line,'k--', 'LineWidth',1.5,'HandleVisibility','off'); 
plot(arr_upper,arr_line,'k--', 'LineWidth',1.5,'HandleVisibility','off'); 
hold off

% % Save figures.
% 
% save_data_name = sprintf('quant_of_int_vs_RPP_whole.fig' );
% save_data_name = strcat('Figures/', save_data_name);
% savefig([g;h;i], save_data_name)

end






























