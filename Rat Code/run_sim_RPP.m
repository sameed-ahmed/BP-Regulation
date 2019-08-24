% This simulates the blood pressure regulation model bp_reg_RPP.m for
% perturbations in renal perfusion pressure.
% 
% Steady state data is calculated by solve_ss_baseline.m or solve_ss_scenario.m.

function run_sim_RPP

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
RPP_per = [-20; 0; 20];
num_per = length(RPP_per);

% Scenarios
% Normal          - normal conditions
% Denerve         - cut off rsna from kidney
% Denerve_&_AT2R- - cut off rsna from kidney and block AT2R
scenario = {'Normal', 'Denerve', 'Denerve & AT2R-'};
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

% x-axis limits
xlower = t0; xupper = tend; 

% y-axis limits
ylower = zeros(num_vars+1); yupper = ylower; 
for i = 1:num_vars+1
    ylower(i) = 0.9*min(min(X_m(i,:,exact_per,exact_scen)), min(X_f(i,:,exact_per,exact_scen)));
    yupper(i) = 1.1*max(max(X_m(i,:,exact_per,exact_scen)), max(X_f(i,:,exact_per,exact_scen)));
    if abs(yupper(i)) < eps*100
        ylower(i) = -10^(-5); yupper(i) = 10^(-5);
    end
end

% Plot all variables vs time. ---------------------------------------------

f  = gobjects(7,1);
s1 = gobjects(7,15);
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
        s1(i,j) = subplot(3,5,j);
        s1(i,j).Position = s1(i,j).Position + [0 0 0.01 0];
        
        plot(s1(i,j), t,X_m((i-1)*15+j,:,exact_per,exact_scen),'b' , ...
                      t,X_f((i-1)*15+j,:,exact_per,exact_scen),'r');
        
%         xlim([xlower, xupper])
        ylim([ylower((i-1)*15 + j), yupper((i-1)*15 + j)])
        
%         Minutes
%         ax = gca;
%         ax.XTick = (tchange : 10 : tend);
%         ax.XTickLabel = {'0'  ,'20' ,'40' ,'60' ,'80' ,'100','120',...
%                          '140','160','180','200','220','140','260',...
%                          '280','300','320','340','360','380','400',...
%                          '420','440','460','480','500','520'};
%         xlabel('$t$ (min)', 'Interpreter','latex')

%         legend('Male', 'Female')
        title(names((i-1)*15 + j), 'Interpreter','latex', 'FontSize',15)
    end
end

% Plot renal perfusion pressure input vs time. ----------------------------

tplot   = [t0:1:tend];
RPPplot = zeros(1,length(tplot));
RPPplot(1        :tchange) = RPP(1,2);
RPPplot(tchange+1:tend+1 ) = RPP(1,2) + RPP_per(exact_per);
g = figure('pos',[100 100 675 450]);
plot(tplot,RPPplot, 'LineWidth',3)
xlabel('$t$ (min)', 'Interpreter','latex')
ylabel('$RPP$'    , 'Interpreter','latex')

% Plot data as in Hilliard 2011. ------------------------------------------

% Time average quantity from 10-30 minutes after perturbation in RPP.
% RPP at 80, 100, 120.
% Phi_rb = var(6), Phi_gfilt = var(7), Phi_u = var(53), Phi_usod = var(26)

% X_m/f = (variables, points, perturbation, scenario)
time_int    = (tchange+10)*ppm+1:(tchange+30)*ppm+1;
time_points = length(time_int);
time_value  = (tchange+150)*ppm+1;
RBF_m  = zeros(num_per,num_scen); RBF_f  = zeros(num_per,num_scen);  
GFR_m  = zeros(num_per,num_scen); GFR_f  = zeros(num_per,num_scen); 
UF_m   = zeros(num_per,num_scen); UF_f   = zeros(num_per,num_scen); 
USOD_m = zeros(num_per,num_scen); USOD_f = zeros(num_per,num_scen); 
for ss = 1:num_scen
    for pp = 1:num_per
        RBF_m (pp,ss) = (sum(X_m(6 , time_int, pp, ss)) / time_points) ...
                      / (sum(X_m(6 , time_int, 2 , ss)) / time_points);
        GFR_m (pp,ss) = (sum(X_m(7 , time_int, pp, ss)) / time_points) ...
                      / (sum(X_m(7 , time_int, 2 , ss)) / time_points);
        UF_m  (pp,ss) = (sum(X_m(63, time_int, pp, ss)) / time_points) ...
                      / (sum(X_m(63, time_int, 2 , ss)) / time_points);
        USOD_m(pp,ss) = (sum(X_m(27, time_int, pp, ss)) / time_points) ...
                      / (sum(X_m(27, time_int, 2 , ss)) / time_points);
        
        RBF_f (pp,ss) = (sum(X_f(6 , time_int, pp, ss)) / time_points) ...
                      / (sum(X_f(6 , time_int, 2 , ss)) / time_points);
        GFR_f (pp,ss) = (sum(X_f(7 , time_int, pp, ss)) / time_points) ...
                      / (sum(X_f(7 , time_int, 2 , ss)) / time_points);
        UF_f  (pp,ss) = (sum(X_f(63, time_int, pp, ss)) / time_points) ...
                      / (sum(X_f(63, time_int, 2 , ss)) / time_points);
        USOD_f(pp,ss) = (sum(X_f(27, time_int, pp, ss)) / time_points) ...
                      / (sum(X_f(27, time_int, 2 , ss)) / time_points);
% % Psuedo steady state value at time_value mintues
%         RBF_m (pp,ss) = X_m(6 , time_value, pp, ss) ...
%                       / X_m(6 , time_value, 2 , ss);
%         GFR_m (pp,ss) = X_m(7 , time_value, pp, ss) ...
%                       / X_m(7 , time_value, 2 , ss);
%         UF_m  (pp,ss) = X_m(63, time_value, pp, ss) ...
%                       / X_m(63, time_value, 2 , ss);
%         USOD_m(pp,ss) = X_m(27, time_value, pp, ss) ...
%                       / X_m(27, time_value, 2 , ss);
%         
%         RBF_f (pp,ss) = X_f(6 , time_value, pp, ss) ...
%                       / X_f(6 , time_value, 2 , ss);
%         GFR_f (pp,ss) = X_f(7 , time_value, pp, ss) ...
%                       / X_f(7 , time_value, 2 , ss);
%         UF_f  (pp,ss) = X_f(63, time_value, pp, ss) ...
%                       / X_f(63, time_value, 2 , ss);
%         USOD_f(pp,ss) = X_f(27, time_value, pp, ss) ...
%                       / X_f(27, time_value, 2 , ss);
    end
end

% RPP
RPP_m = RPP(1,2) + RPP_per; RPP_f = RPP(2,2) + RPP_per; 

% Data --------------------------------------------------------------------

% Yes AT2R
RBFdata_yes_at2r_m  = [0.8894; 1.0000; 1.0609]; RBFdata_yes_at2r_f  = [0.8562; 1.0000; 1.1045]; 
GFRdata_yes_at2r_m  = [0.7816; 1.0000; 1.0083]; GFRdata_yes_at2r_f  = [0.7395; 1.0000; 1.1357]; 
UFdata_yes_at2r_m   = [0.6163; 1.0000; 1.4535]; UFdata_yes_at2r_f   = [0.8097; 1.0000; 2.2327]; 
USODdata_yes_at2r_m = [0.4000; 1.0000; 1.8744]; USODdata_yes_at2r_f = [0.5979; 1.0000; 3.0815]; 
% Block AT2R
RBFdata_blk_at2r_m  = [0.8652; 1.0000; 1.2381]; RBFdata_blk_at2r_f  = [0.5404; 1.0000; 1.0961]; 
GFRdata_blk_at2r_m  = [0.7476; 1.0000; 1.2907]; GFRdata_blk_at2r_f  = [0.2642; 1.0000; 1.3678]; 
UFdata_blk_at2r_m   = [0.7106; 1.0000; 1.7368]; UFdata_blk_at2r_f   = [0.4597; 1.0000; 2.7119]; 
USODdata_blk_at2r_m = [0.0000; 1.0000; 3.5712]; USODdata_blk_at2r_f = [0.5161; 1.0000; 4.2295]; 

% Male combined array for ease of plotting
RBFdata_m  = [zeros(3,1), RBFdata_yes_at2r_m , RBFdata_blk_at2r_m ];
GFRdata_m  = [zeros(3,1), GFRdata_yes_at2r_m , GFRdata_blk_at2r_m ];
UFdata_m   = [zeros(3,1), UFdata_yes_at2r_m  , UFdata_blk_at2r_m  ];
USODdata_m = [zeros(3,1), USODdata_yes_at2r_m, USODdata_blk_at2r_m];
% Female combined array for ease of plotting
RBFdata_f  = [zeros(3,1), RBFdata_yes_at2r_f , RBFdata_blk_at2r_f ];
GFRdata_f  = [zeros(3,1), GFRdata_yes_at2r_f , GFRdata_blk_at2r_f ];
UFdata_f   = [zeros(3,1), UFdata_yes_at2r_f  , UFdata_blk_at2r_f  ];
USODdata_f = [zeros(3,1), USODdata_yes_at2r_f, USODdata_blk_at2r_f];

% y-axis lower limits for uniformity accross scenarios
yRBF_lower  = min( min(min([RBF_m(:,:) ;RBF_f(:,:)]))  , min(min([RBFdata_m(:,2:num_scen) ;RBFdata_f(:,2:num_scen) ])) );
yGFR_lower  = min( min(min([GFR_m(:,:) ;GFR_f(:,:)]))  , min(min([GFRdata_m(:,2:num_scen) ;GFRdata_f(:,2:num_scen) ])) );
yUF_lower   = min( min(min([UF_m(:,:)  ;UF_f(:,:) ]))  , min(min([UFdata_m(:,2:num_scen)  ;UFdata_f(:,2:num_scen)  ])) );
yUSOD_lower = min( min(min([USOD_m(:,:);USOD_f(:,:)])) , min(min([USODdata_m(:,2:num_scen);USODdata_f(:,2:num_scen)])) );
% y-axis upper limits for uniformity accross scenarios
yRBF_upper  = max( max(max([RBF_m(:,:) ;RBF_f(:,:)]))  , max(max([RBFdata_m(:,2:num_scen) ;RBFdata_f(:,2:num_scen) ])) );
yGFR_upper  = max( max(max([GFR_m(:,:) ;GFR_f(:,:)]))  , max(max([GFRdata_m(:,2:num_scen) ;GFRdata_f(:,2:num_scen) ])) );
yUF_upper   = max( max(max([UF_m(:,:)  ;UF_f(:,:) ]))  , max(max([UFdata_m(:,2:num_scen)  ;UFdata_f(:,2:num_scen)  ])) );
yUSOD_upper = max( max(max([USOD_m(:,:);USOD_f(:,:)])) , max(max([USODdata_m(:,2:num_scen);USODdata_f(:,2:num_scen)])) );


% Subplot -----------------------------------------------------------------

g = figure('DefaultAxesFontSize',14);%, 'pos',[100 100 675 450]);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.15, 5]);
s1(1) = subplot(2,2,1); 
s1(2) = subplot(2,2,2); 
s1(3) = subplot(2,2,3);
s1(4) = subplot(2,2,4); 

plot(s1(1), RPP_m,RBF_m     (:,2) ,'x-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim(s1(1), [75,125]); set(s1(1),'XTick', [80,100,120]);
ylim(s1(1), [0.6,1.2])
xlabel(s1(1), 'RPP (mmHg)'); ylabel(s1(1), 'RBF (relative)');
hold(s1(1), 'on')
plot(s1(1), RPP_m,RBFdata_m (:,2) ,'o--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(s1(1), RPP_f,RBF_f     (:,2) ,'x-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(s1(1), RPP_f,RBFdata_f (:,2) ,'o--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8); 
hold(s1(1), 'off')
[~, hobj, ~, ~] = legend(s1(1), {'Male sim','Male data','Female sim','Female data'}, 'FontSize',7,'Location','Southeast');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
title(s1(1), 'A')

plot(s1(2), RPP_m,GFR_m     (:,2) ,'x-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim(s1(2), [75,125]); set(s1(2),'XTick', [80,100,120]);
ylim(s1(2), [0.6,1.2])
xlabel(s1(2), 'RPP (mmHg)'); ylabel(s1(2), 'GFR (relative)');
hold(s1(2), 'on')
plot(s1(2), RPP_m,GFRdata_m (:,2) ,'o--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(s1(2), RPP_f,GFR_f     (:,2) ,'x-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(s1(2), RPP_f,GFRdata_f (:,2) ,'o--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8); 
hold(s1(2), 'off')
title(s1(2), 'B')

plot(s1(3), RPP_m,UF_m      (:,2) ,'x-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim(s1(3), [75,125]); set(s1(3),'XTick', [80,100,120]);
ylim(s1(3), [0.0,3.5])
xlabel(s1(3), 'RPP (mmHg)'); ylabel(s1(3), 'UF (relative)');
hold(s1(3), 'on')
plot(s1(3), RPP_m,UFdata_m  (:,2) ,'o--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(s1(3), RPP_f,UF_f      (:,2) ,'x-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(s1(3), RPP_f,UFdata_f  (:,2) ,'o--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8); 
hold(s1(3), 'off')
title(s1(3), 'C')

plot(s1(4), RPP_m,USOD_m    (:,2) ,'x-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
xlim(s1(4), [75,125]); set(s1(4),'XTick', [80,100,120]);
ylim(s1(4), [0.0,3.5])
xlabel(s1(4), 'RPP (mmHg)'); ylabel(s1(4), 'USOD (relative)');
hold(s1(4), 'on')
plot(s1(4), RPP_m,USODdata_m(:,2) ,'o--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(s1(4), RPP_f,USOD_f    (:,2) ,'x-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(s1(4), RPP_f,USODdata_f(:,2) ,'o--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8); 
hold(s1(4), 'off')
title(s1(4), 'D')

% % Individual plot --------------------------------------------------------

% i(1) = figure('DefaultAxesFontSize',14);
% % i(1) = figure('DefaultAxesFontSize',20, 'pos',[100 450 650 450]);
% plot(RPP_m,RBF_m    (:,2) ,'x-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
% xlim([75,125]); xticks([80,100,120]);
% ylim([0.6,1.2])
% xlabel('RPP (mmHg)'); ylabel('RBF (relative)');
% title('A')
% hold on
% plot(RPP_m,RBFdata_m(:,2) ,'o--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
% plot(RPP_f,RBF_f    (:,2) ,'x-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
% plot(RPP_f,RBFdata_f(:,2) ,'o--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8); 
% legend('Male sim','Male data','Female sim','Female data', 'Location','Southeast')
% hold off
% 
% i(2) = figure('DefaultAxesFontSize',14);
% % i(2) = figure('DefaultAxesFontSize',20, 'pos',[100 450 650 450]);
% plot(RPP_m,GFR_m    (:,2) ,'x-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
% xlim([75,125]); xticks([80,100,120]);
% ylim([0.6,1.2])
% xlabel('RPP (mmHg)'); ylabel('GFR (relative)');
% title('B')
% hold on
% plot(RPP_m,GFRdata_m(:,2) ,'o--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
% plot(RPP_f,GFR_f    (:,2) ,'x-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
% plot(RPP_f,GFRdata_f(:,2) ,'o--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8); 
% legend('Male sim','Male data','Female sim','Female data', 'Location','Southeast')
% hold off
% 
% i(3) = figure('DefaultAxesFontSize',14);
% % i(3) = figure('DefaultAxesFontSize',20, 'pos',[100 450 650 450]);
% plot(RPP_m,UF_m    (:,2) ,'x-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
% xlim([75 ,125]); xticks([80,100,120]);
% ylim([0.0,3.5]); %yticks([0,1,2,3]); yticklabels({'0.0','1.0','2.0','3.0'});
% xlabel('RPP (mmHg)'); ylabel('UF (relative)');
% title('C')
% hold on
% plot(RPP_m,UFdata_m(:,2) ,'o--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
% plot(RPP_f,UF_f    (:,2) ,'x-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
% plot(RPP_f,UFdata_f(:,2) ,'o--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8); 
% legend('Male sim','Male data','Female sim','Female data', 'Location','Northwest')
% hold off
% 
% i(4) = figure('DefaultAxesFontSize',14);
% % i(4) = figure('DefaultAxesFontSize',20, 'pos',[100 450 650 450]);
% plot(RPP_m,USOD_m    (:,2) ,'x-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
% xlim([75,125]); xticks([80,100,120]);
% ylim([0.0,3.5]); %yticks([0,1,2,3]); yticklabels({'0.0','1.0','2.0','3.0'});
% xlabel('RPP (mmHg)'); ylabel('UNa^{+} (relative)');
% title('D')
% hold on
% plot(RPP_m,USODdata_m(:,2) ,'o--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
% plot(RPP_f,USOD_f    (:,2) ,'x-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
% plot(RPP_f,USODdata_f(:,2) ,'o--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8); 
% legend('Male sim','Male data','Female sim','Female data', 'Location','Northwest')
% hold off

% Save figures. -----------------------------------------------------------

% save_data_name = sprintf('quant_of_int_vs_RPP.fig' );
% save_data_name = strcat('Figures/', save_data_name);
% savefig(g, save_data_name)

end






























