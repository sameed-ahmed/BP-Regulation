% This simulates the blood pressure regulation model bp_reg_mod.m for
% perturbations in renal perfusion pressure.
% 
% Steady state data is calculated by solve_ss_scenario.m.

% Input
% scenario: 
% RPP_per:  change in renal perfusion pressure from baseline in mmHg 
% num_points: number of points for plotting resolution
% exact_per: index RRP in RPP_per to plot for all variables

% Output
% plots of all variables vs renal perfusion pressure and vs time

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
% Index of RPP to plot for all variables
exact_per = 3;

% Scenarios
% Normal  - Normal conditions
% m_RAS   - male RAS pars
% m_Reab  - male fractional sodium and water reabsorption
% Pri_Hyp - essential/primary hypertension
scenario = {'Normal'};
% scenario = {'Normal', 'm_RSNA', 'm_AT2R', 'm_RAS', 'm_Reab', ...
%             'm_RAS_m_Reab', 'm_RSNA_m_Reab'};
% scenario = {'Normal', 'm_RSNA', 'm_AT2R', 'm_RAS', 'm_Reab', ...
%             'm_RAS_m_Reab', 'm_RSNA_m_Reab', ...
%             'Pri_Hyp'};
num_scen = length(scenario);
% Index of scenario to plot for all variables
exact_scen = 1;

% Number of points for plotting resolution
num_points = 121;

% Species
spe_ind = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of variables
num_vars = 93;

% Initialize variables.
% X = (variables, points, sex, perturbation, scenario)
X = zeros(num_vars,num_points,2,num_per,num_scen);
% Retrieve male/female. 
% X_m/f = (variables, points, perturbation, scenario)
X_m = zeros(num_vars,num_points,num_per,num_scen);
X_f = zeros(num_vars,num_points,num_per,num_scen);

% Need to store male and female RPP for plotting later.
% RPP = (sex, scenario)
RPP = zeros(2,num_scen);

species = {'human', 'rat'   };
sex     = {'male' , 'female'};

for per_ind = 1:num_per  % perturbation
for sce_ind = 1:num_scen % scenario
for sex_ind = 1:2        % sex

%% Retrieve variables and parameters.

% Initial value
% This initial condition is the steady state data value taken from
% solve_ss_scenario.m.

% Set name for data file to be loaded based upon sex.    
load_data_name = sprintf('%s_%s_ss_data_scenario_%s.mat', ...
                         species{spe_ind},sex{sex_ind},scenario{exact_scen});
% Load data for steady state initial value. 
load(load_data_name, 'SSdata');

% Optional parameters
varargin_input = {'RPP',{RPP_per(per_ind), SSdata}, scenario{sce_ind},true, ...
                  'Denerve',{true, SSdata}, 'Fixed Water Intake',{true, SSdata}};

% Parameter input
pars = get_pars(species{spe_ind}, sex{sex_ind}, varargin_input{:});

% Renal Perfusion Pressure.
RPP(sex_ind,sce_ind) = SSdata(42);

%% Solve DAE

%% Variable names for plotting.
names  = {'$rsna$'; '$\alpha_{map}$'; '$\alpha_{rap}$'; '$R_{r}$'; ...
          '$\beta_{rsna}$'; '$\Phi_{rb}$'; '$\Phi_{gfilt}$'; '$P_{f}$'; ...
          '$P_{gh}$'; '$\Sigma_{tgf}$'; '$\Phi_{filsod}$'; ...
          '$\Phi_{pt-sodreab}$'; '$\eta_{pt-sodreab}$'; ...
          '$\gamma_{filsod}$'; '$\gamma_{at}$'; '$\gamma_{rsna}$'; ...
          '$\Phi_{md-sod}$'; '$\Phi_{dt-sodreab}$'; ...
          '$\eta_{dt-sodreab}$'; '$\psi_{al}$'; '$\Phi_{dt-sod}$'; ...
          '$\Phi_{cd-sodreab}$'; '$\eta_{cd-sodreab}$'; ...
          '$\lambda_{dt}$'; '$\lambda_{anp}$'; '$\lambda_{al}$'; ...
          '$\Phi_{u-sod}$'; '$\Phi_{sodin}$'; '$V_{ecf}$'; '$V_{b}$'; ...
          '$P_{mf}$'; '$\Phi_{vr}$'; '$\Phi_{co}$'; '$P_{ra}$'; ...
          '$vas$'; '$vas_{f}$'; '$vas_{d}$'; '$R_{a}$'; '$R_{ba}$'; ...
          '$R_{vr}$'; '$R_{tp}$'; '$P_{ma}$'; '$\epsilon_{aum}$'; ...
          '$a_{auto}$'; '$a_{chemo}$'; '$a_{baro}$'; '$C_{adh}$'; ...
          '$N_{adh}$'; '$N_{adhs}$'; '$\delta_{ra}$'; ...
          '$M_{sod}$'; '$C_{sod}$'; '$\nu_{md-sod}$'; '$\nu_{rsna}$'; ...
          '$C_{al}$'; '$N_{al}$'; '$N_{als}$'; '$\xi_{k/sod}$'; ...
          '$\xi_{map}$'; '$\xi_{at}$'; '$\hat{C}_{anp}$'; '$AGT$'; ...
          '$\nu_{AT1}$'; '$R_{sec}$'; '$PRC$'; '$PRA$'; '$Ang I$'; ...
          '$Ang II$'; '$Ang II_{AT1R-bound}$'; '$Ang II_{AT2R-bound}$'; ...
          '$Ang (1-7)$'; '$Ang IV$'; '$R_{aa}$'; '$R_{ea}$'; ...
          '$\Sigma_{myo}$'; '$\Psi_{AT1R-AA}$'; '$\Psi_{AT1R-EA}$'; ...
          '$\Psi_{AT2R-AA}$'; '$\Psi_{AT2R-EA}$'; ...
          '$\Phi_{pt-wreab}$'; '$\eta_{pt-wreab}$'; ...
          '$\mu_{pt-sodreab}$'; '$\Phi_{md-u}$'; '$\Phi_{dt-wreab}$'; ...
          '$\eta_{dt-wreab}$'; '$\mu_{dt-sodreab}$'; '$\Phi_{dt-u}$'; ...
          '$\Phi_{cd-wreab}$'; '$\eta_{cd-wreab}$'; ...
          '$\mu_{cd-sodreab}$'; '$\mu_{adh}$'; ...
          '$\Phi_{u}$'; '$\Phi_{win}$'};
%%

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
options = odeset('RelTol',1e-1, 'AbsTol',1e-2, 'MaxStep',1e-2);

% Solve dae
[t,x] = ode15i(@(t,x,x_p) ...
               bp_reg_mod(t,x,x_p,pars,tchange,varargin_input{:}), ...
               tspan, x0, x_p0, options);
t = t'; x = x';

% Store solution.
% X = (variables, points, sex, perturbation, scenario)
X(:,:,sex_ind,per_ind,sce_ind) = x;

end % sex
end % scenario
end % perturbation

%% Retrieve data and visualize

% X_m/f = (variables, points, perturbation, scenario)
X_m(:,:,:,:) = X(:,:,1,:,:);
X_f(:,:,:,:) = X(:,:,2,:,:);

% x-axis limits
xlower = t0; xupper = tend; 

% y-axis limits
ylower = zeros(num_vars); yupper = ylower; 
for i = 1:num_vars
    ylower(i) = 0.9*min(min(X_m(i,:,exact_per,exact_scen)), min(X_f(i,:,exact_per,exact_scen)));
    yupper(i) = 1.1*max(max(X_m(i,:,exact_per,exact_scen)), max(X_f(i,:,exact_per,exact_scen)));
    if abs(yupper(i)) < eps*100
        ylower(i) = -10^(-5); yupper(i) = 10^(-5);
    end
end

% Plot all variables vs time. ---------------------------------------------

f  = gobjects(7,1);
ss1 = gobjects(7,15);
% Loop through each set of subplots.
for i = 1:7
    f(i) = figure('pos',[750 500 650 450]);
    % This is to avoid the empty plots in the last subplot set.
    if i == 7
        last_plot = mod(num_vars, 15);
    else
        last_plot = 15;
    end
    % Loop through each subplot within a set of subplots.
    for j = 1:last_plot
        ss1(i,j) = subplot(3,5,j);
        ss1(i,j).Position = ss1(i,j).Position + [0 0 0.01 0];
        
        plot(ss1(i,j), t,X_m((i-1)*15+j,:,exact_per,exact_scen),'b' , ...
                      t,X_f((i-1)*15+j,:,exact_per,exact_scen),'r');
        
%         xlim([xlower, xupper])
        ylim([ylower((i-1)*15 + j), yupper((i-1)*15 + j)])

%         legend('Male', 'Female')
        title(names((i-1)*15 + j), 'Interpreter','latex', 'FontSize',15)
    end
end

% Plot renal perfusion pressure input vs time. ----------------------------

tplot   = [t0:1:tend];
RPPplot = zeros(1,length(tplot));
RPPplot(1        :tchange) = RPP(1,exact_scen);
RPPplot(tchange+1:tend+1 ) = RPP(1,exact_scen) + RPP_per(exact_per);
g = figure('pos',[100 100 675 450]);
plot(tplot,RPPplot, 'LineWidth',3)
xlabel('$t$ (min)', 'Interpreter','latex')
ylabel('$RPP$'    , 'Interpreter','latex')

% Plot data as in Hilliard 2011. ------------------------------------------

% Time average quantity from 10-30 minutes after perturbation in RPP.
% RPP at 80, 100, 120.
% Phi_rb = var(6), Phi_gfilt = var(7), Phi_u = var(92), Phi_usod = var(26)
% rel = value divided by baseline value; act = actual value

% X_m/f = (variables, points, perturbation, scenario)

time_int    = (tchange+10)*ppm+1:(tchange+30)*ppm+1;
time_points = length(time_int);
RBF_rel_m  = zeros(num_per,num_scen); RBF_rel_f  = zeros(num_per,num_scen);  
GFR_rel_m  = zeros(num_per,num_scen); GFR_rel_f  = zeros(num_per,num_scen); 
UF_rel_m   = zeros(num_per,num_scen); UF_rel_f   = zeros(num_per,num_scen); 
USOD_rel_m = zeros(num_per,num_scen); USOD_rel_f = zeros(num_per,num_scen); 

for sce_ind = 1:num_scen
    for per_ind = 1:num_per
        RBF_rel_m (per_ind,sce_ind) = (sum(X_m(6 , time_int, per_ind, sce_ind)) / time_points) ...
                                    / (sum(X_m(6 , time_int, 2 , sce_ind)) / time_points);
        GFR_rel_m (per_ind,sce_ind) = (sum(X_m(7 , time_int, per_ind, sce_ind)) / time_points) ...
                                    / (sum(X_m(7 , time_int, 2 , sce_ind)) / time_points);
        UF_rel_m  (per_ind,sce_ind) = (sum(X_m(92, time_int, per_ind, sce_ind)) / time_points) ...
                                    / (sum(X_m(92, time_int, 2 , sce_ind)) / time_points);
        USOD_rel_m(per_ind,sce_ind) = (sum(X_m(27, time_int, per_ind, sce_ind)) / time_points) ...
                                    / (sum(X_m(27, time_int, 2 , sce_ind)) / time_points);
        
        RBF_rel_f (per_ind,sce_ind) = (sum(X_f(6 , time_int, per_ind, sce_ind)) / time_points) ...
                                    / (sum(X_f(6 , time_int, 2 , sce_ind)) / time_points);
        GFR_rel_f (per_ind,sce_ind) = (sum(X_f(7 , time_int, per_ind, sce_ind)) / time_points) ...
                                    / (sum(X_f(7 , time_int, 2 , sce_ind)) / time_points);
        UF_rel_f  (per_ind,sce_ind) = (sum(X_f(92, time_int, per_ind, sce_ind)) / time_points) ...
                                    / (sum(X_f(92, time_int, 2 , sce_ind)) / time_points);
        USOD_rel_f(per_ind,sce_ind) = (sum(X_f(27, time_int, per_ind, sce_ind)) / time_points) ...
                                    / (sum(X_f(27, time_int, 2 , sce_ind)) / time_points);
    end
end

% RPP
RPP_m = RPP(1,exact_scen) + RPP_per; RPP_f = RPP(2,exact_scen) + RPP_per; 

% Data --------------------------------------------------------------------

% Yes AT2R
RBFdata_yes_at2r_rel_m  = [0.8894; 1.0000; 1.0609]; RBFdata_yes_at2r_rel_f  = [0.8562; 1.0000; 1.1045]; 
GFRdata_yes_at2r_rel_m  = [0.7816; 1.0000; 1.0083]; GFRdata_yes_at2r_rel_f  = [0.7395; 1.0000; 1.1357]; 
UFdata_yes_at2r_rel_m   = [0.6163; 1.0000; 1.4535]; UFdata_yes_at2r_rel_f   = [0.8097; 1.0000; 2.2327]; 
USODdata_yes_at2r_rel_m = [0.4000; 1.0000; 1.8744]; USODdata_yes_at2r_rel_f = [0.5979; 1.0000; 3.0815]; 
% Block AT2R
RBFdata_blk_at2r_rel_m  = [0.8652; 1.0000; 1.2381]; RBFdata_blk_at2r_rel_f  = [0.5404; 1.0000; 1.0961]; 
GFRdata_blk_at2r_rel_m  = [0.7476; 1.0000; 1.2907]; GFRdata_blk_at2r_rel_f  = [0.2642; 1.0000; 1.3678]; 
UFdata_blk_at2r_rel_m   = [0.7106; 1.0000; 1.7368]; UFdata_blk_at2r_rel_f   = [0.4597; 1.0000; 2.7119]; 
USODdata_blk_at2r_rel_m = [0.0000; 1.0000; 3.5712]; USODdata_blk_at2r_rel_f = [0.5161; 1.0000; 4.2295]; 

% Male combined array for ease of plotting
RBFdata_rel_m  = [zeros(3,1), RBFdata_yes_at2r_rel_m , RBFdata_blk_at2r_rel_m ];
GFRdata_rel_m  = [zeros(3,1), GFRdata_yes_at2r_rel_m , GFRdata_blk_at2r_rel_m ];
UFdata_rel_m   = [zeros(3,1), UFdata_yes_at2r_rel_m  , UFdata_blk_at2r_rel_m  ];
USODdata_rel_m = [zeros(3,1), USODdata_yes_at2r_rel_m, USODdata_blk_at2r_rel_m];
% Female combined array for ease of plotting
RBFdata_rel_f  = [zeros(3,1), RBFdata_yes_at2r_rel_f , RBFdata_blk_at2r_rel_f ];
GFRdata_rel_f  = [zeros(3,1), GFRdata_yes_at2r_rel_f , GFRdata_blk_at2r_rel_f ];
UFdata_rel_f   = [zeros(3,1), UFdata_yes_at2r_rel_f  , UFdata_blk_at2r_rel_f  ];
USODdata_rel_f = [zeros(3,1), USODdata_yes_at2r_rel_f, USODdata_blk_at2r_rel_f];

% Subplot -----------------------------------------------------------------

h(1) = figure('DefaultAxesFontSize',14);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.15, 5]);
s_rel1(1) = subplot(2,2,1); 
s_rel1(2) = subplot(2,2,2); 
s_rel1(3) = subplot(2,2,3);
s_rel1(4) = subplot(2,2,4); 

plot(s_rel1(1), RPP_m,RBF_rel_m     (:,exact_scen) ,'x-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
% xlim(s_rel1(1), [75,125]); set(s_rel1(1),'XTick', [80,100,120]);
% ylim(s_rel1(1), [0.6,1.2])
xlabel(s_rel1(1), 'RPP (mmHg)'); ylabel(s_rel1(1), 'RBF (relative)');
hold(s_rel1(1), 'on')
plot(s_rel1(1), RPP_m,RBFdata_rel_m (:,2         ) ,'o--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(s_rel1(1), RPP_f,RBF_rel_f     (:,exact_scen) ,'x-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(s_rel1(1), RPP_f,RBFdata_rel_f (:,2         ) ,'o--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8); 
hold(s_rel1(1), 'off')
[~, hobj, ~, ~] = legend(s_rel1(1), {'Male sim','Male data','Female sim','Female data'}, 'FontSize',7,'Location','Southeast');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
title(s_rel1(1), 'A')

plot(s_rel1(2), RPP_m,GFR_rel_m     (:,exact_scen) ,'x-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
% xlim(s_rel1(2), [75,125]); set(s_rel1(2),'XTick', [80,100,120]);
% ylim(s_rel1(2), [0.6,1.2])
xlabel(s_rel1(2), 'RPP (mmHg)'); ylabel(s_rel1(2), 'GFR (relative)');
hold(s_rel1(2), 'on')
plot(s_rel1(2), RPP_m,GFRdata_rel_m (:,2         ) ,'o--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(s_rel1(2), RPP_f,GFR_rel_f     (:,exact_scen) ,'x-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(s_rel1(2), RPP_f,GFRdata_rel_f (:,2         ) ,'o--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8); 
hold(s_rel1(2), 'off')
title(s_rel1(2), 'B')

plot(s_rel1(3), RPP_m,UF_rel_m      (:,exact_scen) ,'x-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
% xlim(s_rel1(3), [75,125]); set(s_rel1(3),'XTick', [80,100,120]);
% ylim(s_rel1(3), [0.0,3.5])
xlabel(s_rel1(3), 'RPP (mmHg)'); ylabel(s_rel1(3), 'UF (relative)');
hold(s_rel1(3), 'on')
plot(s_rel1(3), RPP_m,UFdata_rel_m  (:,2         ) ,'o--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(s_rel1(3), RPP_f,UF_rel_f      (:,exact_scen) ,'x-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(s_rel1(3), RPP_f,UFdata_rel_f  (:,2         ) ,'o--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8); 
hold(s_rel1(3), 'off')
title(s_rel1(3), 'C')

plot(s_rel1(4), RPP_m,USOD_rel_m    (:,exact_scen) ,'x-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
% xlim(s_rel1(4), [75,125]); set(s_rel1(4),'XTick', [80,100,120]);
% ylim(s_rel1(4), [0.0,3.5])
xlabel(s_rel1(4), 'RPP (mmHg)'); ylabel(s_rel1(4), 'USOD (relative)');
hold(s_rel1(4), 'on')
plot(s_rel1(4), RPP_m,USODdata_rel_m(:,2         ) ,'o--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(s_rel1(4), RPP_f,USOD_rel_f    (:,exact_scen) ,'x-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(s_rel1(4), RPP_f,USODdata_rel_f(:,2         ) ,'o--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8); 
hold(s_rel1(4), 'off')
title(s_rel1(4), 'D')

%% Save figures.

save_data_name = sprintf('quant_of_int_vs_RPP.fig' );
save_data_name = strcat('Figures/', save_data_name);
savefig(h, save_data_name)

end






























