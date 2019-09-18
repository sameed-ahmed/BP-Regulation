% This simulates the blood pressure regulation model bp_reg.m.
% It is adopted with modifications from Karaaslan 2005 and Leete 2018.
% 
% Steady state data is calculated by solve_ss_baseline.m or solve_ss_scenario.m.

% function [SSdata, f] = run_sim
function run_sim

close all

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

% Scenarios
% Normal         - Normal conditions
% m_RAS          - male RAS pars
% m_Reab         - male fractional sodium and water reabsorption
% m_RAS_&_m_Reab - male RAS pars & fractional sodium and water reabsorption
scenario = {'Normal', 'm_RAS', 'm_Reab', 'm_RAS_&_m_Reab'};
fixed_ss = 1;

num_vars = 92;

gender   = {'male', 'female'};

% Initialize variables and time.
SSDATA   = zeros(num_vars,2);
residual = zeros(num_vars,2);
X        = cell(1,2);
T        = cell(1,2);

% % Jacobian sparsity pattern
% [dfdy_s,dfdy_p_s] = jac_spar;

for gg = 1:2 % gender

%% Parameters

% Parameter input
pars = get_pars(gender{gg}, scenario{fixed_ss});

%% Drugs

% drugs = [Ang II inf rate fmol/(ml min), ACEi target level, ARB target level]

% drugs = [0, 1, 0]; % Total ACEi
% drugs = [0, 0, 1]; % Total ARB
% drugs = [0, 0.78, 0]; % Leete 2018 ACEi
% drugs = [0, 0, 0.67]; % Leete 2018 ARB

drugs = [0, 0, 0]; % No drug

%% Solve DAE

% Initial value
% This initial condition is the steady state data value taken from
% Karaaslan 2005 and Leete 2018. 

% Set name for data file to be loaded based upon gender and scenario.    
load_data_name = sprintf('%s_ss_data_scenario_%s.mat', gender{gg},scenario{fixed_ss});

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
% fact = 2;
% fact = (1-0.625);
% fact = (1-0.95);
% fact = 0.5;
fact = 1;
% fact_var = 'Phi_sodin';
fact_var = 'N_rsna';

% Time at which to keep steady state, change a parameter, etc.
tchange = 1440; days = 13;
% tchange = 10;

% Initial time (min); Final time (min);
t0 = 0*1440; tend = tchange + days*1440;
% t0 = 0; tend = tchange + 1000;

% Time vector
tspan = [t0, tend]; %linspace(t0,tf,N);

% ode options
options = odeset();
% options = odeset('Jacobian',@(t,x,x_p)jac_anal(pars, t,x,x_p));
% options = odeset('JPattern',{ dfdy_s{gg},dfdy_p_s{gg} });
% options = odeset('RelTol',1e-1, 'AbsTol',1e-4); % default is -3, -6
options = odeset('MaxStep',1); % default is 0.1*abs(t0-tf)
% options = odeset('RelTol',1e-2, 'AbsTol',1e-4, 'MaxStep',1e-0);

% Solve dae
[t,x] = ode15i(@(t,x,x_p) ...
               bp_reg_sim(t,x,x_p,pars,fixed_var_pars,SSdata,drugs,...
                          tchange,fact,fact_var,scenario{fixed_ss}), ...
               tspan, x0, x_p0, options);

T{gg} = t';
X{gg} = x';

% residual(:,g) = bp_reg(0,x0,x_p0,pars);

end % gender

%% Plot

% Retrieve male and female.
t_m = T{1}; t_f = T{2};
X_m = X{1}; X_f = X{2};
% t_f = t_m; X_f = X_m; 

% x-axis limits
xlower = t0; xupper = tend; 

% Convert minutes to days for longer simulations.
t_m = t_m/1440; t_f = t_f/1440; tchange = tchange/1440; 
xlower = xlower/1440; xupper = xupper/1440; 

% y-axis limits
ylower = zeros(length(X_m(:,1)),1); yupper = ylower; 
for i = 1:length(ylower)
    ylower(i) = 0.95*min( min(X_m(i,:)), min(X_f(i,:)) );
    yupper(i) = 1.05*max( max(X_m(i,:)), max(X_f(i,:)) );
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
        
        plot(s(i,j), t_m,X_m((i-1)*15 + j,:),'b', t_f,X_f((i-1)*15 + j,:),'r');
        
        xlim([xlower, xupper])
        ylim([ylower((i-1)*15 + j), yupper((i-1)*15 + j)])
        
% %         Minutes
%         xlabel('Time (min)')
%         Days
        ax = gca;
%         ax.XTick = (tchange+0*(1*1440) : 1440 : tchange+days*(1*1440));
        ax.XTick = (tchange+0*(1) : 1 : tchange+days*(1));
        ax.XTickLabel = {'0' ,'1' ,'2' ,'3' ,'4' ,'5' ,'6' ,...
                         '7' ,'8' ,'9' ,'10','11','12','13',...
                         '14','15','16','17','18','19','20',...
                         '21','22','23','24','25','26'};
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

% Save figures.

% if     fact == 1
%     save_data_name = sprintf('all_vars_baseline.fig');
% else
%     if     strcmp(fact_var,'Phi_sodin')
%         save_data_name = sprintf('all_vars_%gxPhisodin.fig',fact);
%     elseif strcmp(fact_var,'N_rsna'   )
%         save_data_name = sprintf('all_vars_%gxNrsna.fig'   ,fact);
%     end
% end
% save_data_name = strcat('Figures/', save_data_name);
% savefig(f, save_data_name)

end






























