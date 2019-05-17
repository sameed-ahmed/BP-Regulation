% This simulates the blood pressure regulation model bp_reg.m.
% 
% Parameters are given by:
% "Long-Term Mathematical Model Involving Renal Sympathetic Nerve Activity,
% Arterial Pressure, and Sodium Excretion" - 2005 - Karaaslan, et. al.
% "Sex-specific Long-term Blood Pressure Regulation: Modeling and Analysis"
% - 2018 - Leete, Layton.
% 
% Steady state data is calculated by solve_ss_numerical.m.

function run_sim_RPP

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Begin user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Renal perfusion pressure perturbation
% Enter postive for increase, negative for decrease.
RPP_per = [-20; 0; 20];
num_per = length(RPP_per);

% Scenarios
% Normal          - normal conditions
% Denerve         - cut off rsna from kidney
% Denerve & AT2R- - cut off rsna from kidney and block AT2R
scenario = {'Normal', 'Denerve', 'Denerve & AT2R-'};
num_scen = length(scenario);

% Number of variables
num_vars   = 92-1;
% Number of points for plotting resolution
num_points = 121;

% Temporary single perfusion pressure at a time until I figure out a good 
% way to plot all three.
exact_per = 3;

% Temporary single scenario at a time until I figure out a good way to plot
% all three.
exact_scen = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gender = {'male', 'female'};

% Initialize variables.
% X = (variables, points, gender, perturbation, scenario)
X = zeros(num_vars+1,num_points,2,num_per,num_scen);
% Retrieve male/female. 
% X_m/f = (variables, points, perturbation, scenario)
X_m = zeros(num_vars+1,num_points,num_per,num_scen);
X_f = zeros(num_vars+1,num_points,num_per,num_scen);

% Need to store male and female RPP for plotting later.
RPP = zeros(2,1);

for pp = 1:num_per  % perturbation
for ss = 1:num_scen % scenario
for gg = 1:1        % gender

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

% Retrieve and replace parameters in fixed variable equations.
if     strcmp(gender{gg}, 'male')
    load(  'male_ss_data_scenario_Normal.mat', 'SSdata');
elseif strcmp(gender{gg}, 'female')
    load('female_ss_data_scenario_Normal.mat', 'SSdata');
end
fixed_ind = [2, 10, 14, 24, 44, 49, 66, 71, 88];
fixed_var_pars = SSdata(fixed_ind);
SF = 4.5*10^(-3)*10^(3);
phicophico = SSdata(33);
fixed_var_pars = [fixed_var_pars; phicophico];

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
% Rat resistance = Human resistance x SF
% Note: This includes conversion from l to ml.
if     strcmp(gender{gg}, 'male')
    SF_R = 0.343;
elseif strcmp(gender{gg}, 'female')
    SF_R = 0.522;
end

N_rsna    = 1;
% R_aass    = 31.67 / SF;   % mmHg min / ml
% R_eass    = 51.66 / SF;   % mmHg min / ml
R_aass    = 10.87;   % mmHg min / ml
R_eass    = 17.74;   % mmHg min / ml
P_B       = 18;           % mmHg
P_go      = 28;           % mmHg
% C_gcf     = 0.00781 * SF;
if     strcmp(gender{gg}, 'male')
    C_gcf     = 0.068;
elseif strcmp(gender{gg}, 'female')
    C_gcf     = 0.068;
end
if     strcmp(gender{gg}, 'male')
    eta_ptsodreab_eq = 0.93; 
    eta_dtsodreab_eq = 0.77; 
    eta_cdsodreab_eq = 0.15;
elseif strcmp(gender{gg}, 'female')
    eta_ptsodreab_eq = 0.5;
    eta_dtsodreab_eq = 0.5; 
    eta_cdsodreab_eq = 0.972;
end
if     strcmp(gender{gg}, 'male')
    eta_ptwreab_eq = 0.86; 
    eta_dtwreab_eq = 0.60; 
    eta_cdwreab_eq = 0.78;
elseif strcmp(gender{gg}, 'female')
    eta_ptwreab_eq = 0.5;
    eta_dtwreab_eq = 0.5; 
    eta_cdwreab_eq = 0.972;
end
K_vd      = 0.00001;
% K_bar     = 16.6 / SF;    % mmHg min / ml
K_bar     = 16.6 * SF_R;    % mmHg min / ml
% R_bv      = 3.4 / SF;     % mmHg min / ml
R_bv      = 3.4 * SF_R;     % mmHg min / ml
T_adh     = 6;            % min
% Phi_sodin = 1.2278;       % microEq / min
Phi_sodin = 2.3875;       % microEq / min
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

pars = [N_rsna; R_aass; R_eass; P_B; P_go; C_gcf; eta_ptsodreab_eq; ...
        eta_dtsodreab_eq; eta_cdsodreab_eq; eta_ptwreab_eq; ...
        eta_dtwreab_eq; eta_cdwreab_eq; K_vd; K_bar; R_bv; T_adh; ...
        Phi_sodin; C_K; T_al; N_rs; X_PRCPRA; h_renin; h_AGT; h_AngI; ...
        h_AngII; h_Ang17; h_AngIV; h_AT1R; h_AT2R; k_AGT; c_ACE; ...
        c_Chym; c_NEP; c_ACE2; c_IIIV; c_AT1R; c_AT2R; AT1R_eq; ...
        AT2R_eq; gen; SF; SF_R];

%% Drugs

% drugs = [Ang II inf rate fmol/(ml min), ACEi target level, ARB target level]
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

if   strcmp(scenario{ss},'Denerve & AT2R-')
    if     strcmp(gender{gg}, 'male')
        load(  'male_ss_data_scenario_AT2R-.mat', 'SSdata');
    elseif strcmp(gender{gg}, 'female')
        load('female_ss_data_scenario_AT2R-.mat', 'SSdata');
    end
else
    if     strcmp(gender{gg}, 'male')
        load(  'male_ss_data_scenario_Normal.mat', 'SSdata');
    elseif strcmp(gender{gg}, 'female')
        load('female_ss_data_scenario_Normal.mat', 'SSdata');
    end
    fixed_ind = [2, 10, 14, 24, 44, 49, 66, 71, 88];
    SSdata(fixed_ind) = 1;
end

% if     strcmp(gender{gg}, 'male')
%     load(  'male_ss_data_new_Phitwreab.mat', 'SSdata');
% elseif strcmp(gender{gg}, 'female')
%     load('female_ss_data_new_Phitwreab.mat', 'SSdata');
% end

% Store water intake as an input and delete it as a variable.
Phi_win_input = SSdata(28);
SSdata(28) = '';

% Input Renal Perfusion Pressure.
RPP(gg) = SSdata(42-1);

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
                bp_reg_sim_RPP(t,x,x_p,pars,fixed_var_pars,Phi_win_input,...
                               tchange,drugs,RPP(gg),RPP_per(pp),SSdata,scenario{ss}), ...
                tspan, x0, x_p0, options);
t = t'; x = x';

% Add in Phi_win where it originally was.
Phi_win = Phi_win_input*ones(1,length(t));
x = [x(1:27,:); Phi_win; x(28:end,:)];
% Store solution.
% X = (variables, points, gender, perturbation, scenario)
X(:,:,gg,pp,ss) = x;

% % Check this for mu_Na value
% fprintf('%s %s %d \n', scenario{ss},gender{gg},RPP_per(pp))
% if     strcmp(gender{gg},'male')
%     temp = [66.3766,  8.2376,  7.5987];
% %     temp = 3/2 * [27.3219, 13.6740, 13.2629]; % female
% elseif strcmp(gender{gg},'female')
%     temp = [27.3219, 13.6740, 13.2629];
% %     temp = 2/3*[66.3766,  8.2376,  7.5987]; % male
% end
% mu_Na = (X(12,(tchange+10)*ppm+1,gg,pp,ss) + X(18,10*ppm+1,gg,pp,ss) + X(22,10*ppm+1,gg,pp,ss)) / (temp(1) + temp(2) + temp(3))

end % gender
end % scenario
end % perturbation

%% Retrieve data and visualize

% X_m/f = (variables, points, perturbation, scenario)
X_m(:,:,:,:) = X(:,:,1,:,:);
X_f(:,:,:,:) = X(:,:,2,:,:); X_f = X_m;

% x-axis limits
xlower = t0; xupper = tend; 

% % Convert minutes to days for longer simulations.
% t = t/1440; tchange = tchange/1440; 
% xlower = xlower/1440; xupper = xupper/1440; 

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

% Plot renal perfusion pressure input vs time. ----------------------------
tplot   = [t0:1:tend];
RPPplot = zeros(1,length(tplot));
RPPplot(1        :tchange) = RPP(1);
RPPplot(tchange+1:tend+1 ) = RPP(1) + RPP_per(exact_per);
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
    end
end

% RPP
RPP_m = RPP(1) + RPP_per; RPP_f = RPP(2) + RPP_per; 

% Data

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

% h  = gobjects(num_scen,1);
% s2 = gobjects(num_scen,4);
% % Loop through each set of subplots.
% for ss = 1:num_scen
%     h(ss)    = figure('DefaultAxesFontSize',22, 'pos',[100 450 1100 700]); 
%     s2(ss,1) = subplot(2,2,1); 
%     s2(ss,2) = subplot(2,2,2); 
%     s2(ss,3) = subplot(2,2,3); 
%     s2(ss,4) = subplot(2,2,4); 
%     
%     % For the "Normal" secenario, there is no data.
%     if ss == 1
%     plot(s2(ss,1), RPP_m,RBF_m    (:,ss) ,'bx-' , ...
%                    RPP_f,RBF_f    (:,ss) ,'rx-' ,'LineWidth',3,'MarkerSize',10);
%     xlim(s2(ss,1), [75,125]); set(s2(ss,1),'XTick', [80,100,120]);
% %     ylim(s2(ss,1), [yRBF_lower;yRBF_upper])
%     ylim(s2(ss,1), [0.6,1.2])
%     xlabel(s2(ss,1), 'RPP (mmHg)'); ylabel(s2(ss,1), 'RBF (relative)');
%     legend(s2(ss,1), 'male sim','female sim', 'Location','Southeast')
%     title(s2(ss,1), 'A')
% 
%     plot(s2(ss,2), RPP_m,GFR_m    (:,ss) ,'bx-' , ...
%                    RPP_f,GFR_f    (:,ss) ,'rx-' ,'LineWidth',3,'MarkerSize',10);
%     xlim(s2(ss,2), [75,125]); set(s2(ss,2),'XTick', [80,100,120]);
% %     ylim(s2(ss,2), [yGFR_lower;yGFR_upper])
%     ylim(s2(ss,2), [0.6,1.2])
%     xlabel(s2(ss,2), 'RPP (mmHg)'); ylabel(s2(ss,2), 'GFR (relative)');
%     legend(s2(ss,2), 'male sim','female sim', 'Location','Southeast')
%     title(s2(ss,2), 'B')
% 
%     plot(s2(ss,3), RPP_m,UF_m    (:,ss) ,'bx-' , ...
%                    RPP_f,UF_f    (:,ss) ,'rx-' ,'LineWidth',3,'MarkerSize',10);
%     xlim(s2(ss,3), [75,125]); set(s2(ss,3),'XTick', [80,100,120]);
% %     ylim(s2(ss,3), [yUF_lower;yUF_upper])
%     ylim(s2(ss,3), [0.0,3.5])
%     xlabel(s2(ss,3), 'RPP (mmHg)'); ylabel(s2(ss,3), 'UF (relative)');
%     legend(s2(ss,3), 'male sim','female sim', 'Location','Northwest')
%     title(s2(ss,3), 'C')
% 
%     plot(s2(ss,4), RPP_m,USOD_m    (:,ss) ,'bx-' , ...
%                    RPP_f,USOD_f    (:,ss) ,'rx-' ,'LineWidth',3,'MarkerSize',10);
%     xlim(s2(ss,4), [75,125]); set(s2(ss,4),'XTick', [80,100,120]);
% %     ylim(s2(ss,4), [yUSOD_lower;yUSOD_upper])
%     ylim(s2(ss,4), [0.0,3.5])
%     xlabel(s2(ss,4), 'RPP (mmHg)'); ylabel(s2(ss,4), 'USOD (relative)');
%     legend(s2(ss,4), 'male sim','female sim', 'Location','Northwest')
%     title(s2(ss,4), 'D')
%     
% %     suptitle(scenario(ss))
%     else
%     plot(s2(ss,1), RPP_m,RBF_m    (:,ss) ,'bx-' , ...
%                    RPP_m,RBFdata_m(:,ss) ,'bo--', ...
%                    RPP_f,RBF_f    (:,ss) ,'rx-' , ...
%                    RPP_f,RBFdata_f(:,ss) ,'ro--','LineWidth',3,'MarkerSize',10);
%     xlim(s2(ss,1), [75,125]); set(s2(ss,1),'XTick', [80,100,120]);
% %     ylim(s2(ss,1), [yRBF_lower;yRBF_upper])
%     ylim(s2(ss,1), [0.6,1.2])
%     xlabel(s2(ss,1), 'RPP (mmHg)'); ylabel(s2(ss,1), 'RBF (relative)');
%     legend(s2(ss,1), 'male sim','male data','female sim','female data', 'Location','Southeast')
%     title(s2(ss,1), 'A')
% 
%     plot(s2(ss,2), RPP_m,GFR_m    (:,ss) ,'bx-' , ...
%                    RPP_m,GFRdata_m(:,ss) ,'bo--', ...
%                    RPP_f,GFR_f    (:,ss) ,'rx-' , ...
%                    RPP_f,GFRdata_f(:,ss) ,'ro--','LineWidth',3,'MarkerSize',10);
%     xlim(s2(ss,2), [75,125]); set(s2(ss,2),'XTick', [80,100,120]);
% %     ylim(s2(ss,2), [yGFR_lower;yGFR_upper])
%     ylim(s2(ss,2), [0.6,1.2])
%     xlabel(s2(ss,2), 'RPP (mmHg)'); ylabel(s2(ss,2), 'GFR (relative)');
%     legend(s2(ss,2), 'male sim','male data','female sim','female data', 'Location','Southeast')
%     title(s2(ss,2), 'B')
% 
%     plot(s2(ss,3), RPP_m,UF_m    (:,ss) ,'bx-' , ...
%                    RPP_m,UFdata_m(:,ss) ,'bo--', ...
%                    RPP_f,UF_f    (:,ss) ,'rx-' , ...
%                    RPP_f,UFdata_f(:,ss) ,'ro--','LineWidth',3,'MarkerSize',10);
%     xlim(s2(ss,3), [75,125]); set(s2(ss,3),'XTick', [80,100,120]);
% %     ylim(s2(ss,3), [yUF_lower;yUF_upper])
%     ylim(s2(ss,3), [0.0,3.5])
%     xlabel(s2(ss,3), 'RPP (mmHg)'); ylabel(s2(ss,3), 'UF (relative)');
%     legend(s2(ss,3), 'male sim','male data','female sim','female data', 'Location','Northwest')
%     title(s2(ss,3), 'C')
% 
%     plot(s2(ss,4), RPP_m,USOD_m    (:,ss) ,'bx-' , ...
%                    RPP_m,USODdata_m(:,ss) ,'bo--', ...
%                    RPP_f,USOD_f    (:,ss) ,'rx-' , ...
%                    RPP_f,USODdata_f(:,ss) ,'ro--','LineWidth',3,'MarkerSize',10);
%     xlim(s2(ss,4), [75,125]); set(s2(ss,4),'XTick', [80,100,120]);
% %     ylim(s2(ss,4), [yUSOD_lower;yUSOD_upper])
%     ylim(s2(ss,4), [0.0,3.5])
%     xlabel(s2(ss,4), 'RPP (mmHg)'); ylabel(s2(ss,4), 'USOD (relative)');
%     legend(s2(ss,4), 'male sim','male data','female sim','female data', 'Location','Northwest')
%     title(s2(ss,4), 'D')
% 
% %     suptitle(scenario(ss))
%     end
% end

i(1) = figure('DefaultAxesFontSize',30);
% i(1) = figure('DefaultAxesFontSize',30, 'pos',[100 450 650 450]);
plot(RPP_m,RBF_m    (:,2) ,'x-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',5, 'MarkerSize',12);
xlim([75,125]); %xticks([80,100,120]);
ylim([0.6,1.2])
xlabel('RPP (mmHg)'); ylabel('RBF (relative)');
title('A')
hold on
plot(RPP_m,RBFdata_m(:,2) ,'o--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',5, 'MarkerSize',12);
plot(RPP_f,RBF_f    (:,2) ,'x-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',5, 'MarkerSize',12);
plot(RPP_f,RBFdata_f(:,2) ,'o--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',5, 'MarkerSize',12); 
legend('Male sim','Male data','Female sim','Female data', 'Location','Southeast')
hold off

i(2) = figure('DefaultAxesFontSize',30);
% i(2) = figure('DefaultAxesFontSize',30, 'pos',[100 450 650 450]);
plot(RPP_m,GFR_m    (:,2) ,'x-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',5, 'MarkerSize',12);
xlim([75,125]); %xticks([80,100,120]);
ylim([0.6,1.2])
xlabel('RPP (mmHg)'); ylabel('GFR (relative)');
title('B')
hold on
plot(RPP_m,GFRdata_m(:,2) ,'o--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',5, 'MarkerSize',12);
plot(RPP_f,GFR_f    (:,2) ,'x-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',5, 'MarkerSize',12);
plot(RPP_f,GFRdata_f(:,2) ,'o--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',5, 'MarkerSize',12); 
legend('Male sim','Male data','Female sim','Female data', 'Location','Southeast')
hold off

i(3) = figure('DefaultAxesFontSize',30);
% i(3) = figure('DefaultAxesFontSize',30, 'pos',[100 450 650 450]);
plot(RPP_m,UF_m    (:,2) ,'x-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',5, 'MarkerSize',12);
xlim([75 ,125]); %xticks([80,100,120]);
ylim([0.0,3.5]); %yticks([0,1,2,3]); yticklabels({'0.0','1.0','2.0','3.0'});
xlabel('RPP (mmHg)'); ylabel('UF (relative)');
title('C')
hold on
plot(RPP_m,UFdata_m(:,2) ,'o--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',5, 'MarkerSize',12);
plot(RPP_f,UF_f    (:,2) ,'x-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',5, 'MarkerSize',12);
plot(RPP_f,UFdata_f(:,2) ,'o--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',5, 'MarkerSize',12); 
legend('Male sim','Male data','Female sim','Female data', 'Location','Northwest')
hold off

i(4) = figure('DefaultAxesFontSize',30);
% i(4) = figure('DefaultAxesFontSize',30, 'pos',[100 450 650 450]);
plot(RPP_m,USOD_m    (:,2) ,'x-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',5, 'MarkerSize',12);
xlim([75,125]); %xticks([80,100,120]);
ylim([0.0,3.5]); %yticks([0,1,2,3]); yticklabels({'0.0','1.0','2.0','3.0'});
xlabel('RPP (mmHg)'); ylabel('UNa^{+} (relative)');
title('D')
hold on
plot(RPP_m,USODdata_m(:,2) ,'o--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',5, 'MarkerSize',12);
plot(RPP_f,USOD_f    (:,2) ,'x-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',5, 'MarkerSize',12);
plot(RPP_f,USODdata_f(:,2) ,'o--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',5, 'MarkerSize',12); 
legend('Male sim','Male data','Female sim','Female data', 'Location','Northwest')
hold off

% % Save figures.
% 
% savefig(f, 'all_vars_RPP_inc.fig')
% savefig(i, 'quant_of_int_vs_RPP.fig')

end






























