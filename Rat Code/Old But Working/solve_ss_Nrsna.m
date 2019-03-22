% This script calculates the steady state solution to the system using
% fsolve. Some previous values are used as an initial guess. These are
% taken from Jessica, which are taken in part from some paper (Karaaslan
% 2005?).

function solve_ss_Nrsna

close all

% P_ma_range = (70:0.05:72);
% iteration = length(P_ma_range);
N_rsna_range = (-0.09:0.02:2.45);
iteration = length(N_rsna_range);

gender = {'male', 'female'};

% X = (variable, iteration, gender)
X = zeros(82,iteration,2);

for gg   = 1:1         % gender
for iter = 1:iteration % range

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

% N_rsna    = 1;
N_rsna    = N_rsna_range(iter);
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
C_K       = 5;               % mEq / l 
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
        c_IIIV; c_AT1R; c_AT2R; AT1R_eq; AT2R_eq; gen; SF];

%% Variables initial guess

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

% Load data for steady state initial value. 
% Need to first run transform_data.m on Jessica's data files.
% if     strcmp(gender{gg}, 'male')
%     SSdataIG = csvread(  'male_ss_data_IG.txt');
% elseif strcmp(gender{gg}, 'female')
%     SSdataIG = csvread('female_ss_data_IG.txt');
% end
if     strcmp(gender{gg}, 'male')
    load(  'male_ss_data_IG.mat', 'SSdataIG');
elseif strcmp(gender{gg}, 'female')
    load('female_ss_data_IG.mat', 'SSdataIG');
end

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

% Initial guess for the variables.
% Find the steady state solution, so the derivative is 0.
% Arbitrary value for time to input.
x0 = SSdataIG; x_p0 = zeros(82,1); t = 0;

% P_ma_ss = P_ma_range(iter);

%% Find steady state solution

% options = optimset(); %options = optimset('MaxFunEvals',8200+10000);
options = optimset('Display','off');
% [SSdata, residual, ...
%  exitflag, output] = fsolve(@(x) bp_reg_solve_pert_Pma(t,x,x_p0,pars,P_ma_ss), ...
%                             x0, options);
[SSdata, residual, ...
 exitflag, output] = fsolve(@(x) bp_reg_solve_Nrsna(t,x,x_p0,pars), ...
                            x0, options);

% Check for solver convergence.
if exitflag == 0
    disp('Solver did not converge.')
    disp(output)
end

% Check for imaginary solution.
if not (isreal(SSdata))
    disp('Imaginary number returned.')
end

% Set any values that are within machine precision of 0 equal to 0.
for i = 1:length(SSdata)
    if abs(SSdata(i)) < eps*100
        SSdata(i) = 0;
    end
end

X(:,iter,gg) = SSdata;
fprintf('Iteration = %s \n', num2str(N_rsna_range(iter)))

end % range
end % gender

%% Retrieve data and visualize

xplot = 'rsna';
% xplot = 'P_ma';

% Retrieve male and female.
% X = (variable, iteration, gender)
X_m    = X (:,:,1); X_f    = X (:,:,2);
if     strcmp(xplot, 'rsna')
    rsna_m = X_m(1 ,:); rsna_f = X_f(1 ,:);
elseif strcmp(xplot, 'P_ma')
    P_ma_m = X_m(41,:); P_ma_f = X_f(41,:);
end

% x-axis limits
if     strcmp(xplot, 'rsna')
    xlower = rsna_m(1); xupper = rsna_m(end); 
elseif strcmp(xplot, 'P_ma')
    xlower = P_ma_m(1); xupper = P_ma_m(end); 
end

% y-axis limits
X_f = X_m;
ylower = zeros(length(X_m(:,1)),1); yupper = ylower; 
for i = 1:length(ylower)
    ylower(i) = 0.9*min(min(X_m(i,:)), min(X_f(i,:)));
    yupper(i) = 1.1*max(max(X_m(i,:)), max(X_f(i,:)));
    if ylower(i) == yupper(i)
        ylower(i) = ylower(i) - 10^(-5); yupper(i) = yupper(i) + 10^(-5);
    end
end
X_f = zeros(size(X_f));

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
        
        if     strcmp(xplot, 'rsna')
            plot(s(i,j), rsna_m,X_m((i-1)*15 + j,:), rsna_f,X_f((i-1)*15 + j,:));
        elseif strcmp(xplot, 'P_ma')
            plot(s(i,j), P_ma_m,X_m((i-1)*15 + j,:), P_ma_f,X_f((i-1)*15 + j,:));
        end
        
        xlim([xlower, xupper])
        ylim([ylower((i-1)*15 + j), yupper((i-1)*15 + j)])
        
        if     strcmp(xplot, 'rsna')
            xlabel(names(1 ), 'Interpreter','latex', 'FontSize',15)
        elseif strcmp(xplot, 'P_ma')
            xlabel(names(41), 'Interpreter','latex', 'FontSize',15)
        end
%         legend('Male', 'Female')
        title(names((i-1)*15 + j), 'Interpreter','latex', 'FontSize',15)
    end
end

% savefig(f, 'all_vars_vs_Pma.fig')
% savefig(f, 'all_vars_vs_Pma_new_sigmamyo.fig')
% savefig(f, 'all_vars_vs_Pma_new_sigmamyo_no_betarsna.fig')

% % Plot GFR vs BP
% 
% g = figure('pos',[100 100 675 450]);
% plot(P_ma_m,X_m(7,:))
% % xlim([xlower, xupper])
% 
% set(gca,'FontSize',14)
% xlabel('BP (mmHg)'  ,'FontSize',14, 'FontWeight','bold')
% ylabel('GFR (l/min)','FontSize',14, 'FontWeight','bold')

% str = {'Normal', 'conditions', '\downarrow'};
% annotation('textbox',[.15 .895 0 0],'String',str,...
%            'FontSize',18,'FitBoxToText','on');

% savefig(g, 'GFRvsBP_cal.fig')

end






























