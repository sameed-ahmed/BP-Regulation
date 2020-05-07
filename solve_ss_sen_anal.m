% This script performs local sensitivty analysis by perturbing parameters 
% and then calculating the steady state values of the blood pressure 
% regulation model bp_reg_solve_sen_anal.m for each inputted value. 
% 
% Steady state data for the intial guess is calculated by solve_ss_baseline.m.
% 
% Certain variables are then plotted as relative change versus each
% perturbed parameter.

function solve_ss_sen_anal

tic

close all

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Begin user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input of interest
% N_rsna       - par 3 ; -30%, +30%
% N_als        - par 18; -30%, +30%
% N_rs         - par 21; -30%, +30%
% R_aass       - par 4 ; -50%, +250%
% R_eass       - par 5 ; -50%, +250%
% c_ACE        - par 32; -50%, +50%
% c_AT1R       - par 37; -50%, +50%
% c_AT2R       - par 38; -50%, +50%
% % Psi_AT2RAA - par 42; -50%, +50%
% % Psi_AT2REA - par 43; -50%, +50%
% Parameter indices
in_ind      = [3;18;21;4;5;32;37;38];
% Range for percent change
range_lower = [30;30;30;50 ;50 ;50;50;50];
range_upper = [30;30;30;250;250;50;50;50];
par_iter    = length(in_ind);
% Parameter labels for bar graph
in_label = flip({'RSNA\newline-30% to +30%', ...
                 'Aldosterone Secretion\newline-30% to +30%', ...
                 'Renin Secretion\newline-30% to +30%', ...
                 'Afferent Arteriole Resistance\newline-50% to +250%', ...
                 'Efferent Arteriole Resistance\newline-50% to +250%', ...
                 'ACE Reaction Rate\newline-50% to +50%', ...
                 'AT1R Binding Rate\newline-50% to +50%', ...
                 'AT2R Binding Rate\newline-50% to +50%'});
% ycat = flip({'N_{rsna} -30% to +30%', 'N_{als} -30% to +30%', 'N_{rs} -30% to +30%', ...
%         'R_{aa} -50% to +250%', 'R_{ea} -50% to +250%', ...
%         'c_{ACE} -50% to +50%', 'c_{AT1R} -50% to +50%', 'c_{AT2R} -50% to +50%', ...
%         '\Psi_{AT2R-AA} -50% to +50%','\Psi_{AT2R-EA} -50% to +50%'});

% Output of interest
% Phi_gfilt         - var 7
% P_ma              - var 42
% Phi_u             - var 92
% R_aa              - var 73
% R_ea              - var 74
% eta_alpha-sodreab - var 13, 19, 23 (alpha = pt, dt, cd)
% eta_alpha-wreab   - var 81, 85, 89 (alpha = pt, dt, cd)
out_ind = [7;42;92;73;74;13;19;23;81;85;89];
out_num = length(out_ind);
% Variable labels for bar graph
out_label = {'GFR', 'MAP', 'UF', 'R_{AA}', 'R_{EA}', ...
             'F_{Na}^{pt}', 'F_{Na}^{dt}', 'F_{Na}^{cd}', ...
             'F_{W}^{pt}', 'F_{W}^{dt}', 'F_{W}^{cd}'};

% Scenarios
% Normal - Normal
% AngII  - Ang II infusion
% ACEi   - Angiotensin convernting enzyme inhibitor
% ARB    - Angiotensin receptor blocker
scenario = {'Normal', ...
            'AngII', 'ACEi', 'ARB', ...
            'Pri_Hyp'};
% Index of scenario to plot for all variables
fixed_ss = 5;

% Species
spe_ind = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of variables
num_vars = 93;

% Initialize variables.
% X = (variable, parameter, change, sex)
X = zeros(num_vars,par_iter,2,2);

species = {'human'   , 'rat'     };
sex     = {'male'    , 'female'  };
change  = {'decrease', 'increase'};

for iter    = 1:par_iter % par
for cha_ind = 1:2        % change
for sex_ind = 1:2        % sex

%% Parameters

varargin_input = {scenario{fixed_ss},true};

% Parameter input
pars = get_pars(species{spe_ind}, sex{sex_ind}, varargin_input{:});

% Change parameter of interest.
if     strcmp(change{cha_ind}, 'decrease')
    pars(in_ind(iter)) = pars(in_ind(iter)) * (1 - range_lower(iter)/100);
elseif strcmp(change{cha_ind}, 'increase')
    pars(in_ind(iter)) = pars(in_ind(iter)) * (1 + range_upper(iter)/100);
end

%% Drugs

% drugs = [Ang II inf rate fmol/(ml min), ACEi target level, ARB target level]
if     strcmp(scenario{fixed_ss}, 'AngII')
    if     strcmp(sex{sex_ind}, 'male'  )
        varargin_input = {'AngII',2022}; % Sampson 2008
    elseif strcmp(sex{sex_ind}, 'female')
        varargin_input = {'AngII',2060}; % Sampson 2008
    end
elseif strcmp(scenario{fixed_ss}, 'ACEi' )
        varargin_input = {'ACEi' ,0.78 }; % Leete 2018
elseif strcmp(scenario{fixed_ss}, 'ARB'  )
        varargin_input = {'ARB'  ,0.67 }; % Leete 2018
end

%% Variables initial guess

% Load data for steady state initial guess. 
% Set name for data file to be loaded based upon sex.    
load_data_name = sprintf('%s_%s_ss_data_scenario_%s.mat', ...
                         species{spe_ind},sex{sex_ind},scenario{fixed_ss});
load(load_data_name, 'SSdata');
SSdataIG = SSdata;
clear SSdata

% Order
% x  = [rsna; alpha_map; alpha_rap; R_r; beta_rsna; Phi_rb; Phi_gfilt; ...
%       P_f; P_gh; Sigma_tgf; Phi_filsod; Phi_ptsodreab; eta_ptsodreab; ...
%       gamma_filsod; gamma_at; gamma_rsna; Phi_mdsod; Phi_dtsodreab; ...
%       eta_dtsodreab; psi_al; Phi_dtsod; Phi_cdsodreab; eta_cdsodreab; ...
%       lambda_dt; lambda_anp; lambda_al; Phi_usod; Phi_sodin; V_ecf; ...
%       V_b; P_mf; Phi_vr; Phi_co; P_ra; vas; vas_f; vas_d; R_a; R_ba; ...
%       R_vr; R_tp; P_ma; epsilon_aum; a_auto; a_chemo; a_baro; C_adh; ...
%       N_adh; N_adhs; delta_ra; ...
%       M_sod; C_sod; nu_mdsod; nu_rsna; C_al; N_al; N_als; xi_ksod; ...
%       xi_map; xi_at; hatC_anp; AGT; nu_AT1; R_sec; PRC; PRA; AngI; ...
%       AngII; AT1R; AT2R; Ang17; AngIV; R_aa; R_ea; Sigma_myo; ...
%       Psi_AT1RAA; Psi_AT1REA; Psi_AT2RAA; Psi_AT2REA; ...
%       Phi_ptwreab; eta_ptwreab; mu_ptsodreab; ...
%       Phi_mdu; Phi_dtwreab; eta_dtwreab; mu_dtsodreab; Phi_dtu; ...
%       Phi_cdwreab; eta_cdwreab; mu_cdsodreab; mu_adh; Phi_u; Phi_win];

% Initial guess for the variables.
% Find the steady state solution, so the derivative is 0.
% Arbitrary value for time to input, greater than tchange + deltat.
x0 = SSdataIG; x_p0 = zeros(num_vars,1); t = 30;

% Time at which to change place holder.
tchange = 0;

%% Find steady state solution

options = optimset('Display','off');
[SSdata, residual, ...
 exitflag, output] = fsolve(@(x) ...
                            bp_reg_mod(t,x,x_p0,pars,tchange,varargin_input{:}), ...
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

% Relative percent change
SSdata = (SSdata - SSdataIG) ./ SSdataIG * 100;

% Store solution.
X(:,iter,cha_ind,sex_ind) = SSdata;

end % sex
end % change
end % par

%% Retrieve data and visualize

% Retrieve male/female.
% X_m/f = (variable, parameter, change)
X_m = X(:,:,:,1);
X_f = X(:,:,:,2);

% Plot individual. --------------------------------------------------------

f = gobjects(out_num,1);
for i = 1:out_num
    f(i) = figure('DefaultAxesFontSize',20, 'pos',[550 450 700 560]);
    bdec = barh([par_iter:-1:1], [X_m(out_ind(i),:,1)',X_f(out_ind(i),:,1)']);
    hold on
    binc = barh([par_iter:-1:1], [X_m(out_ind(i),:,2)',X_f(out_ind(i),:,2)']);
    
%     ax = gca;
%     ax.YTickLabel = in_label;
    
%     set(gca,'yticklabel',ycat)
%     ytickangle(45)
    
    yticklabels(in_label)
    ylim([1-0.5,par_iter+0.5])

    bdec(1).FaceColor = [0.203, 0.592, 0.835]; bdec(2).FaceColor = [0.835, 0.203, 0.576]; 
    binc(1).FaceColor = [0.203, 0.592, 0.835]; binc(2).FaceColor = [0.835, 0.203, 0.576];
    bdec(1).FaceAlpha = 0.5; bdec(2).FaceAlpha = 0.5;
    
    xlabel(strcat(out_label(i),' Percent change'))
    
    bars = get(gca, 'Children');
    legend(bars([4,2,3,1]), {'Male Dec','Male Inc','Female Dec','Female Inc'}, ...
           'FontSize',15,'Location','Southeast');

end

% Plot subplot. -----------------------------------------------------------

g = gobjects(2,1);

g(1) = figure('DefaultAxesFontSize',10);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.15, 4.5]);
s1(1) = subplot(1,2,1); 
s1(2) = subplot(1,2,2); 
s1(1).Position = s1(1).Position + [0.14, 0, 0.00, 0];
s1(2).Position = s1(2).Position + [0.06, 0, 0.00, 0];

bdec1 = barh(s1(1), [par_iter:-1:1], [X_m(out_ind(2),:,1)',X_f(out_ind(2),:,1)']);
hold(s1(1), 'on')
binc1 = barh(s1(1), [par_iter:-1:1], [X_m(out_ind(2),:,2)',X_f(out_ind(2),:,2)']);
% ax = gca;
% ax.YTickLabel = in_label;
yticklabels(s1(1), in_label)
ylim(s1(1), [1-0.5,par_iter+0.5])
bdec1(1).FaceColor = [0.203, 0.592, 0.835]; bdec1(2).FaceColor = [0.835, 0.203, 0.576]; 
binc1(1).FaceColor = [0.203, 0.592, 0.835]; binc1(2).FaceColor = [0.835, 0.203, 0.576];
bdec1(1).FaceAlpha = 0.5; bdec1(2).FaceAlpha = 0.5;
xlabel(s1(1), strcat(out_label(2),' Percent change'))
hold(s1(1), 'off')
title(s1(1), 'A')

bdec2 = barh(s1(2), [par_iter:-1:1], [X_m(out_ind(1),:,1)',X_f(out_ind(1),:,1)']);
hold(s1(2), 'on')
binc2 = barh(s1(2), [par_iter:-1:1], [X_m(out_ind(1),:,2)',X_f(out_ind(1),:,2)']);
set(gca,'YTick',[]);
ylim(s1(2), [1-0.5,par_iter+0.5])
bdec2(1).FaceColor = [0.203, 0.592, 0.835]; bdec2(2).FaceColor = [0.835, 0.203, 0.576]; 
binc2(1).FaceColor = [0.203, 0.592, 0.835]; binc2(2).FaceColor = [0.835, 0.203, 0.576];
bdec2(1).FaceAlpha = 0.5; bdec2(2).FaceAlpha = 0.5;
xlabel(s1(2), strcat(out_label(1),' Percent change'))
bars2 = get(gca, 'Children');
legend(s1(2), bars2([4,2,3,1]), {'Male Dec','Male Inc','Female Dec','Female Inc'}, ...
       'FontSize',7,'Position',[0.52 0.19 0.0 0.0]);
hold(s1(2), 'off')
title(s1(2), 'B')

% ---

g(2) = figure('DefaultAxesFontSize',10);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.15, 4.5]);
s2(1) = subplot(1,2,1); 
s2(2) = subplot(1,2,2); 
s2(1).Position = s2(1).Position + [0.14, 0, 0.00, 0];
s2(2).Position = s2(2).Position + [0.06, 0, 0.00, 0];

bdec1 = barh(s2(1), [par_iter:-1:1], [X_m(out_ind(4),:,1)',X_f(out_ind(4),:,1)']);
hold(s2(1), 'on')
binc1 = barh(s2(1), [par_iter:-1:1], [X_m(out_ind(4),:,2)',X_f(out_ind(4),:,2)']);
% ax = gca;
% ax.YTickLabel = in_label;
yticklabels(s2(1), in_label)
ylim(s2(1), [1-0.5,par_iter+0.5])
bdec1(1).FaceColor = [0.203, 0.592, 0.835]; bdec1(2).FaceColor = [0.835, 0.203, 0.576]; 
binc1(1).FaceColor = [0.203, 0.592, 0.835]; binc1(2).FaceColor = [0.835, 0.203, 0.576];
bdec1(1).FaceAlpha = 0.5; bdec1(2).FaceAlpha = 0.5;
xlabel(s2(1), strcat(out_label(4),' Percent change'))
hold(s2(1), 'off')
title(s2(1), 'A')

bdec2 = barh(s2(2), [par_iter:-1:1], [X_m(out_ind(3),:,1)',X_f(out_ind(3),:,1)']);
hold(s2(2), 'on')
binc2 = barh(s2(2), [par_iter:-1:1], [X_m(out_ind(3),:,2)',X_f(out_ind(3),:,2)']);
set(gca,'YTick',[]);
ylim(s2(2), [1-0.5,par_iter+0.5])
bdec2(1).FaceColor = [0.203, 0.592, 0.835]; bdec2(2).FaceColor = [0.835, 0.203, 0.576]; 
binc2(1).FaceColor = [0.203, 0.592, 0.835]; binc2(2).FaceColor = [0.835, 0.203, 0.576];
bdec2(1).FaceAlpha = 0.5; bdec2(2).FaceAlpha = 0.5;
xlabel(s2(2), strcat(out_label(3),' Percent change'))
bars2 = get(gca, 'Children');
legend(s2(2), bars2([4,2,3,1]), {'Male Dec','Male Inc','Female Dec','Female Inc'}, ...
       'FontSize',7,'Position',[0.52 0.845 0.0 0.0]);
hold(s2(2), 'off')
title(s2(2), 'B')

% % Save figures ------------------------------------------------------------
% 
% save_data_name = sprintf('quant_of_int_sen_analTEST.fig' );
% save_data_name = strcat('Figures/', save_data_name);
% savefig([f;g], save_data_name)

toc

end





























