% This script performs local sensitivty analysis by perturbing parameters 
% and then calculating the steady state values of the blood pressure 
% regulation model bp_reg_solve_sen_anal.m for each inputted value. 
% 
% Steady state data for the intial guess is calculated by solve_ss_baseline.m.
% 
% Certain variables are then plotted as relative change versus each
% perturbed parameter.

function solve_ss_sen_anal

close all

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Begin user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input of interest
% N_rsna       - par 1 ; -30%, +30%
% N_als        - par 19; -30%, +30%
% N_rs         - par 22; -30%, +30%
% R_aass       - par 2 ; -50%, +250%
% R_eass       - par 3 ; -50%, +250%
% c_ACE        - par 33; -50%, +50%
% c_AT1R       - par 38; -50%, +50%
% c_AT2R       - par 39; -50%, +50%
% % Psi_AT2RAA - par 42; -50%, +50%
% % Psi_AT2REA - par 43; -50%, +50%
% Parameter indices
in_ind      = [1;19;22;2;3;33;38;39];
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
% Phi_u             - var 63
% R_aa              - var 86
% R_ea              - var 87
% eta_alpha-sodreab - var 13, 19, 23 (alpha = pt, dt, cd)
% eta_alpha-wreab   - var 52, 56, 60 (alpha = pt, dt, cd)
out_ind = [7;42;63;86;87;13;19;23;52;56;60];
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
% AT2R-  - Block AT2R through decay
% RHyp   - Renal hypertension due to increased afferent arteriolar resistance
scenario = {'Normal', 'AngII', 'ACEi', 'ARB', 'AT2R-', 'RHyp'};
% Index of scenario to plot for all variables
fixed_ss = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of variables
num_vars = 92;

% Initialize variables.
% X = (variable, parameter, change, gender)
X = zeros(num_vars,par_iter,2,2);

gender = {'male'    , 'female'  };
change = {'decrease', 'increase'};

for iter = 1:par_iter % par
for cc   = 1:2        % change
for gg   = 1:2        % gender

%% Parameters

% Parameter input
pars = get_pars(gender{gg}, scenario{fixed_ss});

% Change parameter of interest.
if     strcmp(change{cc}, 'decrease')
    pars(in_ind(iter)) = pars(in_ind(iter)) * (1 - range_lower(iter)/100);
elseif strcmp(change{cc}, 'increase')
    pars(in_ind(iter)) = pars(in_ind(iter)) * (1 + range_upper(iter)/100);
end

%% Drugs

% drugs = [Ang II inf rate fmol/(ml min), ACEi target level, ARB target level, AT2R decay rate]
if     strcmp(scenario{fixed_ss}, 'Normal') || strcmp(scenario{fixed_ss}, 'RHyp')
    drugs = [0, 0, 0, 0];
elseif strcmp(scenario{fixed_ss}, 'AngII')
    if     strcmp(gender{gg}, 'male')
        drugs = [2022, 0, 0, 0]; % Sampson 2008
    elseif strcmp(gender{gg}, 'female')
        drugs = [2060, 0, 0, 0]; % Sampson 2008
    end
elseif strcmp(scenario{fixed_ss}, 'ACEi')
    drugs = [0, 0.78, 0, 0]; % Leete 2018
elseif strcmp(scenario{fixed_ss}, 'ARB')
    drugs = [0, 0, 0.67, 0]; % Leete 2018
elseif strcmp(scenario{fixed_ss}, 'AT2R-')
    drugs = [0, 0, 0, 10];
end

%% Variables initial guess

% Retrieve and replace parameters in fixed variable equations.
% Load data for steady state initial guess. 
% Set name for data file to be loaded based upon gender.    
load_data_name = sprintf('%s_ss_data_scenario_Normal.mat', gender{gg});
load(load_data_name, 'SSdata');
fixed_ind = [2, 10, 14, 24, 44, 49, 66, 71, 88];
fixed_var_pars = SSdata(fixed_ind);
phicophico = SSdata(33); cadhcadh = SSdata(47);
fixed_var_pars = [fixed_var_pars; cadhcadh; phicophico];
SSdata(fixed_ind) = 1;
SSdataIG = SSdata;
clear SSdata

% Ang II infusion steady state value is sensitive to such a huge change.
% Thus it is done incrementally and iteratively.
if     strcmp(scenario{fixed_ss}, 'AngII')
    load_data_name = sprintf('%s_ss_data_scenario_AngII.mat', gender{gg});
    load(load_data_name, 'SSdata');
    SSdataIG = SSdata;
    clear SSdata;
end

% Initial guess for the variables.
% Find the steady state solution, so the derivative is 0.
% Arbitrary value for time to input.
x0 = SSdataIG; x_p0 = zeros(num_vars,1); t = 0;

%% Find steady state solution

% options = optimset(); %options = optimset('MaxFunEvals',num_vars*100+10000);
options = optimset('Display','off');
[SSdata, residual, ...
 exitflag, output] = fsolve(@(x) bp_reg_solve_sen_anal(t,x,x_p0,pars ,...
                                                       fixed_var_pars,...
                                                       drugs)        ,...
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
X(:,iter,cc,gg) = SSdata;

% % Sanity check to see script's progress. Also a check for where to
% % troubleshoot in case the solver does not converge.
% if     strcmp(gender{gg}, 'male')
%     fprintf('  male %s iteration = %s out of %s \n', ...
%             change{cc},num2str(iter),num2str(par_iter))
% elseif strcmp(gender{gg}, 'female')
%     fprintf('female %s iteration = %s out of %s \n', ...
%             change{cc},num2str(iter),num2str(par_iter))
% end

end % gender
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
    
%     set(gca,'yticklabel',ycat)
    yticklabels(in_label)
%     ytickangle(45)
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
       'FontSize',7,'Location','Southeast');
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
       'FontSize',7,'Location','Southeast');
hold(s2(2), 'off')
title(s2(2), 'B')

% Save figures ------------------------------------------------------------

% save_data_name = sprintf('quant_of_int_sen_anal.fig' );
% save_data_name = strcat('Figures/', save_data_name);
% savefig([f;g], save_data_name)

end





























