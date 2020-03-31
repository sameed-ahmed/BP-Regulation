% This simulates the blood pressure regulation model bp_reg.m for drug administration.
% 
% Steady state data is calculated by solve_ss_scenario.m.

function solve_ss_drugs_dist

close all

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

%% Variable names for plotting.
var_names  = {'$rsna$'; '$\alpha_{map}$'; '$\alpha_{rap}$'; '$R_{r}$'; ...
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
              '$\Phi_{u}$'; '$\Phi_{win}$'; '$R_{ea}/R_{r}$'};
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Begin user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Physiological scenarios
% Normal  - Normal conditions
% m_RAS   - male RAS pars
% m_Reab  - male fractional sodium and water reabsorption
% Pri_Hyp - essential/primary hypertension
scenario1 = {'Normal', 'm_RSNA', 'm_AT2R', 'm_RAS', 'm_Reab', ...
             'm_RAS_m_Reab', 'm_RSNA_m_Reab'};
% scenario1 = {'Normal', 'm_RSNA', 'm_AT2R', 'm_RAS', 'm_Reab', ...
%              'm_RAS_m_Reab', 'm_RSNA_m_Reab', ...
%              'Pri_Hyp'};
fixed_ss1 = 1;
num_scen = length(scenario1);
% Drug scenarios
% Normal - Normal conditions
% AngII  - Ang II infusion fmol/(ml min)
% ACEi   - Angiotensin converting enzyme inhibitor %
% ARB1   - Angiotensin receptor 1 blocker %
% ARB2   - Angiotensin receptor 2 blocker %
% DRI    - Direct renin inhibitor %
% MRB    - Aldosterone blocker (MR?) %
% RSS    - Renin secretion stimulator (thiazide?) % % NOT COMPLETE
scenario2 = {'Normal', 'AngII', 'ACEi', 'ARB1', 'ARB2', 'DRI', 'MRB', 'RSS'};
fixed_ss2 = [3];


% Species
spe_ind = 2;

% Bootstrap replicate sample number
% sample_num = random('Discrete Uniform',1000)
% sample_num = 208
% sample_num = 655
% fixed_sample = 655;
num_samples = 1000;
% num_samples = 5;
fixed_sample = 1;

% % Number of intervals for dose
% num_doses = 21;
% % Range for dose reponse
% drug_dose = linspace(0,1.00,num_doses);
% % drug_dose = linspace(0,0.95,num_iter);
% inhibit = 0.95;
% fixed_dose = round((inhibit - 0) / ((1-0)/(21-1)));

% Drug dose
drug_dose = 0.95;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

species = {'human', 'rat'   };
sex     = {'male' , 'female'};

% Number of variables
num_vars = 93;

% Initialize variables.
% X = (variable, sample, sex, scenario)
X_ss = zeros(num_vars,num_samples,2,num_scen);
X_bl = zeros(num_vars,num_samples,2,num_scen);

for sce_ind = fixed_ss1:fixed_ss1 % scenario
for sex_ind = 1:2        % sex

%% Load bootstrap replicate parameters & variables created by create_par_bs_rep.m.

% Parameters
load_data_name_pars = sprintf('%s_%s_pars_scenario_Pri_Hyp_bs_rep1000.mat', ...
                              species{spe_ind},sex{sex_ind});
load(load_data_name_pars, 'pars_rep');
num_pars = size(pars_rep,1);

% Variables
load_data_name_vars = sprintf('%s_%s_ss_data_scenario_Pri_Hyp_bs_rep1000.mat', ...
                              species{spe_ind},sex{sex_ind});
load(load_data_name_vars, 'SSdata_rep');
num_vars = size(SSdata_rep,1);
% Store baseline value to computer relative change.
X_bl(:,:,sex_ind,sce_ind) = SSdata_rep(:,1:num_samples);
% X_bl(:,:,sex_ind,sce_ind) = SSdata_rep(:,fixed_sample);

for sam_iter = 1:num_samples % samples
% for sam_iter = fixed_sample:fixed_sample % samples

%% Drugs

varargin_input = {scenario1{sce_ind},true};

for i = 1:length(fixed_ss2)
    if     strcmp(scenario2{fixed_ss2(i)}, 'AngII')
        if     strcmp(sex{sex_ind}, 'male')
            varargin_input = [varargin_input, 'AngII',910]; % Sullivan 2010
        elseif strcmp(sex{sex_ind}, 'female')
            varargin_input = [varargin_input, 'AngII',505]; % Sullivan 2010
        end
    elseif strcmp(scenario2{fixed_ss2(i)}, 'ACEi' )
            varargin_input = [varargin_input, 'ACEi' ,drug_dose]; % 
    elseif strcmp(scenario2{fixed_ss2(i)}, 'ARB1' )
            varargin_input = [varargin_input, 'ARB1' ,drug_dose]; % 
    elseif strcmp(scenario2{fixed_ss2(i)}, 'ARB2' )
            varargin_input = [varargin_input, 'ARB2' ,drug_dose]; % 
    elseif strcmp(scenario2{fixed_ss2(i)}, 'DRI'  )
            varargin_input = [varargin_input, 'DRI'  ,drug_dose]; % 
    elseif strcmp(scenario2{fixed_ss2(i)}, 'MRB'  )
            varargin_input = [varargin_input, 'MRB'  ,drug_dose]; % 
    elseif strcmp(scenario2{fixed_ss2(i)}, 'RSS'  )
            varargin_input = [varargin_input, 'RSS'  ,drug_dose]; % 
    end
end

%% Solve system steady state

%% Initialization

% Initial guess for the variables.
% Find the steady state solution, so the derivative is 0.
% Arbitrary value for time to input, sufficiently greater than tchange so
% that drug dose reaches steady state value (tanh).
x0_ss = SSdata_rep(:,sam_iter); x_p0_ss = zeros(num_vars,1); t_ss = 2000;

% Time at which to change place holder.
tchange_ss = 0;

% Solver options
options_ss = optimset('Display','off');
% Solve system
[SSdata, residual, ...
 exitflag, output] = fsolve(@(x) ...
                            bp_reg_mod(t_ss,x,x_p0_ss,...
                                       pars_rep(:,sam_iter),tchange_ss,...
                                       varargin_input{:}), ...
                            x0_ss, options_ss);

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

% Store solution.
% X = (variable, sample, sex, scenario)
X_ss(:,sam_iter,sex_ind,sce_ind) = SSdata;

% Sanity check to see script's progress. Also a check for where to
% troubleshoot in case the solver does not converge.
fprintf('%s %s sample = %s out of %s \n', ...
        scenario1{sce_ind},sex{sex_ind},num2str(sam_iter),num2str(num_samples))

end % samples
end % sex
end % scenario

%% Plot

% Retrieve male and female.
% X_m/f = (variable, sample, scenario)
X_ss_m = reshape(X_ss(:,:,1,:), [num_vars,num_samples,num_scen]); 
X_ss_f = reshape(X_ss(:,:,2,:), [num_vars,num_samples,num_scen]); 
X_bl_m = reshape(X_bl(:,:,1,:), [num_vars,num_samples,num_scen]); 
X_bl_f = reshape(X_bl(:,:,2,:), [num_vars,num_samples,num_scen]); 

% Compute relative change.
X_rel_m = (X_ss_m(:,:,fixed_ss1) - X_bl_m(:,:,fixed_ss1)) ...
          ./ X_bl_m(:,:,fixed_ss1) * 100;
X_rel_f = (X_ss_f(:,:,fixed_ss1) - X_bl_f(:,:,fixed_ss1)) ...
          ./ X_bl_f(:,:,fixed_ss1) * 100;

% X_bl_m(7,:,fixed_ss1)'
% X_ss_m(7,:,fixed_ss1)'
% X_bl_f(7,:,fixed_ss1)'
% X_ss_f(7,:,fixed_ss1)'

% x-axis
xscale = drug_dose * 100;
xlower = drug_dose(1); xupper = drug_dose(end); 

% y-axis limits
ylower = zeros(num_vars,1); yupper = ylower; 
for i = 1:length(ylower)
    ylower(i) = 0.95*min( min(X_ss_m(i,:,fixed_ss1)), min(X_ss_f(i,:,fixed_ss1)) );
    yupper(i) = 1.05*max( max(X_ss_m(i,:,fixed_ss1)), max(X_ss_f(i,:,fixed_ss1)) );
    if abs(yupper(i)) < eps*100
        ylower(i) = -10^(-5); yupper(i) = 10^(-5);
    end
end

%% Plot all vars distribution. --------------------------------------------

f1 = gobjects(7,1);
s1 = gobjects(7,15);
% Loop through each set of subplots.
for i = 1:7
    f1(i) = figure('pos',[750 500 650 450]);
    % This is to avoid the empty plots in the last subplot set.
    if i == 7
        last_plot = mod(num_vars, 15);
    else
        last_plot = 15;
    end
    % Loop through each subplot within a set of subplots.
    for j = 1:last_plot
        s1(i,j) = subplot(3,5,j);
%         s(i,j).Position = s(i,j).Position + [0 0 0.01 0];
        
%         h1 = histogram(s1(i,j),X_ss_m((i-1)*15 + j,:,fixed_ss1));
        h1 = histogram(s1(i,j),X_rel_m((i-1)*15 + j,:),10);
%         hold on
%         h2 = histogram(s2(i,j),X_ss_f((i-1)*15 + j,fixed_dose,:,fixed_ss1), 3);
%         hold off
        
        h1.Normalization = 'probability'; 
%         h1.BinWidth = 1.0; h2.BinWidth = 1.0; 
        h1.FaceColor = [0.203, 0.592, 0.835];

%         xlabel_name = strcat(vars_names((i-1)*15 + j), ' (', num2str(pars_hyp_bl(i),3), pars_units(i), ')');
        xlabel_name = strcat(var_names((i-1)*15 + j));
        xlabel(s1(i,j), xlabel_name, 'Interpreter','latex', 'FontSize',16)

%         legend('Male', 'Female')
    end
    hist_title = sprintf('Male, ACEi %s%%',num2str(drug_dose*100));
    sgtitle(hist_title, 'FontSize',14)
end

f2 = gobjects(7,1);
s2 = gobjects(7,15);
% Loop through each set of subplots.
for i = 1:7
    f2(i) = figure('pos',[750 500 650 450]);
    % This is to avoid the empty plots in the last subplot set.
    if i == 7
        last_plot = mod(num_vars, 15);
    else
        last_plot = 15;
    end
    % Loop through each subplot within a set of subplots.
    for j = 1:last_plot
        s2(i,j) = subplot(3,5,j);
%         s(i,j).Position = s(i,j).Position + [0 0 0.01 0];
        
%         h1 = histogram(s2(i,j),X_ss_f((i-1)*15 + j,:,fixed_ss1));
        h1 = histogram(s2(i,j),X_rel_f((i-1)*15 + j,:),10);
%         hold on
%         h2 = histogram(s2(i,j),X_ss_f((i-1)*15 + j,fixed_dose,:,fixed_ss1), 3);
%         hold off
        
        h1.Normalization = 'probability'; 
%         h1.BinWidth = 1.0; h2.BinWidth = 1.0; 
        h1.FaceColor = [0.835, 0.203, 0.576];

%         xlabel_name = strcat(vars_names((i-1)*15 + j), ' (', num2str(pars_hyp_bl(i),3), pars_units(i), ')');
        xlabel_name = strcat(var_names((i-1)*15 + j));
        xlabel(s2(i,j), xlabel_name, 'Interpreter','latex', 'FontSize',16)

%         legend('Male', 'Female')
    end
    hist_title = sprintf('Female, ACEi %s%%',num2str(drug_dose*100));
    sgtitle(hist_title, 'FontSize',14)
end

%% Plot interesting variables. --------------------------------------------

% Interesting variables to plot.
var_ind = [30;33;41;42;9;73;74;6;7;92;69;70]; num_vars_sub = length(var_ind);

f3 = figure('pos',[000 000 600 600], 'DefaultAxesFontSize',12);
s3 = gobjects(1,num_vars_sub);
% Loop through each subplot within a set of subplots.
for j = 1:num_vars_sub
    s3(j) = subplot(4,3,j);
%     if     mod(j,3) == 1
%         hshift = -0.05;
%     elseif mod(j,3) == 0
%         hshift = 0.05;
%     else
%         hshift = 0;
%     end
%     s3(j).Position = s3(j).Position + [hshift 0 0.01 0.01];
    
        h1 = histogram(s3(j),X_rel_m(var_ind(j),:),10);
%         hold(s3(j), 'on')
%         h2 = histogram(s3(j),X_rel_f(var_ind(j),:),10);
%         hold(s3(j), 'off')
        
        h1.Normalization = 'probability'; %h2.Normalization = 'probability'; 
%         h1.BinWidth = 1.0; h2.BinWidth = 1.0; 
        h1.FaceColor = [0.203, 0.592, 0.835]; %h2.FaceColor = [0.835, 0.203, 0.576];

%     plot(s3(j), xscale,X_ss_m(var_ind(j),:,fixed_ss1), 'Color',[0.203, 0.592, 0.835], 'LineWidth',2.5);
%     hold(s3(j), 'on')
%     plot(s3(j), xscale,X_ss_f(var_ind(j),:,fixed_ss1), 'Color',[0.835, 0.203, 0.576], 'LineWidth',2.5);
%     hold(s3(j), 'off')

    xlabel(s3(j), var_names(var_ind(j)), 'Interpreter','latex', 'FontSize',16)
end
% legend(s3(1),'Male','Female', 'Location','east')
hist_title = sprintf('Male, ACEi %s%%',num2str(drug_dose*100));
sgtitle(hist_title, 'FontSize',14)

f4 = figure('pos',[000 000 600 600], 'DefaultAxesFontSize',12);
s4 = gobjects(1,num_vars_sub);
% Loop through each subplot within a set of subplots.
for j = 1:num_vars_sub
    s4(j) = subplot(4,3,j);
%     if     mod(j,3) == 1
%         hshift = -0.05;
%     elseif mod(j,3) == 0
%         hshift = 0.05;
%     else
%         hshift = 0;
%     end
%     s4(j).Position = s4(j).Position + [hshift 0 0.01 0.01];
    
        h1 = histogram(s4(j),X_rel_f(var_ind(j),:),10);
%         hold(s4(j), 'on')
%         h2 = histogram(s4(j),X_rel_f(var_ind(j),:),10);
%         hold(s4(j), 'off')
        
        h1.Normalization = 'probability'; %h2.Normalization = 'probability'; 
%         h1.BinWidth = 1.0; h2.BinWidth = 1.0; 
        h1.FaceColor = [0.835, 0.203, 0.576];

%     plot(s4(j), xscale,X_ss_m(var_ind(j),:,fixed_ss1), 'Color',[0.203, 0.592, 0.835], 'LineWidth',2.5);
%     hold(s4(j), 'on')
%     plot(s4(j), xscale,X_ss_f(var_ind(j),:,fixed_ss1), 'Color',[0.835, 0.203, 0.576], 'LineWidth',2.5);
%     hold(s4(j), 'off')

    xlabel(s4(j), var_names(var_ind(j)), 'Interpreter','latex', 'FontSize',16)
end
% legend(s4(1),'Male','Female', 'Location','east')
hist_title = sprintf('Female, ACEi %s%%',num2str(drug_dose*100));
sgtitle(hist_title, 'FontSize',14)

%% Plot variables that explain mechanims. ---------------------------------

R_ss_m = reshape(X_ss_m(74,:,:) ./ X_ss_m(4,:,:), [num_samples,num_scen]);
R_ss_f = reshape(X_ss_f(74,:,:) ./ X_ss_f(4,:,:), [num_samples,num_scen]);
R_bl_m = reshape(X_bl_m(74,:,:) ./ X_bl_m(4,:,:), [num_samples,num_scen]);
R_bl_f = reshape(X_bl_f(74,:,:) ./ X_bl_f(4,:,:), [num_samples,num_scen]);

R_rel_m = (R_ss_m(:,fixed_ss1) - R_bl_m(:,fixed_ss1)) ...
          ./ R_bl_m(:,fixed_ss1) * 100;
R_rel_f = (R_ss_f(:,fixed_ss1) - R_bl_f(:,fixed_ss1)) ...
          ./ R_bl_f(:,fixed_ss1) * 100;

FRNA_ss_m = reshape((X_ss_m(11,:,:) - X_ss_m(27,:,:)) ./ X_ss_m(11,:,:), [num_samples,num_scen]) * 100;
FRNA_ss_f = reshape((X_ss_f(11,:,:) - X_ss_f(27,:,:)) ./ X_ss_f(11,:,:), [num_samples,num_scen]) * 100;
FRNA_bl_m = reshape((X_bl_m(11,:,:) - X_bl_m(27,:,:)) ./ X_bl_m(11,:,:), [num_samples,num_scen]) * 100;
FRNA_bl_f = reshape((X_bl_f(11,:,:) - X_bl_f(27,:,:)) ./ X_bl_f(11,:,:), [num_samples,num_scen]) * 100;

FRNA_rel_m = (FRNA_ss_m(:,fixed_ss1) - FRNA_bl_m(:,fixed_ss1)) ...
             ./ FRNA_bl_m(:,fixed_ss1) * 100;
FRNA_rel_f = (FRNA_ss_f(:,fixed_ss1) - FRNA_bl_f(:,fixed_ss1)) ...
             ./ FRNA_bl_f(:,fixed_ss1) * 100;

FRW_ss_m = reshape((X_ss_m( 7,:,:) - X_ss_m(92,:,:)) ./ X_ss_m( 7,:,:), [num_samples,num_scen]) * 100;
FRW_ss_f = reshape((X_ss_f( 7,:,:) - X_ss_f(92,:,:)) ./ X_ss_f( 7,:,:), [num_samples,num_scen]) * 100;
FRW_bl_m = reshape((X_bl_m( 7,:,:) - X_bl_m(92,:,:)) ./ X_bl_m( 7,:,:), [num_samples,num_scen]) * 100;
FRW_bl_f = reshape((X_bl_f( 7,:,:) - X_bl_f(92,:,:)) ./ X_bl_f( 7,:,:), [num_samples,num_scen]) * 100;

FRW_rel_m = (FRW_ss_m(:,fixed_ss1) - FRW_bl_m(:,fixed_ss1)) ...
            ./ FRW_bl_m(:,fixed_ss1) * 100;
FRW_rel_f = (FRW_ss_f(:,fixed_ss1) - FRW_bl_f(:,fixed_ss1)) ...
            ./ FRW_bl_f(:,fixed_ss1) * 100;

g1 = figure('DefaultAxesFontSize',14);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 12, 4]);
s_mech(1) = subplot(1,3,1); 
s_mech(2) = subplot(1,3,2); 
s_mech(3) = subplot(1,3,3); 

h1 = histogram(s_mech(1),R_rel_m(:),10);
hold(s_mech(1), 'on')
h2 = histogram(s_mech(1),R_rel_f(:),10);
hold(s_mech(1), 'off')
h1.Normalization = 'probability'; h2.Normalization = 'probability';  
h1.BinWidth = 5.0; h2.BinWidth = 5.0; 
h1.FaceColor = [0.203, 0.592, 0.835]; h2.FaceColor = [0.835, 0.203, 0.576];
legend(s_mech(1), [h1, h2],{'Male','Female'}, 'FontSize',10,'Location','Northeast');
xlabel(s_mech(1), 'R_{EA}/R_R');
title(s_mech(1), 'A')

h1 = histogram(s_mech(2),FRNA_rel_m(:),10);
hold(s_mech(2), 'on')
h2 = histogram(s_mech(2),FRNA_rel_f(:),10);
hold(s_mech(2), 'off')
h1.Normalization = 'probability'; h2.Normalization = 'probability';  
h1.BinWidth = 0.05; h2.BinWidth = 0.05; 
h1.FaceColor = [0.203, 0.592, 0.835]; h2.FaceColor = [0.835, 0.203, 0.576];
xlabel(s_mech(2), 'FR_{Na^+}');
title(s_mech(2), 'B')

h1 = histogram(s_mech(3),FRW_rel_m(:),10);
hold(s_mech(3), 'on')
h2 = histogram(s_mech(3),FRW_rel_f(:),10);
hold(s_mech(3), 'off')
h1.Normalization = 'probability'; h2.Normalization = 'probability';  
h1.BinWidth = 0.05; h2.BinWidth = 0.05; 
h1.FaceColor = [0.203, 0.592, 0.835]; h2.FaceColor = [0.835, 0.203, 0.576];
xlabel(s_mech(3), 'FR_{U}');
title(s_mech(3), 'C')

%% Plot Mean Arterial Pressure distribution. ------------------------------

% Actual, change, and % change in MAP.
% X_m/f = (variable, sample, scenario)
MAP_ac_m = zeros(num_samples,num_scen); MAP_ac_f = zeros(num_samples,num_scen);
MAP_ch_m = zeros(num_samples,num_scen); MAP_ch_f = zeros(num_samples,num_scen);
MAP_pc_m = zeros(num_samples,num_scen); MAP_pc_f = zeros(num_samples,num_scen);
for i = 1:num_scen
    MAP_ac_m(:,i) = (X_ss_m(42,:,i)                 )                        ;
    MAP_ac_f(:,i) = (X_ss_f(42,:,i)                 )                        ;
    
    MAP_ch_m(:,i) = (X_ss_m(42,:,i) - X_bl_m(42,:,i))                        ;
    MAP_ch_f(:,i) = (X_ss_f(42,:,i) - X_bl_f(42,:,i))                        ;
    
    MAP_pc_m(:,i) = (X_ss_m(42,:,i) - X_bl_m(42,:,i)) ./ X_bl_m(42,1,i) * 100;
    MAP_pc_f(:,i) = (X_ss_f(42,:,i) - X_bl_f(42,:,i)) ./ X_bl_f(42,1,i) * 100;
end

g2 = figure('DefaultAxesFontSize',14);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 12, 4]);
s_map1(1) = subplot(1,3,1); 
s_map1(2) = subplot(1,3,2); 
s_map1(3) = subplot(1,3,3);

h1 = histogram(s_map1(1),MAP_ac_m(:,fixed_ss1),10);
hold(s_map1(1), 'on')
h2 = histogram(s_map1(1),MAP_ac_f(:,fixed_ss1),10);
hold(s_map1(1), 'off')
h1.Normalization = 'probability'; h2.Normalization = 'probability';  
h1.BinWidth = 5.0; h2.BinWidth = 5.0; 
h1.FaceColor = [0.203, 0.592, 0.835]; h2.FaceColor = [0.835, 0.203, 0.576];
legend(s_map1(1), [h1, h2],{'Male','Female'}, 'FontSize',10,'Location','Northwest');
xlabel(s_map1(1), 'MAP (mmHg)');
title(s_map1(1), 'A')

h1 = histogram(s_map1(2),MAP_ch_m(:,fixed_ss1),10);
hold(s_map1(2), 'on')
h2 = histogram(s_map1(2),MAP_ch_f(:,fixed_ss1),10);
hold(s_map1(2), 'off')
h1.Normalization = 'probability'; h2.Normalization = 'probability';  
h1.BinWidth = 5.0; h2.BinWidth = 5.0; 
h1.FaceColor = [0.203, 0.592, 0.835]; h2.FaceColor = [0.835, 0.203, 0.576];
xlabel(s_map1(2), '\DeltaMAP (mmHg)');
title(s_map1(2), 'B')

h1 = histogram(s_map1(3),MAP_pc_m(:,fixed_ss1),10);
hold(s_map1(3), 'on')
h2 = histogram(s_map1(3),MAP_pc_f(:,fixed_ss1),10);
hold(s_map1(3), 'off')
h1.Normalization = 'probability'; h2.Normalization = 'probability';  
h1.BinWidth = 5.0; h2.BinWidth = 5.0; 
h1.FaceColor = [0.203, 0.592, 0.835]; h2.FaceColor = [0.835, 0.203, 0.576];
xlabel(s_map1(3), '% \DeltaMAP');
title(s_map1(3), 'C')

%% Save figures and data. -------------------------------------------------
 
save_data_name = sprintf('dose_distribution_%s%s%%.fig', ...
                         scenario2{fixed_ss2},num2str(drug_dose*100));
save_data_name = strcat('Figures/', save_data_name);
savefig([f1;f2;f3;f4;g1;g2], save_data_name)

save_data_name = sprintf('%s_male_ss_data_scenario_Pri_Hyp_%s%s%%.mat'  , ...
                         species{spe_ind},scenario2{fixed_ss2},num2str(drug_dose*100));
save_data_name = strcat('Data/', save_data_name);
save(save_data_name, 'X_ss_m', 'X_bl_m')

save_data_name = sprintf('%s_female_ss_data_scenario_Pri_Hyp_%s%s%%.mat', ...
                         species{spe_ind},scenario2{fixed_ss2},num2str(drug_dose*100));
save_data_name = strcat('Data/', save_data_name);
save(save_data_name, 'X_ss_f', 'X_bl_f')

end






























