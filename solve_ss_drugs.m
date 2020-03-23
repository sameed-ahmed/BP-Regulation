% This simulates the blood pressure regulation model bp_reg.m for drug administration.
% 
% Steady state data is calculated by solve_ss_scenario.m.

function solve_ss_drugs

close all

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

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
sample_num = 655

% Number of intervals for dose
num_iter = 21;
% drug_dose = linspace(0,1.00,num_iter);
drug_dose = linspace(0,0.95,num_iter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

species = {'human', 'rat'   };
sex     = {'male' , 'female'};

% Number of variables
num_vars = 93;

% Initialize variables.
% X = (variable, iteration, sex, scenario)
X_ss = zeros(num_vars,num_iter,2,num_scen);

for sce_ind = fixed_ss1:fixed_ss1 % scenario
for sex_ind = 1:2        % sex

%% Load bootstrap replicate parameters & variables created by create_par_bs_rep.m.

% Parameters
load_data_name_pars = sprintf('%s_%s_pars_scenario_Pri_Hyp_bs_rep1000.mat', ...
                              species{spe_ind},sex{sex_ind});
load(load_data_name_pars, 'pars_rep');
num_pars   = size(pars_rep,1);
num_sample = size(pars_rep,2);
pars_rep = pars_rep(:,sample_num);

% Variables
load_data_name_vars = sprintf('%s_%s_ss_data_scenario_Pri_Hyp_bs_rep1000.mat', ...
                              species{spe_ind},sex{sex_ind});
load(load_data_name_vars, 'SSdata_rep');
num_vars   = size(SSdata_rep,1);
SSdata_rep = SSdata_rep(:,sample_num);

for iter = 1:num_iter % range

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
            varargin_input = [varargin_input, 'ACEi' ,drug_dose(iter)]; % 
    elseif strcmp(scenario2{fixed_ss2(i)}, 'ARB1' )
            varargin_input = [varargin_input, 'ARB1' ,drug_dose(iter)]; % 
    elseif strcmp(scenario2{fixed_ss2(i)}, 'ARB2' )
            varargin_input = [varargin_input, 'ARB2' ,drug_dose(iter)]; % 
    elseif strcmp(scenario2{fixed_ss2(i)}, 'DRI'  )
            varargin_input = [varargin_input, 'DRI'  ,drug_dose(iter)]; % 
    elseif strcmp(scenario2{fixed_ss2(i)}, 'MRB'  )
            varargin_input = [varargin_input, 'MRB'  ,drug_dose(iter)]; % 
    elseif strcmp(scenario2{fixed_ss2(i)}, 'RSS'  )
            varargin_input = [varargin_input, 'RSS'  ,drug_dose(iter)]; % 
    end
end

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

%% Solve system steady state

% Initial guess for the variables.
% Find the steady state solution, so the derivative is 0.
% Arbitrary value for time to input, sufficiently greater than tchange.
x0_ss = SSdata_rep; x_p0_ss = zeros(num_vars,1); t_ss = 2000;

% Time at which to change place holder.
tchange_ss = 0;

% Solver options
options_ss = optimset('Display','off');
% Solve system
[SSdata, residual, ...
 exitflag, output] = fsolve(@(x) ...
                            bp_reg_mod(t_ss,x,x_p0_ss,pars_rep,tchange_ss,varargin_input{:}), ...
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

X_ss(42,iter,sex_ind,sce_ind);
% Store solution.
% X = (variable, iteration, sex, scenario)
X_ss(:,iter,sex_ind,sce_ind) = SSdata;
X_ss(42,iter,sex_ind,sce_ind);

SSdata_rep = SSdata;

% Sanity check to see script's progress. Also a check for where to
% troubleshoot in case the solver does not converge.
fprintf('%s %s iteration = %s out of %s \n', ...
        scenario1{sce_ind},sex{sex_ind},num2str(iter),num2str(num_iter))

end % range
end % sex
end % scenario

%% Plot

% Retrieve male and female.
% X_m/f = (variable, iteration, scenario)
X_ss_m = reshape(X_ss(:,:,1,:), [num_vars,num_iter,num_scen]); 
X_ss_f = reshape(X_ss(:,:,2,:), [num_vars,num_iter,num_scen]); 

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

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

% Interesting variables to plot.
var_ind = [33;41;42;9;73;74;6;7;27;92;93;29]; num_vars_sub = length(var_ind);

% Plot all vars vs time. --------------------------------------------------

f1 = gobjects(7,1);
s = gobjects(7,15);
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
        s(i,j) = subplot(3,5,j);
%         s(i,j).Position = s(i,j).Position + [0 0 0.01 0];
        
        plot(s(i,j), xscale,X_ss_m((i-1)*15 + j,:,fixed_ss1),'b', ...
                     xscale,X_ss_f((i-1)*15 + j,:,fixed_ss1),'r');
        
%         xlim([xlower, xupper])
        ylim([ylower((i-1)*15 + j), yupper((i-1)*15 + j)])

        xlabel_name = strcat(scenario2{fixed_ss2}, ' %');
        xlabel(xlabel_name, 'FontSize',12)
        title(names((i-1)*15 + j), 'Interpreter','latex', 'FontSize',15)
%         legend('Male', 'Female')
    end
end

% Plot interesting variables. ---------------------------------------------

f2 = figure('pos',[000 000 600 600], 'DefaultAxesFontSize',12);
s2 = gobjects(1,num_vars_sub);
% Loop through each subplot within a set of subplots.
for j = 1:num_vars_sub
    s2(j) = subplot(4,3,j);
    if     mod(j,3) == 1
        hshift = -0.05;
    elseif mod(j,3) == 0
        hshift = 0.05;
    else
        hshift = 0;
    end
    s2(j).Position = s2(j).Position + [hshift 0 0.01 0.01];

    plot(s2(j), xscale,X_ss_m(var_ind(j),:,fixed_ss1), 'Color',[0.203, 0.592, 0.835], 'LineWidth',2.5);
    hold(s2(j), 'on')
    plot(s2(j), xscale,X_ss_f(var_ind(j),:,fixed_ss1), 'Color',[0.835, 0.203, 0.576], 'LineWidth',2.5);
    hold(s2(j), 'off')

%     xlim([xlower, xupper])
    ylim([ylower(var_ind(j)), yupper(var_ind(j))])
    
%     if j == 10 || j == 11
%         ax = gca;
%         ax.YAxis.Exponent = -3;
%     end

    ylabel(names(var_ind(j)), 'Interpreter','latex', 'FontSize',16)
end
legend(s2(1),'Male','Female', 'Location','east')
xlh = xlabel(s2(11),xlabel_name);
xlh.Position(2) = xlh.Position(2) - 0.0005;

% Plot Mean Arterial Pressure vs Dose. ------------------------------------

% Actual, change, and % change in MAP.
% X_m/f = (variable, points, scenario)
MAP_ac_m = zeros(num_iter,num_scen); MAP_ac_f = zeros(num_iter,num_scen);
MAP_ch_m = zeros(num_iter,num_scen); MAP_ch_f = zeros(num_iter,num_scen);
MAP_pc_m = zeros(num_iter,num_scen); MAP_pc_f = zeros(num_iter,num_scen);
for i = 1:num_scen
    MAP_ac_m(:,i) = (X_ss_m(42,:,i)                 )                        ;
    MAP_ac_f(:,i) = (X_ss_f(42,:,i)                 )                        ;
    
    MAP_ch_m(:,i) = (X_ss_m(42,:,i) - X_ss_m(42,1,i))                        ;
    MAP_ch_f(:,i) = (X_ss_f(42,:,i) - X_ss_f(42,1,i))                        ;
    
    MAP_pc_m(:,i) = (X_ss_m(42,:,i) - X_ss_m(42,1,i)) ./ X_ss_m(42,1,i) * 100;
    MAP_pc_f(:,i) = (X_ss_f(42,:,i) - X_ss_f(42,1,i)) ./ X_ss_f(42,1,i) * 100;
end

g = figure('DefaultAxesFontSize',14);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 12, 4]);
s_map(1) = subplot(1,3,1); 
s_map(2) = subplot(1,3,2); 
s_map(3) = subplot(1,3,3);

plot(s_map(1), xscale,MAP_ac_m(:,fixed_ss1), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
% xlim(s_map(1), [lower, upper]);
set(s_map(1), 'XTick', [0, 20, 40, 60, 80, 100]);
% ylim(s_map(1), [0.99,1.03])
xlabel(s_map(1), xlabel_name); ylabel(s_map(1), 'MAP (mmHg)');
hold(s_map(1), 'on')
plot(s_map(1), xscale,MAP_ac_f(:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold(s_map(1), 'off')
[~, hobj, ~, ~] = legend(s_map(1), {'Male','Female'}, 'FontSize',10,'Location','Southwest');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
title(s_map(1), 'A')

plot(s_map(2), xscale,MAP_ch_m(:,fixed_ss1), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
% xlim(s_map(2), [lower, upper]);
set(s_map(2), 'XTick', [0, 20, 40, 60, 80, 100]);
% ylim(s_map(2), [0.99,1.03])
xlabel(s_map(2), xlabel_name); ylabel(s_map(2), '\DeltaMAP (mmHg)');
hold(s_map(2), 'on')
plot(s_map(2), xscale,MAP_ch_f(:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold(s_map(2), 'off')
title(s_map(2), 'B')

plot(s_map(3), xscale,MAP_pc_m(:,fixed_ss1), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3,'MarkerSize',8);
% xlim(s_map(3), [lower, upper]);
set(s_map(3), 'XTick', [0, 20, 40, 60, 80, 100]);
% ylim(s_map(3), [0.99,1.03])
xlabel(s_map(3), xlabel_name); ylabel(s_map(3), '% \DeltaMAP');
hold(s_map(3), 'on')
plot(s_map(3), xscale,MAP_pc_f(:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold(s_map(3), 'off')
title(s_map(3), 'C')

% %% Save figures. -----------------------------------------------------------
%  
% save_data_name = sprintf('dose_response_%s%s.fig', scenario2{fixed_ss2},num2str(sample_num));
% save_data_name = strcat('Figures/', save_data_name);
% savefig([f1;f2;g], save_data_name)

end






























