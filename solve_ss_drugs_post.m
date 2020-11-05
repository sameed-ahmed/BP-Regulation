% This script performs post processing visualization on the result of
% administering a drug on the model for the entire inhibition range of the
% drug and the entire virtual population.

% Input
% fixed_ss2      : index of which drug to administer
% fixed_drug_dose: inhibition level to fix for distribution plots
% Output
% -plots and saves all variables relative change mean and 95% confidence 
%  interval versus drug dose inhibition level
% -plots and saves  all variables relative change distribution for a fixed 
%  drug dose inhibition level

function solve_ss_drugs_post

close all

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))
mypath = pwd;
mypath = strcat(mypath, '/Data/Large Files');
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

% Drug scenarios
% Normal - Normal conditions
% ACEi   - Angiotensin converting enzyme inhibitor % 
% ARB1   - Angiotensin receptor 1 blocker % 
% CCB    - Calcium channel blocker % 
% TZD    - Thiazide diuretic % 
scenario2 = {'Normal', 'ACEi', 'ARB1', 'CCB', 'TZD'};
fixed_ss2 = [2];

% Species
spe_ind = 2;

% Fixed drug dose to plot
fixed_drug_dose = 93+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Drug dose range
num_dose = 100;
drug_dose = linspace(0,0.99,num_dose);

% Pathophysiologically perturbed parameters
pars_ind = [13;14;4;21;18;3;15;41];
pars_hyp_num = length(pars_ind);

species = {'human', 'rat'   };
sex     = {'male' , 'female'};

% Number of variables
num_vars = 93;

%% Load bootstrap replicate parameters & variables before and after drug dose.

% Parameters
load_data_name_pars = sprintf(...
      '%s_male_pars_scenario_Pri_Hyp_bs_rep1000.mat', species{spe_ind});
load(load_data_name_pars, 'pars_rep');
num_pars = size(pars_rep,1);
pars_hyp_m = pars_rep(pars_ind,:);
load_data_name_pars = sprintf(...
    '%s_female_pars_scenario_Pri_Hyp_bs_rep1000.mat', species{spe_ind});
load(load_data_name_pars, 'pars_rep');
num_pars = size(pars_rep,1);
pars_hyp_f = pars_rep(pars_ind,:);
% Bootstrap replicate sample number
num_samples = size(pars_rep,2);

% Baseline parameters for relative change
varargin_input = {'Normal',true};
pars_bl_m = get_pars(species{spe_ind}, 'male'  , varargin_input{:}); 
pars_bl_f = get_pars(species{spe_ind}, 'female', varargin_input{:});
pars_bl_m = pars_bl_m(pars_ind); pars_bl_f = pars_bl_f(pars_ind); 

% Variables after drug dose and hypertensiv baseline before drug dose
% X_m/f = (variable, sample, scenario)
load_data_name_vars = sprintf('%s_male_ss_data_scenario_Pri_Hyp_%s.mat'  , ...
                              species{spe_ind},scenario2{fixed_ss2});
load(load_data_name_vars, ...
     'X_ss_m' , 'X_ss_mean_m' , 'X_ss_std_m' , ...
     'X_bl_m' , 'X_bl_mean_m' , 'X_bl_std_m' , ...
     'X_rel_m', 'X_rel_mean_m', 'X_rel_std_m');
load_data_name_vars = sprintf('%s_female_ss_data_scenario_Pri_Hyp_%s.mat', ...
                              species{spe_ind},scenario2{fixed_ss2});
load(load_data_name_vars, ...
     'X_ss_f' , 'X_ss_mean_f' , 'X_ss_std_f' , ...
     'X_bl_f' , 'X_bl_mean_f' , 'X_bl_std_f' , ...
     'X_rel_f', 'X_rel_mean_f', 'X_rel_std_f');

% Compute 95% confidence interval
X_bl_lower_m  = X_bl_mean_m  - 2*X_bl_std_m ; X_bl_upper_m  = X_bl_mean_m  + 2*X_bl_std_m ;
X_ss_lower_m  = X_ss_mean_m  - 2*X_ss_std_m ; X_ss_upper_m  = X_ss_mean_m  + 2*X_ss_std_m ;
X_rel_lower_m = X_rel_mean_m - 2*X_rel_std_m; X_rel_upper_m = X_rel_mean_m + 2*X_rel_std_m;
% ---
X_bl_lower_f  = X_bl_mean_f  - 2*X_bl_std_f ; X_bl_upper_f  = X_bl_mean_f  + 2*X_bl_std_f ;
X_ss_lower_f  = X_ss_mean_f  - 2*X_ss_std_f ; X_ss_upper_f  = X_ss_mean_f  + 2*X_ss_std_f ;
X_rel_lower_f = X_rel_mean_f - 2*X_rel_std_f; X_rel_upper_f = X_rel_mean_f + 2*X_rel_std_f;

%% Plot all vars vs dose. 

% x-axis
xscale = drug_dose * 100;
xlower = drug_dose(1) * 100; xupper = drug_dose(end) * 100; 

% y-axis limits
ylower = zeros(num_vars,1); yupper = ylower; 
for i = 1:length(ylower)
    ylower(i) = 0.95*min( min(X_rel_mean_m(i,:)), min(X_rel_mean_f(i,:)) );
    yupper(i) = 1.05*max( max(X_rel_mean_m(i,:)), max(X_rel_mean_f(i,:)) );
    if abs(yupper(i)) < eps*100
        ylower(i) = -10^(-5); yupper(i) = 10^(-5);
    end
end

f_dose_res1 = gobjects(7,1);
s1 = gobjects(7,15);
% Loop through each set of subplots.
for i = 1:7
    f_dose_res1(i) = figure('pos',[750 500 650 450]);
    % This is to avoid the empty plots in the last subplot set.
    if i == 7
        last_plot = mod(num_vars, 15);
    else
        last_plot = 15;
    end
    % Loop through each subplot within a set of subplots.
    for j = 1:last_plot
        s1(i,j) = subplot(3,5,j);

        plot(s1(i,j), xscale,X_rel_mean_m ((i-1)*15 + j,:),'-', 'Color',[0.203, 0.592, 0.835]);
        hold(s1(i,j), 'on')
        plot(s1(i,j), xscale,X_rel_mean_f ((i-1)*15 + j,:),'-', 'Color',[0.835, 0.203, 0.576]);
        plot(s1(i,j), xscale,X_rel_lower_m((i-1)*15 + j,:),':', 'Color',[0.203, 0.592, 0.835]);
        plot(s1(i,j), xscale,X_rel_lower_f((i-1)*15 + j,:),':', 'Color',[0.835, 0.203, 0.576]);
        plot(s1(i,j), xscale,X_rel_upper_m((i-1)*15 + j,:),':', 'Color',[0.203, 0.592, 0.835]);
        plot(s1(i,j), xscale,X_rel_upper_f((i-1)*15 + j,:),':', 'Color',[0.835, 0.203, 0.576]);
        hold(s1(i,j), 'off')
        
%         xlim([xlower, xupper])
%         ylim([ylower((i-1)*15 + j), yupper((i-1)*15 + j)])

        xlabel_name = strcat(scenario2{fixed_ss2}, ' %');
        xlabel(xlabel_name, 'FontSize',12)
        title(var_names((i-1)*15 + j), 'Interpreter','latex', 'FontSize',15)
%         legend('Male', 'Female')
    end
end

%% Plot interesting variables. 

% Interesting variables to plot.
var_ind = [30;33;41;42;9;73;74;6;7;92;69;70]; num_vars_sub = length(var_ind);

f_dose_res2 = figure('pos',[000 000 600 600], 'DefaultAxesFontSize',12);
s2 = gobjects(1,num_vars_sub);
% Loop through each subplot within a set of subplots.
for j = 1:num_vars_sub
    s2(j) = subplot(4,3,j);
    
    plot(s2(j), xscale,X_rel_mean_m (var_ind(j),:),'-', 'Color',[0.203, 0.592, 0.835]);
    hold(s2(j), 'on')
    plot(s2(j), xscale,X_rel_mean_f (var_ind(j),:),'-', 'Color',[0.835, 0.203, 0.576]);
    plot(s2(j), xscale,X_rel_lower_m(var_ind(j),:),':', 'Color',[0.203, 0.592, 0.835]);
    plot(s2(j), xscale,X_rel_lower_f(var_ind(j),:),':', 'Color',[0.835, 0.203, 0.576]);
    plot(s2(j), xscale,X_rel_upper_m(var_ind(j),:),':', 'Color',[0.203, 0.592, 0.835]);
    plot(s2(j), xscale,X_rel_upper_f(var_ind(j),:),':', 'Color',[0.835, 0.203, 0.576]);
    hold(s2(j), 'off')


    xlabel_name = strcat(scenario2{fixed_ss2}, ' %');
    xlabel(xlabel_name, 'FontSize',12)
    title(var_names(var_ind(j)), 'Interpreter','latex', 'FontSize',15)
%         legend('Male', 'Female')
end

%% Plot all vars distribution for given dose.

f_dose_dist1 = gobjects(7,1);
ss1 = gobjects(7,15);
% Loop through each set of subplots.
for i = 1:7
    f_dose_dist1(i) = figure('pos',[750 500 650 450]);
    % This is to avoid the empty plots in the last subplot set.
    if i == 7
        last_plot = mod(num_vars, 15);
    else
        last_plot = 15;
    end
    % Loop through each subplot within a set of subplots.
    for j = 1:last_plot
        ss1(i,j) = subplot(3,5,j);
        
        h1 = histogram(ss1(i,j),X_rel_m((i-1)*15 + j,:,fixed_drug_dose),10);
        
        h1.Normalization = 'probability'; 
%         h1.BinWidth = 1.0; 
        h1.FaceColor = [0.203, 0.592, 0.835]; 

        xlabel_name = strcat(var_names((i-1)*15 + j));
        xlabel(ss1(i,j), xlabel_name, 'Interpreter','latex', 'FontSize',16)

    end
    hist_title = sprintf('Male, ACEi %s%%',num2str(fixed_drug_dose-1));
    sgtitle(hist_title, 'FontSize',14)
end

f_dose_dist2 = gobjects(7,1);
ss2 = gobjects(7,15);
% Loop through each set of subplots.
for i = 1:7
    f_dose_dist2(i) = figure('pos',[750 500 650 450]);
    % This is to avoid the empty plots in the last subplot set.
    if i == 7
        last_plot = mod(num_vars, 15);
    else
        last_plot = 15;
    end
    % Loop through each subplot within a set of subplots.
    for j = 1:last_plot
        ss2(i,j) = subplot(3,5,j);
        h1 = histogram(ss2(i,j),X_rel_f((i-1)*15 + j,:,fixed_drug_dose),10);
        
        h1.Normalization = 'probability'; 
%         h1.BinWidth = 1.0; 
        h1.FaceColor = [0.835, 0.203, 0.576];

        xlabel_name = strcat(var_names((i-1)*15 + j));
        xlabel(ss2(i,j), xlabel_name, 'Interpreter','latex', 'FontSize',16)

    end
    hist_title = sprintf('Female, ACEi %s%%',num2str(fixed_drug_dose-1));
    sgtitle(hist_title, 'FontSize',14)
end

%% Plot interesting variables. 

% Interesting variables to plot.
var_ind = [30;33;41;42;9;73;74;6;7;92;69;70]; num_vars_sub = length(var_ind);

f_dose_dist3 = figure('pos',[000 000 600 600], 'DefaultAxesFontSize',12);
ss3 = gobjects(1,num_vars_sub);
% Loop through each subplot within a set of subplots.
for j = 1:num_vars_sub
    ss3(j) = subplot(4,3,j);
    
        h1 = histogram(ss3(j),X_rel_m(var_ind(j),:,fixed_drug_dose),10);
        
        h1.Normalization = 'probability'; 
%         h1.BinWidth = 1.0; 
        h1.FaceColor = [0.203, 0.592, 0.835]; 

    xlabel(ss3(j), var_names(var_ind(j)), 'Interpreter','latex', 'FontSize',16)
end
hist_title = sprintf('Male, ACEi %s%%',num2str(fixed_drug_dose-1));
sgtitle(hist_title, 'FontSize',14)

f_dose_dist4 = figure('pos',[000 000 600 600], 'DefaultAxesFontSize',12);
ss4 = gobjects(1,num_vars_sub);
% Loop through each subplot within a set of subplots.
for j = 1:num_vars_sub
    ss4(j) = subplot(4,3,j);
    
        h1 = histogram(ss4(j),X_rel_f(var_ind(j),:,fixed_drug_dose),10);
        
        h1.Normalization = 'probability';
%         h1.BinWidth = 1.0;
        h1.FaceColor = [0.835, 0.203, 0.576];

    xlabel(ss4(j), var_names(var_ind(j)), 'Interpreter','latex', 'FontSize',16)
end
hist_title = sprintf('Female, ACEi %s%%',num2str(fixed_drug_dose-1));
sgtitle(hist_title, 'FontSize',14)

%% Save figures.

% Save all vars rel change vs dose. ---------------------------------------
save_data_name = sprintf('TESTdose_resp_rel_%s.fig', scenario2{fixed_ss2});
save_data_name = strcat('Figures/', save_data_name);
savefig([f_dose_res1;f_dose_res2], save_data_name)

% Save all vars distribution for fixed dose. ------------------------------
save_data_name = sprintf('TESTdose_dist_rel_%s.fig', scenario2{fixed_ss2});
save_data_name = strcat('Figures/', save_data_name);
savefig([f_dose_dist1;f_dose_dist2;f_dose_dist3;f_dose_dist4], save_data_name)

end






























