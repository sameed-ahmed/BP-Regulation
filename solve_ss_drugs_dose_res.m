% This simulates the blood pressure regulation model bp_reg.m for drug administration.
% 
% Steady state data is calculated by solve_ss_scenario.m.

function solve_ss_drugs_dose_res

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
% ACEi   - Angiotensin converting enzyme inhibitor % 95
% ARB1   - Angiotensin receptor 1 blocker % 94
% CCB    - Calcium channel blocker % 84
% TZD    - Thiazide diuretic % 0.5 1?
% ARB2   - Angiotensin receptor 2 blocker %
% DRI    - Direct renin inhibitor %
% MRB    - Aldosterone blocker (MR?) %
% RSS    - Renin secretion stimulator (thiazide?) % % NOT COMPLETE
% AngII  - Ang II infusion fmol/(ml min)
scenario2 = {'Normal', 'ACEi', 'ARB1', 'CCB', 'TZD', ...
             'ARB2'  , 'DRI' , 'MRB' , 'RSS', 'AngII'};
fixed_ss2 = [5];

% Species
spe_ind = 2;

% Bootstrap replicate sample number
% sample_num = random('Discrete Uniform',1000)
% sample_num = 208
% sample_num = 655
% fixed_sample = 655;
% num_sample = 1000;
% num_sample = 5;
% fixed_sample = 1;
num_sample = 1000;

% Number of intervals for dose
% num_iter = 21;
% drug_dose = linspace(0,1.00,num_iter);
num_dose = 100;
% num_dose = 5;
drug_dose = linspace(0,0.99,num_dose);
% drug_dose_vaso = 0.2 * drug_dose;
drug_dose_vaso = zeros(1,num_dose);
% a = 3; b = 1;
a = 11/9; b = 1/9;
% drug_dose_rsec = 2*drug_dose
drug_dose_rsec = a * drug_dose ./ (b + drug_dose);
% drug_dose_rsec = zeros(1,num_iter)
% drug_dose_rsec = ones(1,num_iter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

species = {'human', 'rat'   };
sex     = {'male' , 'female'};

% Number of variables
num_vars = 93;

% Initialize variables.
% X = (variable, sample, iteration, sex, scenario)
X_ss = zeros(num_vars,num_sample,num_dose,2,num_scen);
X_bl = zeros(num_vars,num_sample,2,num_scen);

parpool
for sce_ind = fixed_ss1:fixed_ss1 % scenario
for sex_ind = 1:2        % sex

%% Load bootstrap replicate parameters & variables created by create_par_bs_rep.m.

% Parameters
load_data_name_pars = sprintf('%s_%s_pars_scenario_Pri_Hyp_bs_rep1000NEWNEW.mat', ...
                              species{spe_ind},sex{sex_ind});
load(load_data_name_pars, 'pars_rep');
num_pars   = size(pars_rep,1);

% Variables
load_data_name_vars = sprintf('%s_%s_ss_data_scenario_Pri_Hyp_bs_rep1000NEWNEW.mat', ...
                              species{spe_ind},sex{sex_ind});
load(load_data_name_vars, 'SSdata_rep');
num_vars = size(SSdata_rep,1);
% Store baseline value to compute relative change.
X_bl(:,:,sex_ind,sce_ind) = SSdata_rep(:,1:num_sample);

parfor sam_iter = 1:num_sample % samples
% for sam_iter = 1:num_sample % samples
% for sam_iter = 1:10

for dose_iter = 1:num_dose % range

%% Drugs

varargin_input = {scenario1{sce_ind},true};

for i = 1:length(fixed_ss2)
    if     strcmp(scenario2{fixed_ss2(i)}, 'ACEi' )
            varargin_input = [varargin_input, 'ACEi' ,drug_dose(dose_iter)]; % 
    elseif strcmp(scenario2{fixed_ss2(i)}, 'ARB1' )
            varargin_input = [varargin_input, 'ARB1' ,drug_dose(dose_iter)]; % 
    elseif strcmp(scenario2{fixed_ss2(i)}, 'CCB'  )
            varargin_input = [varargin_input, 'CCB'  ,[drug_dose(dose_iter),2/3]]; % 
    elseif strcmp(scenario2{fixed_ss2(i)}, 'TZD'  )
        if     strcmp(sex{sex_ind}, 'male')
            varargin_input = [varargin_input, 'TZD'  ,[drug_dose(dose_iter)/1.0,drug_dose_vaso(dose_iter),drug_dose_rsec(dose_iter)]]; % 
        elseif strcmp(sex{sex_ind}, 'female')
            varargin_input = [varargin_input, 'TZD'  ,[drug_dose(dose_iter)/1.0,drug_dose_vaso(dose_iter),drug_dose_rsec(dose_iter)]]; % 
        end
    elseif strcmp(scenario2{fixed_ss2(i)}, 'ARB2' )
            varargin_input = [varargin_input, 'ARB2' ,drug_dose(dose_iter)]; % 
    elseif strcmp(scenario2{fixed_ss2(i)}, 'DRI'  )
            varargin_input = [varargin_input, 'DRI'  ,drug_dose(dose_iter)]; % 
    elseif strcmp(scenario2{fixed_ss2(i)}, 'MRB'  )
            varargin_input = [varargin_input, 'MRB'  ,drug_dose(dose_iter)]; % 
    elseif strcmp(scenario2{fixed_ss2(i)}, 'RSS'  )
            varargin_input = [varargin_input, 'RSS'  ,drug_dose(dose_iter)]; % 

    elseif strcmp(scenario2{fixed_ss2(i)}, 'AngII')
        if     strcmp(sex{sex_ind}, 'male')
            varargin_input = [varargin_input, 'AngII',910]; % Sullivan 2010
        elseif strcmp(sex{sex_ind}, 'female')
            varargin_input = [varargin_input, 'AngII',505]; % Sullivan 2010
        end
    end
end

%% Solve system steady state

% Initialization

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
% X = (variable, sample, iteration, sex, scenario)
X_ss(:,sam_iter,dose_iter,sex_ind,sce_ind) = SSdata;

end % range

% Sanity check to see script's progress. Also a check for where to
% troubleshoot in case the solver does not converge.
fprintf('%s %s sample = %s out of %s \n', ...
        scenario1{sce_ind},sex{sex_ind},num2str(sam_iter),num2str(num_sample))

end % samples
end % sex
end % scenario
delete(gcp)

%% Post processing

% Retrieve male and female data. Find relative change, mean, standard
% deviation, etc. in order to save data for loading in post 
% processing script.

% Retrieve male and female.
% Delete other scenarios for now.
X_ss_m = reshape(X_ss(:,:,:,1,1), [num_vars,num_sample,num_dose]); 
X_ss_f = reshape(X_ss(:,:,:,2,1), [num_vars,num_sample,num_dose]); 
X_bl_m = reshape(X_bl(:,:,  1,1), [num_vars,num_sample         ]); 
X_bl_f = reshape(X_bl(:,:,  2,1), [num_vars,num_sample         ]); 
% X_ss_m = reshape(X_ss(:,1:10,:,1,1), [num_vars,10,num_dose]); 
% X_ss_f = reshape(X_ss(:,1:10,:,2,1), [num_vars,10,num_dose]); 
% X_bl_m = reshape(X_bl(:,1:10,  1,1), [num_vars,10         ]); 
% X_bl_f = reshape(X_bl(:,1:10,  2,1), [num_vars,10         ]); 

% Compute relative change.
X_rel_m = zeros(num_vars,num_sample,num_dose);
X_rel_f = zeros(num_vars,num_sample,num_dose);
% X_rel_m = zeros(num_vars,10,num_dose);
% X_rel_f = zeros(num_vars,10,num_dose);
for i = 1:num_dose
    X_rel_m(:,:,i) = (X_ss_m(:,:,i) - X_bl_m) ./ X_bl_m * 100;
    X_rel_f(:,:,i) = (X_ss_f(:,:,i) - X_bl_f) ./ X_bl_f * 100;
end

% Compute mean and standard deviation.
X_bl_mean_m  = mean(X_bl_m ,  2); 
X_bl_mean_f  = mean(X_bl_f ,  2);
X_bl_std_m   = std (X_bl_m ,0,2); 
X_bl_std_f   = std (X_bl_f ,0,2);
% ---
X_ss_mean_m  = reshape(mean(X_ss_m ,  2), [num_vars,num_dose]); 
X_ss_mean_f  = reshape(mean(X_ss_f ,  2), [num_vars,num_dose]); 
X_ss_std_m   = reshape(std (X_ss_m ,0,2), [num_vars,num_dose]); 
X_ss_std_f   = reshape(std (X_ss_f ,0,2), [num_vars,num_dose]); 
% ---
X_rel_mean_m = reshape(mean(X_rel_m,  2), [num_vars,num_dose]); 
X_rel_mean_f = reshape(mean(X_rel_f,  2), [num_vars,num_dose]); 
X_rel_std_m  = reshape(std (X_rel_m,0,2), [num_vars,num_dose]); 
X_rel_std_f  = reshape(std (X_rel_f,0,2), [num_vars,num_dose]); 

% Compute 95% confidence interval
X_bl_lower_m  = X_bl_mean_m  - 2*X_bl_std_m ; X_bl_upper_m  = X_bl_mean_m  + 2*X_bl_std_m ;
X_bl_lower_f  = X_bl_mean_f  - 2*X_bl_std_f ; X_bl_upper_f  = X_bl_mean_f  + 2*X_bl_std_f ;
% ---
X_ss_lower_m  = X_ss_mean_m  - 2*X_ss_std_m ; X_ss_upper_m  = X_ss_mean_m  + 2*X_ss_std_m ;
X_ss_lower_f  = X_ss_mean_f  - 2*X_ss_std_f ; X_ss_upper_f  = X_ss_mean_f  + 2*X_ss_std_f ;
% ---
X_rel_lower_m = X_rel_mean_m - 2*X_rel_std_m; X_rel_upper_m = X_rel_mean_m + 2*X_rel_std_m;
X_rel_lower_f = X_rel_mean_f - 2*X_rel_std_f; X_rel_upper_f = X_rel_mean_f + 2*X_rel_std_f;

%% Plot

% x-axis
xscale = drug_dose * 100;
xlower = drug_dose(1) * 100; xupper = drug_dose(end) * 100; 

% y-axis limits
ylower = zeros(num_vars,1); yupper = ylower; 
for i = 1:length(ylower)
    ylower(i) = 0.95*min( min(X_ss_mean_m(i,:)), min(X_ss_mean_f(i,:)) );
    yupper(i) = 1.05*max( max(X_ss_mean_m(i,:)), max(X_ss_mean_f(i,:)) );
    if abs(yupper(i)) < eps*100
        ylower(i) = -10^(-5); yupper(i) = 10^(-5);
    end
end

%% Plot all vars vs dose. 

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
        
        plot(s1(i,j), xscale,X_ss_mean_m ((i-1)*15 + j,:),'b-' , ...
                      xscale,X_ss_mean_f ((i-1)*15 + j,:),'r-' );
        hold(s1(i,j), 'on')
        plot(s1(i,j), xscale,X_ss_lower_m((i-1)*15 + j,:),'b--', ...
                      xscale,X_ss_lower_f((i-1)*15 + j,:),'r--');
        plot(s1(i,j), xscale,X_ss_upper_m((i-1)*15 + j,:),'b--', ...
                      xscale,X_ss_upper_f((i-1)*15 + j,:),'r--');
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

f2 = figure('pos',[000 000 600 600], 'DefaultAxesFontSize',12);
s2 = gobjects(1,num_vars_sub);
% Loop through each subplot within a set of subplots.
for j = 1:num_vars_sub
    s2(j) = subplot(4,3,j);
%     if     mod(j,3) == 1
%         hshift = -0.05;
%     elseif mod(j,3) == 0
%         hshift = 0.05;
%     else
%         hshift = 0;
%     end
%     s3(j).Position = s3(j).Position + [hshift 0 0.01 0.01];

    plot(s2(j), xscale,X_ss_mean_m (var_ind(j),:,fixed_ss1), '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',2.5);
    hold(s2(j), 'on')
    plot(s2(j), xscale,X_ss_mean_f (var_ind(j),:,fixed_ss1), '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',2.5);
    plot(s2(j), xscale,X_ss_lower_m(var_ind(j),:,fixed_ss1), '--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',2.5);
    plot(s2(j), xscale,X_ss_lower_f(var_ind(j),:,fixed_ss1), '--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',2.5);
    plot(s2(j), xscale,X_ss_upper_m(var_ind(j),:,fixed_ss1), '--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',2.5);
    plot(s2(j), xscale,X_ss_upper_f(var_ind(j),:,fixed_ss1), '--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',2.5);
    hold(s2(j), 'off')

    xlabel(xlabel_name)
%     ylabel(var_names(var_ind(j)), 'Interpreter','latex')
    title(var_names(var_ind(j)), 'Interpreter','latex', 'FontSize',16)
        
%     xlabel(s2(j), var_names(var_ind(j)), 'Interpreter','latex', 'FontSize',16)
end

%% Plot variables that explain mechanims. 

% Efferent to total arteriolar resistance ---------------------------------
R_ss_m = reshape(X_ss_m(74,:,:) ./ X_ss_m(4,:,:), [num_sample,num_dose]);
R_ss_f = reshape(X_ss_f(74,:,:) ./ X_ss_f(4,:,:), [num_sample,num_dose]);
R_bl_m = reshape(X_bl_m(74,:  ) ./ X_bl_m(4,:  ), [num_sample,1       ]);
R_bl_f = reshape(X_bl_f(74,:  ) ./ X_bl_f(4,:  ), [num_sample,1       ]);
% R_ss_m = reshape(X_ss_m(74,:,:) ./ X_ss_m(4,:,:), [10,num_dose]);
% R_ss_f = reshape(X_ss_f(74,:,:) ./ X_ss_f(4,:,:), [10,num_dose]);
% R_bl_m = reshape(X_bl_m(74,:  ) ./ X_bl_m(4,:  ), [10,1       ]);
% R_bl_f = reshape(X_bl_f(74,:  ) ./ X_bl_f(4,:  ), [10,1       ]);
% ---
R_rel_m = (R_ss_m - R_bl_m) ./ R_bl_m * 100;
R_rel_f = (R_ss_f - R_bl_f) ./ R_bl_f * 100;
% ---
R_rel_mean_m  = reshape(mean(R_rel_m,  1), [1,num_dose]); 
R_rel_mean_f  = reshape(mean(R_rel_f,  1), [1,num_dose]); 
R_rel_std_m   = reshape(std (R_rel_m,0,1), [1,num_dose]); 
R_rel_std_f   = reshape(std (R_rel_f,0,1), [1,num_dose]); 
% Compute confidence interval
R_rel_lower_m = R_rel_mean_m - R_rel_std_m; R_rel_upper_m = R_rel_mean_m + R_rel_std_m;
R_rel_lower_f = R_rel_mean_f - R_rel_std_f; R_rel_upper_f = R_rel_mean_f + R_rel_std_f;

% Total fractional sodium reabsorption ------------------------------------
FRNA_ss_m = reshape((X_ss_m(11,:,:) - X_ss_m(27,:,:)) ./ X_ss_m(11,:,:), [num_sample,num_dose]) * 100;
FRNA_ss_f = reshape((X_ss_f(11,:,:) - X_ss_f(27,:,:)) ./ X_ss_f(11,:,:), [num_sample,num_dose]) * 100;
FRNA_bl_m = reshape((X_bl_m(11,:  ) - X_bl_m(27,:  )) ./ X_bl_m(11,:  ), [num_sample,1       ]) * 100;
FRNA_bl_f = reshape((X_bl_f(11,:  ) - X_bl_f(27,:  )) ./ X_bl_f(11,:  ), [num_sample,1       ]) * 100;
% FRNA_ss_m = reshape((X_ss_m(11,:,:) - X_ss_m(27,:,:)) ./ X_ss_m(11,:,:), [10,num_dose]) * 100;
% FRNA_ss_f = reshape((X_ss_f(11,:,:) - X_ss_f(27,:,:)) ./ X_ss_f(11,:,:), [10,num_dose]) * 100;
% FRNA_bl_m = reshape((X_bl_m(11,:  ) - X_bl_m(27,:  )) ./ X_bl_m(11,:  ), [10,1       ]) * 100;
% FRNA_bl_f = reshape((X_bl_f(11,:  ) - X_bl_f(27,:  )) ./ X_bl_f(11,:  ), [10,1       ]) * 100;
% ---
FRNA_rel_m = (FRNA_ss_m - FRNA_bl_m) ./ FRNA_bl_m * 100;
FRNA_rel_f = (FRNA_ss_f - FRNA_bl_f) ./ FRNA_bl_f * 100;
% ---
FRNA_rel_mean_m  = reshape(mean(FRNA_rel_m,  1), [1,num_dose]); 
FRNA_rel_mean_f  = reshape(mean(FRNA_rel_f,  1), [1,num_dose]); 
FRNA_rel_std_m   = reshape(std (FRNA_rel_m,0,1), [1,num_dose]); 
FRNA_rel_std_f   = reshape(std (FRNA_rel_f,0,1), [1,num_dose]); 
% Compute confidence interval
FRNA_rel_lower_m = FRNA_rel_mean_m - FRNA_rel_std_m; FRNA_rel_upper_m = FRNA_rel_mean_m + FRNA_rel_std_m;
FRNA_rel_lower_f = FRNA_rel_mean_f - FRNA_rel_std_f; FRNA_rel_upper_f = FRNA_rel_mean_f + FRNA_rel_std_f;

% Total fractional water reabsorption -------------------------------------
FRW_ss_m = reshape((X_ss_m(7,:,:) - X_ss_m(92,:,:)) ./ X_ss_m(7,:,:), [num_sample,num_dose]) * 100;
FRW_ss_f = reshape((X_ss_f(7,:,:) - X_ss_f(92,:,:)) ./ X_ss_f(7,:,:), [num_sample,num_dose]) * 100;
FRW_bl_m = reshape((X_bl_m(7,:  ) - X_bl_m(92,:  )) ./ X_bl_m(7,:  ), [num_sample,1       ]) * 100;
FRW_bl_f = reshape((X_bl_f(7,:  ) - X_bl_f(92,:  )) ./ X_bl_f(7,:  ), [num_sample,1       ]) * 100;
% FRW_ss_m = reshape((X_ss_m(7,:,:) - X_ss_m(92,:,:)) ./ X_ss_m(7,:,:), [10,num_dose]) * 100;
% FRW_ss_f = reshape((X_ss_f(7,:,:) - X_ss_f(92,:,:)) ./ X_ss_f(7,:,:), [10,num_dose]) * 100;
% FRW_bl_m = reshape((X_bl_m(7,:  ) - X_bl_m(92,:  )) ./ X_bl_m(7,:  ), [10,1       ]) * 100;
% FRW_bl_f = reshape((X_bl_f(7,:  ) - X_bl_f(92,:  )) ./ X_bl_f(7,:  ), [10,1       ]) * 100;
% ---
FRW_rel_m = (FRW_ss_m - FRW_bl_m) ./ FRW_bl_m * 100;
FRW_rel_f = (FRW_ss_f - FRW_bl_f) ./ FRW_bl_f * 100;
% ---
FRW_rel_mean_m  = reshape(mean(FRW_rel_m,  1), [1,num_dose]); 
FRW_rel_mean_f  = reshape(mean(FRW_rel_f,  1), [1,num_dose]); 
FRW_rel_std_m   = reshape(std (FRW_rel_m,0,1), [1,num_dose]); 
FRW_rel_std_f   = reshape(std (FRW_rel_f,0,1), [1,num_dose]); 
% Compute confidence interval
FRW_rel_lower_m = FRW_rel_mean_m - FRW_rel_std_m; FRW_rel_upper_m = FRW_rel_mean_m + FRW_rel_std_m;
FRW_rel_lower_f = FRW_rel_mean_f - FRW_rel_std_f; FRW_rel_upper_f = FRW_rel_mean_f + FRW_rel_std_f;

g1 = figure('DefaultAxesFontSize',14);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 4]);
t1 = tiledlayout(1,3,'TileSpacing','Normal','Padding','Compact');

nexttile
plot(xscale,R_rel_mean_m , '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
% xlim([xlower, xupper]);
xticks(0:20:100);
xlabel(xlabel_name); ylabel('EAR/RVR (relative)');
hold on
plot(xscale,R_rel_mean_f , '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(xscale,R_rel_lower_m, '--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(xscale,R_rel_lower_f, '--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(xscale,R_rel_upper_m, '--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(xscale,R_rel_upper_f, '--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold off
legend('Male', 'Female')
% [~, hobj, ~, ~] = legend({'Male','Female'}, 'FontSize',7,'Location','Southeast');
% hl = findobj(hobj,'type','line');
% set(hl,'LineWidth',1.5);
title('A')

nexttile
plot(xscale,FRNA_rel_mean_m , '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
% xlim([xlower, xupper]);
xticks(0:20:100);
xlabel(xlabel_name); ylabel('FR_{Na^+} (relative)');
hold on
plot(xscale,FRNA_rel_mean_f , '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(xscale,FRNA_rel_lower_m, '--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(xscale,FRNA_rel_lower_f, '--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(xscale,FRNA_rel_upper_m, '--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(xscale,FRNA_rel_upper_f, '--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold off
title('B')

nexttile
plot(xscale,FRW_rel_mean_m , '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
% xlim([xlower, xupper]);
xticks(0:20:100);
xlabel(xlabel_name); ylabel('FR_{W} (relative)');
hold on
plot(xscale,FRW_rel_mean_f , '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(xscale,FRW_rel_lower_m, '--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(xscale,FRW_rel_lower_f, '--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(xscale,FRW_rel_upper_m, '--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(xscale,FRW_rel_upper_f, '--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold off
title('C')

%% Plot Mean Arterial Pressure distribution. 

% Actual, change, and % change in MAP. ------------------------------------
MAP_ac_m = reshape((X_ss_m(42,:,:)               )                      , ...
    [num_sample,num_dose]);
%     [10,num_dose]);
MAP_ac_f = reshape((X_ss_f(42,:,:)               )                      , ...
    [num_sample,num_dose]);
%     [10,num_dose]);
% ---
MAP_ch_m = reshape((X_ss_m(42,:,:) - X_bl_m(42,:))                      , ...
    [num_sample,num_dose]);
%     [10,num_dose]);
MAP_ch_f = reshape((X_ss_f(42,:,:) - X_bl_f(42,:))                      , ...
    [num_sample,num_dose]);
%     [10,num_dose]);
% ---
MAP_pc_m = reshape((X_ss_m(42,:,:) - X_bl_m(42,:)) ./ X_bl_m(42,:) * 100, ...
    [num_sample,num_dose]);
%     [10,num_dose]);
MAP_pc_f = reshape((X_ss_f(42,:,:) - X_bl_f(42,:)) ./ X_bl_f(42,:) * 100, ...
    [num_sample,num_dose]);
%     [10,num_dose]);
% ---
MAP_ac_mean_m  = reshape(mean(MAP_ac_m,  1), [1,num_dose]); 
MAP_ac_mean_f  = reshape(mean(MAP_ac_f,  1), [1,num_dose]); 
MAP_ac_std_m   = reshape(std (MAP_ac_m,0,1), [1,num_dose]); 
MAP_ac_std_f   = reshape(std (MAP_ac_f,0,1), [1,num_dose]); 
% ---
MAP_ch_mean_m  = reshape(mean(MAP_ch_m,  1), [1,num_dose]); 
MAP_ch_mean_f  = reshape(mean(MAP_ch_f,  1), [1,num_dose]); 
MAP_ch_std_m   = reshape(std (MAP_ch_m,0,1), [1,num_dose]); 
MAP_ch_std_f   = reshape(std (MAP_ch_f,0,1), [1,num_dose]); 
% ---
MAP_pc_mean_m  = reshape(mean(MAP_pc_m,  1), [1,num_dose]); 
MAP_pc_mean_f  = reshape(mean(MAP_pc_f,  1), [1,num_dose]); 
MAP_pc_std_m   = reshape(std (MAP_pc_m,0,1), [1,num_dose]); 
MAP_pc_std_f   = reshape(std (MAP_pc_f,0,1), [1,num_dose]); 
% Compute confidence interval
MAP_ac_lower_m = MAP_ac_mean_m - MAP_ac_std_m; MAP_ac_upper_m = MAP_ac_mean_m + MAP_ac_std_m;
MAP_ac_lower_f = MAP_ac_mean_f - MAP_ac_std_f; MAP_ac_upper_f = MAP_ac_mean_f + MAP_ac_std_f;
% ---
MAP_ch_lower_m = MAP_ch_mean_m - MAP_ch_std_m; MAP_ch_upper_m = MAP_ch_mean_m + MAP_ch_std_m;
MAP_ch_lower_f = MAP_ch_mean_f - MAP_ch_std_f; MAP_ch_upper_f = MAP_ch_mean_f + MAP_ch_std_f;
% ---
MAP_pc_lower_m = MAP_pc_mean_m - MAP_pc_std_m; MAP_pc_upper_m = MAP_pc_mean_m + MAP_pc_std_m;
MAP_pc_lower_f = MAP_pc_mean_f - MAP_pc_std_f; MAP_pc_upper_f = MAP_pc_mean_f + MAP_pc_std_f;

g2 = figure('DefaultAxesFontSize',14);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 4]);
t2 = tiledlayout(1,3,'TileSpacing','Normal','Padding','Compact');

nexttile
plot(xscale,MAP_ac_mean_m , '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
% xlim([xlower, xupper]);
xticks(0:20:100);
xlabel(xlabel_name); ylabel('MAP (mmHg)');
hold on
plot(xscale,MAP_ac_mean_f , '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(xscale,MAP_ac_lower_m, '--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(xscale,MAP_ac_lower_f, '--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(xscale,MAP_ac_upper_m, '--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(xscale,MAP_ac_upper_f, '--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold off
legend('Male', 'Female')
% [~, hobj, ~, ~] = legend({'Male','Female'}, 'FontSize',7,'Location','Southeast');
% hl = findobj(hobj,'type','line');
% set(hl,'LineWidth',1.5);
title('A')

nexttile
plot(xscale,MAP_ch_mean_m , '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
% xlim([xlower, xupper]);
xticks(0:20:100);
xlabel(xlabel_name); ylabel('\DeltaMAP (mmHg)');
hold on
plot(xscale,MAP_ch_mean_f , '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(xscale,MAP_ch_lower_m, '--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(xscale,MAP_ch_lower_f, '--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(xscale,MAP_ch_upper_m, '--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(xscale,MAP_ch_upper_f, '--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold off
title('B')

nexttile
plot(xscale,MAP_pc_mean_m , '-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
% xlim([xlower, xupper]);
xticks(0:20:100);
xlabel(xlabel_name); ylabel('% \DeltaMAP');
hold on
plot(xscale,MAP_pc_mean_f , '-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(xscale,MAP_pc_lower_m, '--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(xscale,MAP_pc_lower_f, '--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
plot(xscale,MAP_pc_upper_m, '--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',3, 'MarkerSize',8);
plot(xscale,MAP_pc_upper_f, '--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',3, 'MarkerSize',8);
hold off
title('C')

%% % Save figures and data. 
%  
% save_data_name = sprintf('dose_response_%sNEW.fig', ...
%                          scenario2{fixed_ss2});
% save_data_name = strcat('Figures/', save_data_name);
% savefig([f1;f2;g1;g2], save_data_name)
% 
% save_data_name = sprintf('%s_male_ss_data_scenario_Pri_Hyp_%sNEW.mat'  , ...
%                          species{spe_ind},scenario2{fixed_ss2});
% save_data_name = strcat('Data/', save_data_name);
% save(save_data_name, 'X_bl_m' , 'X_bl_mean_m' , 'X_bl_std_m' , ...
%                      'X_ss_m' , 'X_ss_mean_m' , 'X_ss_std_m' , ...
%                      'X_rel_m', 'X_rel_mean_m', 'X_rel_std_m')
% 
% save_data_name = sprintf('%s_female_ss_data_scenario_Pri_Hyp_%sNEW.mat', ...
%                          species{spe_ind},scenario2{fixed_ss2});
% save_data_name = strcat('Data/', save_data_name);
% save(save_data_name, 'X_bl_f' , 'X_bl_mean_f' , 'X_bl_std_f' , ...
%                      'X_ss_f' , 'X_ss_mean_f' , 'X_ss_std_f' , ...
%                      'X_rel_f', 'X_rel_mean_f', 'X_rel_std_f')

end






























