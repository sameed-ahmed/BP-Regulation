% This simulates the blood pressure regulation model bp_reg.m for drug administration.
% 
% Steady state data is calculated by solve_ss_scenario.m.

function solve_ss_drugs_post

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
fixed_ss2 = [4];

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

% Drug dose
drug_dose = 0.95;

% Mean arterial pressure threshold
MAP_th = -20;
% Glomerular filtration rate threshold
GFR_th = 25;

% Pathophysiologically perturbed parameters
pars_ind     = [13;14;4;21;18;3];
pars_hyp_num = length(pars_ind);
pars_names   = {'K_{bar}$'            , 'R_{bv}$'             , ...
                'R_{aa-ss}$'          , 'N_{rs}$'             , ...
                'N_{als}^{eq}$'       , 'N_{rsna}$'           };
pars_units   = {'$\frac{mmHg}{ml/min}$', '$\frac{mmHg}{ml/min}$', ...
                '$\frac{mmHg}{ml/min}$', '$-$'                  , ...
                '$-$'                  , '$-$'                  };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

species = {'human', 'rat'   };
sex     = {'male' , 'female'};

% Number of variables
num_vars = 93;

%% Load bootstrap replicate parameters & variables before and after drug dose.

% Parameters
load_data_name_pars = sprintf('%s_male_pars_scenario_Pri_Hyp_bs_rep1000.mat'  , ...
                              species{spe_ind});
load(load_data_name_pars, 'pars_rep');
num_pars = size(pars_rep,1);
pars_hyp_m = pars_rep(pars_ind,:);
load_data_name_pars = sprintf('%s_female_pars_scenario_Pri_Hyp_bs_rep1000.mat', ...
                              species{spe_ind});
load(load_data_name_pars, 'pars_rep');
num_pars = size(pars_rep,1);
pars_hyp_f = pars_rep(pars_ind,:);

% Baseline parameters for relative change
varargin_input = {'Normal',true};
pars_bl_m = get_pars(species{spe_ind}, 'male'  , varargin_input{:}); 
pars_bl_f = get_pars(species{spe_ind}, 'female', varargin_input{:});
pars_bl_m = pars_bl_m(pars_ind); pars_bl_f = pars_bl_f(pars_ind); 

% Variables after drug dose and hypertensiv baseline before drug dose
% X_m/f = (variable, sample, scenario)
load_data_name_vars = sprintf('%s_male_ss_data_scenario_Pri_Hyp_%s%s%%.mat'  , ...
                         species{spe_ind},scenario2{fixed_ss2},num2str(drug_dose*100));
load(load_data_name_vars, 'X_ss_m', 'X_bl_m');
load_data_name_vars = sprintf('%s_female_ss_data_scenario_Pri_Hyp_%s%s%%.mat', ...
                         species{spe_ind},scenario2{fixed_ss2},num2str(drug_dose*100));
load(load_data_name_vars, 'X_ss_f', 'X_bl_f');

%% Post processing

% Compute relative change in parameters.
pars_rel_m = (pars_hyp_m(:,:) - pars_bl_m(:)) ...
          ./ pars_bl_m(:) * 100;
pars_rel_f = (pars_hyp_f(:,:) - pars_bl_f(:)) ...
          ./ pars_bl_f(:) * 100;

% Compute relative change in variables.
% X_rel_m/f = (variable, sample)
X_rel_m = (X_ss_m(:,:,fixed_ss1) - X_bl_m(:,:,fixed_ss1)) ...
          ./ X_bl_m(:,:,fixed_ss1) * 100;
X_rel_f = (X_ss_f(:,:,fixed_ss1) - X_bl_f(:,:,fixed_ss1)) ...
          ./ X_bl_f(:,:,fixed_ss1) * 100;

% MAP success threshold indices
MAP_rel_m = X_rel_m(42,:); MAP_rel_f = X_rel_f(42,:);
MAP_success_ind_m = find(MAP_rel_m - MAP_th <= 0);
MAP_failure_ind_m = 1:num_samples; MAP_failure_ind_m(MAP_success_ind_m) = '';
num_success_m = length(MAP_success_ind_m); num_failure_m = length(MAP_failure_ind_m); 
MAP_success_ind_f = find(MAP_rel_f - MAP_th <= 0);
MAP_failure_ind_f = 1:num_samples; MAP_failure_ind_f(MAP_success_ind_f) = '';
num_success_f = length(MAP_success_ind_f); num_failure_f = length(MAP_failure_ind_f); 

% Success/failure relative change in parameters from normotensive to hypertensive.
pars_success_m = pars_rel_m(:,MAP_success_ind_m);
pars_success_f = pars_rel_f(:,MAP_success_ind_f);
pars_failure_m = pars_rel_m(:,MAP_failure_ind_m);
pars_failure_f = pars_rel_f(:,MAP_failure_ind_f);

% Success/failure hypertensive variables before drug dose.
X_success_m = X_bl_m(:,MAP_success_ind_m,fixed_ss1);
X_success_f = X_bl_f(:,MAP_success_ind_f,fixed_ss1);
X_failure_m = X_bl_m(:,MAP_failure_ind_m,fixed_ss1);
X_failure_f = X_bl_f(:,MAP_failure_ind_f,fixed_ss1);

% Compute statistical significance of difference between success and
% failure variables before treatment.
[h_m,p_m] = ttest2(X_success_m',X_failure_m');
[h_f,p_f] = ttest2(X_success_f',X_failure_f');
sig_m = [p_m',h_m']; sig_f = [p_f',h_f'];

%% % Plot hypertensive perturbed parameters success and failure
% 
% fp1 = figure('DefaultAxesFontSize',14);
% sp1 = gobjects(pars_hyp_num);
% for i = 1:pars_hyp_num
%     sp1(i) = subplot(3,2,i);
%     h1 = histogram(sp1(i),pars_success_m(i,:));
%     hold(sp1(i), 'on')
%     h2 = histogram(sp1(i),pars_failure_m(i,:));
%     hold(sp1(i), 'off')
%     
%     h1.Normalization = 'probability'; h2.Normalization = 'probability';  
%     h1.BinWidth = 10; h2.BinWidth = 10; 
%     h1.FaceColor = [0.070, 0.886, 0.333]; h2.FaceColor = [0.917, 0.121, 0.121];
% 
%     xlabel_name = strcat('$ \% \Delta ', {' '}, pars_names(i));
%     xlabel(sp1(i), xlabel_name, 'Interpreter','latex', 'FontSize',16)    
% end
% hist_title = sprintf('Male Parameters');
% sgtitle(hist_title, 'FontSize',16)
% 
% fp2 = figure('DefaultAxesFontSize',14);
% sp2 = gobjects(pars_hyp_num);
% for i = 1:pars_hyp_num
%     sp2(i) = subplot(3,2,i);
%     h1 = histogram(sp2(i),pars_success_f(i,:));
%     hold(sp2(i), 'on')
%     h2 = histogram(sp2(i),pars_failure_f(i,:));
%     hold(sp2(i), 'off')
%     
%     h1.Normalization = 'probability'; h2.Normalization = 'probability';  
%     h1.BinWidth = 10; h2.BinWidth = 10; 
%     h1.FaceColor = [0.070, 0.886, 0.333]; h2.FaceColor = [0.917, 0.121, 0.121];
% 
%     xlabel_name = strcat('$ \% \Delta ', {' '}, pars_names(i));
%     xlabel(sp2(i), xlabel_name, 'Interpreter','latex', 'FontSize',16)    
% end
% hist_title = sprintf('Female Parameters');
% sgtitle(hist_title, 'FontSize',16)

%% % Plot all variables success and failure
% 
% fv1 = gobjects(7,1);
% sv1 = gobjects(7,15);
% % Loop through each set of subplots.
% for i = 1:7
%     fv1(i) = figure('pos',[750 500 650 450]);
%     % This is to avoid the empty plots in the last subplot set.
%     if i == 7
%         last_plot = mod(num_vars, 15);
%     else
%         last_plot = 15;
%     end
%     % Loop through each subplot within a set of subplots.
%     for j = 1:last_plot
%         sv1(i,j) = subplot(3,5,j);
% %         s(i,j).Position = s(i,j).Position + [0 0 0.01 0];
%         
%         h1 = histogram(sv1(i,j),X_success_m((i-1)*15 + j,:),10);
%         hold on
%         h2 = histogram(sv1(i,j),X_failure_m((i-1)*15 + j,:),10);
%         hold off
%         
%         h1.Normalization = 'probability'; h2.Normalization = 'probability'; 
% %         h1.BinWidth = 1.0; h2.BinWidth = 1.0; 
%         h1.FaceColor = [0.070, 0.886, 0.333]; h2.FaceColor = [0.917, 0.121, 0.121];
% 
%         xlabel_name = strcat(var_names((i-1)*15 + j));
%         xlabel(sv1(i,j), xlabel_name, 'Interpreter','latex', 'FontSize',16)
% 
% %         legend('Male', 'Female')
%     end
%     hist_title = sprintf('Male Variables %s %s%%',scenario2{fixed_ss2},num2str(drug_dose*100));
%     sgtitle(hist_title, 'FontSize',14)
% end
% 
% fv2 = gobjects(7,1);
% sv2 = gobjects(7,15);
% % Loop through each set of subplots.
% for i = 1:7
%     fv2(i) = figure('pos',[750 500 650 450]);
%     % This is to avoid the empty plots in the last subplot set.
%     if i == 7
%         last_plot = mod(num_vars, 15);
%     else
%         last_plot = 15;
%     end
%     % Loop through each subplot within a set of subplots.
%     for j = 1:last_plot
%         sv2(i,j) = subplot(3,5,j);
% %         s(i,j).Position = s(i,j).Position + [0 0 0.01 0];
%         
%         h1 = histogram(sv2(i,j),X_success_f((i-1)*15 + j,:),10);
%         hold on
%         h2 = histogram(sv2(i,j),X_failure_f((i-1)*15 + j,:),10);
%         hold off
%         
%         h1.Normalization = 'probability'; h2.Normalization = 'probability'; 
% %         h1.BinWidth = 1.0; h2.BinWidth = 1.0; 
%         h1.FaceColor = [0.070, 0.886, 0.333]; h2.FaceColor = [0.917, 0.121, 0.121];
% 
%         xlabel_name = strcat(var_names((i-1)*15 + j));
%         xlabel(sv2(i,j), xlabel_name, 'Interpreter','latex', 'FontSize',16)
% 
% %         legend('Male', 'Female')
%     end
%     hist_title = sprintf('Female Variables %s %s%%',scenario2{fixed_ss2},num2str(drug_dose*100));
%     sgtitle(hist_title, 'FontSize',14)
% end

%% % Plot some interesting variables
% 
% R_bl_m = reshape(X_bl_m(74,:,:) ./ X_bl_m(4,:,:), [num_samples,num_scen]);
% R_bl_f = reshape(X_bl_f(74,:,:) ./ X_bl_f(4,:,:), [num_samples,num_scen]);
% size(R_bl_m(:,fixed_ss1))
% 
% FRNA_bl_m = reshape((X_bl_m(11,:,:) - X_bl_m(27,:,:)) ./ X_bl_m(11,:,:), [num_samples,num_scen]) * 100;
% FRNA_bl_f = reshape((X_bl_f(11,:,:) - X_bl_f(27,:,:)) ./ X_bl_f(11,:,:), [num_samples,num_scen]) * 100;
% 
% FRW_bl_m = reshape((X_bl_m( 7,:,:) - X_bl_m(92,:,:)) ./ X_bl_m( 7,:,:), [num_samples,num_scen]) * 100;
% FRW_bl_f = reshape((X_bl_f( 7,:,:) - X_bl_f(92,:,:)) ./ X_bl_f( 7,:,:), [num_samples,num_scen]) * 100;
% 
% g1 = figure('DefaultAxesFontSize',14);
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 12, 4]);
% s_mech(1) = subplot(1,3,1); 
% s_mech(2) = subplot(1,3,2); 
% s_mech(3) = subplot(1,3,3); 
% 
% h1 = histogram(s_mech(1),R_bl_m(:,fixed_ss1),10);
% hold(s_mech(1), 'on')
% h2 = histogram(s_mech(1),R_bl_f(:,fixed_ss1),10);
% hold(s_mech(1), 'off')
% h1.Normalization = 'probability'; h2.Normalization = 'probability';  
% h1.BinWidth = 0.01; h2.BinWidth = 0.01; 
% h1.FaceColor = [0.203, 0.592, 0.835]; h2.FaceColor = [0.835, 0.203, 0.576];
% legend(s_mech(1), [h1, h2],{'Male','Female'}, 'FontSize',10,'Location','Northeast');
% xlabel(s_mech(1), 'R_{EA}/R_R');
% title(s_mech(1), 'A')
% 
% h1 = histogram(s_mech(2),FRNA_bl_m(:,fixed_ss1),10);
% hold(s_mech(2), 'on')
% h2 = histogram(s_mech(2),FRNA_bl_f(:,fixed_ss1),10);
% hold(s_mech(2), 'off')
% h1.Normalization = 'probability'; h2.Normalization = 'probability';  
% h1.BinWidth = 0.01; h2.BinWidth = 0.01; 
% h1.FaceColor = [0.203, 0.592, 0.835]; h2.FaceColor = [0.835, 0.203, 0.576];
% xlabel(s_mech(2), 'FR_{Na^+}');
% title(s_mech(2), 'B')
% 
% h1 = histogram(s_mech(3),FRW_bl_m(:,fixed_ss1),10);
% hold(s_mech(3), 'on')
% h2 = histogram(s_mech(3),FRW_bl_f(:,fixed_ss1),10);
% hold(s_mech(3), 'off')
% h1.Normalization = 'probability'; h2.Normalization = 'probability';  
% h1.BinWidth = 0.02; h2.BinWidth = 0.02; 
% h1.FaceColor = [0.203, 0.592, 0.835]; h2.FaceColor = [0.835, 0.203, 0.576];
% xlabel(s_mech(3), 'FR_{U}');
% title(s_mech(3), 'C')

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

mean(FRNA_rel_m)
mean(FRNA_rel_f)

FRW_ss_m = reshape((X_ss_m( 7,:,:) - X_ss_m(92,:,:)) ./ X_ss_m( 7,:,:), [num_samples,num_scen]) * 100;
FRW_ss_f = reshape((X_ss_f( 7,:,:) - X_ss_f(92,:,:)) ./ X_ss_f( 7,:,:), [num_samples,num_scen]) * 100;
FRW_bl_m = reshape((X_bl_m( 7,:,:) - X_bl_m(92,:,:)) ./ X_bl_m( 7,:,:), [num_samples,num_scen]) * 100;
FRW_bl_f = reshape((X_bl_f( 7,:,:) - X_bl_f(92,:,:)) ./ X_bl_f( 7,:,:), [num_samples,num_scen]) * 100;

FRW_rel_m = (FRW_ss_m(:,fixed_ss1) - FRW_bl_m(:,fixed_ss1)) ...
            ./ FRW_bl_m(:,fixed_ss1) * 100;
FRW_rel_f = (FRW_ss_f(:,fixed_ss1) - FRW_bl_f(:,fixed_ss1)) ...
            ./ FRW_bl_f(:,fixed_ss1) * 100;

mean(FRW_rel_m)
mean(FRW_rel_f)

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
xlabel(s_mech(3), 'FR_{W}');
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

%% Save figures and data.

% save_data_name = sprintf('success_failure_distribution_%s%s%%.fig', ...
%                          scenario2{fixed_ss2},num2str(drug_dose*100));
% save_data_name = strcat('Figures/', save_data_name);
% savefig([fp1;fp2;fv1;fv2], save_data_name)
% 
% save_data_name = sprintf('%s_male_signif_scenario_Pri_Hyp_%s%s%%.mat'  , ...
%                          species{spe_ind},scenario2{fixed_ss2},num2str(drug_dose*100));
% save_data_name = strcat('Data/', save_data_name);
% save(save_data_name, 'sig_m', 'num_success_m', 'num_failure_m')
% 
% save_data_name = sprintf('%s_female_signif_scenario_Pri_Hyp_%s%s%%.mat', ...
%                          species{spe_ind},scenario2{fixed_ss2},num2str(drug_dose*100));
% save_data_name = strcat('Data/', save_data_name);
% save(save_data_name, 'sig_f', 'num_success_f', 'num_failure_f')

end






























