function create_vp_dist

close all

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Begin user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters to perturb
% K_bar            - par 13; - 0%, +200%
% R_bv             - par 14; - 0%, +200%
% % C_gcf            - par 8 ; -20%
% R_aass           - par 4 ; - 0%, +100%
% N_rs             - par 21; - 0%, +100%
% N_als_eq         - par 18; - 0%, +100%
% N_rsna           - par 3 ; - 0%, +100%

% Range for parameters
pars_range_lower = [0  ;0  ;0  ;0  ;0  ;0  ]/100;
pars_range_upper = [200;200;100;100;100;100]/100;

pars_ind     = [13;14;4;21;18;3];
pars_hyp_num = length(pars_ind);
pars_names   = {'$K_{bar}$', '$R_{bv}$'      , '$R_{aa-ss}$', ...
                '$N_{rs}$' , '$N_{als}^{eq}$', '$N_{rsna}$' };

vars_ind     = [42;33;41;29;30;52;6;7;92];
vars_hyp_num = length(vars_ind);
vars_names   = {'$P_{ma}$'        , '$\Phi_{co}$'     , '$R_{tp}$'             , ...
                '$V_{ecf}$'       , '$V_{b}$'         , '$C_{sod}$'            , ...
                '$\Phi_{rb}$'     , '$\Phi_{gfilt}$'  , '$\Phi_{u}$'           };
vars_units   = {'$mmHg$'          , '$\frac{ml}{min}$', '$\frac{mmHg}{ml/min}$', ...
                '$ml$'            , '$ml$'            , '$\frac{\mu eq}{ml}$'  , ...
                '$\frac{ml}{min}$', '$\frac{ml}{min}$', '$\frac{ml}{min}$'     };

% Scenario
scenario = {'Normal'};

% Species
spe_ind = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

species = {'human', 'rat'   };
sex     = {'male' , 'female'};

for sex_ind = 1:2 % sex

%% Load bootstrap replicate parameters & variables created by create_par_bs_rep.m.

% Parameters
load_data_name_pars = sprintf('%s_%s_pars_scenario_Pri_Hyp_bs_rep1000.mat', ...
                              species{spe_ind},sex{sex_ind});
load(load_data_name_pars, 'pars_rep');
pars_hyp = pars_rep(pars_ind,:);

num_pars   = size(pars_rep,1);
num_sample = size(pars_rep,2);

% Variables
load_data_name_vars = sprintf('%s_%s_ss_data_scenario_Pri_Hyp_bs_rep1000.mat', ...
                              species{spe_ind},sex{sex_ind});
load(load_data_name_vars, 'SSdata_rep');
vars_hyp = SSdata_rep(vars_ind,:);

num_vars   = size(SSdata_rep,1);

%% Load baseline parameters and variables for comparison.

% Parameters
varargin_input = {scenario{1},true};
pars_bl = get_pars(species{spe_ind}, sex{sex_ind}, varargin_input{:});
pars_hyp_bl = pars_bl(pars_ind,:);

% Variables
load_data_name_vars_bl = sprintf('%s_%s_ss_data_scenario_%s.mat', ...
                         species{spe_ind},sex{sex_ind},scenario{1});
load(load_data_name_vars_bl, 'SSdata');
SSdata_bl = SSdata; clear SSdata;
vars_hyp_bl = SSdata_bl(vars_ind,:);

%% Plot parameter and variable distributions.

% % Parameter bin edges.
% edges_pars = [];

% Plot parameters
f(sex_ind,1) = figure('DefaultAxesFontSize',14);
s1   = gobjects(pars_hyp_num);
for i = 1:pars_hyp_num
    s1(i) = subplot(3,2,i);
    histogram(s1(i),pars_hyp(i,:),10)
    xlabel(s1(i), pars_names(i), 'Interpreter','latex', 'FontSize',16)
end
hist_title = sprintf('%s',sex{sex_ind});
sgtitle(hist_title, 'FontSize',16)

% % Variable bin edges.
% edges_vars = [];

% Variable number of bins.
% P_ma_edge      = linspace(min(vars_hyp(1,:)),max(vars_hyp(1,:)),7);
P_ma_edge      = (floor(min(vars_hyp(1,:))) : 2 : floor(min(vars_hyp(1,:)))+12);
% Phi_co_edge    = linspace(min(vars_hyp(2,:)),max(vars_hyp(2,:)),4);
% Phi_co_edge    = [0.95*vars_hyp_bl(2), 1.0*vars_hyp_bl(2), 1.05*vars_hyp_bl(2)];
% Phi_co_edge    = [0.91*vars_hyp_bl(2), 0.97*vars_hyp_bl(2), 1.03*vars_hyp_bl(2), 1.09*vars_hyp_bl(2)];
Phi_co_edge    = [0.90*vars_hyp_bl(2), 0.95*vars_hyp_bl(2), ...
                  1.00*vars_hyp_bl(2), 1.05*vars_hyp_bl(2), 1.10*vars_hyp_bl(2)];
% C_sod_edge     = [0.95*vars_hyp_bl(6), 1.0*vars_hyp_bl(6), 1.05*vars_hyp_bl(6)];
C_sod_edge     = [0.90*vars_hyp_bl(6), 0.95*vars_hyp_bl(6), ...
                  1.00*vars_hyp_bl(6), 1.05*vars_hyp_bl(6), 1.10*vars_hyp_bl(6)];
% Phi_rb_edge    = [0.95*vars_hyp_bl(7), 1.0*vars_hyp_bl(7), 1.05*vars_hyp_bl(7)];
Phi_rb_edge    = [0.90*vars_hyp_bl(7), 0.95*vars_hyp_bl(7), ...
                  1.00*vars_hyp_bl(7), 1.05*vars_hyp_bl(7), 1.10*vars_hyp_bl(7)];
% Phi_gfilt_edge = [0.95*vars_hyp_bl(8), 1.0*vars_hyp_bl(8), 1.05*vars_hyp_bl(8)];
Phi_gfilt_edge = [0.90*vars_hyp_bl(8), 0.95*vars_hyp_bl(8), ...
                  1.00*vars_hyp_bl(8), 1.05*vars_hyp_bl(8), 1.10*vars_hyp_bl(8)];
% Phi_u_edge     = [0.95*vars_hyp_bl(9), 1.0*vars_hyp_bl(9), 1.05*vars_hyp_bl(9)];
Phi_u_edge     = [0.90*vars_hyp_bl(9), 0.95*vars_hyp_bl(9), ...
                  1.00*vars_hyp_bl(9), 1.05*vars_hyp_bl(9), 1.10*vars_hyp_bl(9)];
bins_vars = {P_ma_edge;Phi_co_edge;6;6;6;C_sod_edge;Phi_rb_edge;Phi_gfilt_edge;Phi_u_edge};
% bins_vars = {6;3;10;10;10;3;3;3;3};

% Variable bin height for 

% Plot variables
f(sex_ind,2) = figure('DefaultAxesFontSize',14);
s2   = gobjects(vars_hyp_num);
for i = 1:vars_hyp_num
    s2(i) = subplot(3,3,i);
    h = histogram(s2(i),vars_hyp(i,:),bins_vars{i});
    
%     hold(s2(i),'on')
%     plot(s2(i),[vars_hyp_bl(i),vars_hyp_bl(i)],[0,max(h.Values)])
%     hold(s2(i),'off')
    
    xlabel_name = strcat(vars_names(i), ' (', num2str(vars_hyp_bl(i),3), vars_units(i), ')');
    xlabel(s2(i), xlabel_name, 'Interpreter','latex', 'FontSize',16)
end
hist_title = sprintf('%s',sex{sex_ind});
sgtitle(hist_title, 'FontSize',16)

end % sex

%% Save figures.

save_data_name = sprintf('par_var_dist1000.fig');
save_data_name = strcat('Figures/', save_data_name);
savefig([f(1,:),f(2,:)], save_data_name)

end % create_vp_dist


























