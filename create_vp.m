function create_vp

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
% Indices
pars_ind = [13;14;4;21;18;3];
pars_hyp_num = length(pars_ind);
% Range for parameters
pars_range_lower = [0  ;0  ;0  ;0  ;0  ;0  ]/100;
pars_range_upper = [200;200;100;100;100;100]/100;

pars_names = {'$K_{bar}$', '$R_{bv}$'      , '$R_{aa-ss}$', ...
              '$N_{rs}$' , '$N_{als}^{eq}$', '$N_{rsna}$' };

% Scenario
scenario = {'Normal'};

% Species
spe_ind = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of variables.
num_vars = 93;

species = {'human', 'rat'   };
sex     = {'male' , 'female'};

for sex_ind = 1:2 % sex

%% Load bootstrap replicate parameters created by create_par_bs_rep.m.

load_data_name_pars = sprintf('%s_%s_pars_scenario_Pri_Hyp_bs_rep10.mat', ...
                              species{spe_ind},sex{sex_ind});
load(load_data_name_pars, 'pars_rep');
pars_hyp = pars_rep(pars_ind,:);

pars_num   = size(pars_rep,1);
num_sample = size(pars_rep,2);
% num_sample = 10;

%% Plot parameter distribution.

% Plot
f1  = figure('DefaultAxesFontSize',14);
s1 = gobjects(pars_hyp_num);
for i = 1:pars_hyp_num
    s1(i) = subplot(3,2,i);
    histogram(s1(i),pars_hyp(i,:))
    xlabel(s1(i), pars_names(i), 'Interpreter','latex', 'FontSize',16)
end
hist_title = sprintf('%s',sex{sex_ind});
sgtitle(hist_title, 'FontSize',16)

%% Create virtual population corresponding to each parameter set.

tic
SSdata_rep = zeros(num_vars, num_sample);
for j = 1:num_sample
    SSdata_rep(:,j) = solve_ss_scenario(pars_rep(:,j));
    fprintf('%s iteration = %s out of %s \n', ...
            sex{sex_ind},num2str(j),num2str(num_sample))
end
bs_rep_solve_time = toc

%% Plot variables of interest distribution.

vars_ind   = [42;33;41;29;30;52;6;7;92];
vars_names = {'$P_{ma}$'   , '$\Phi_{co}$'   , '$R_{tp}$'  , ...
              '$V_{ecf}$'  , '$V_{b}$'       , '$C_{sod}$' , ...
              '$\Phi_{rb}$', '$\Phi_{gfilt}$', '$\Phi_{u}$'};
vars_hyp_num = length(vars_ind);
vars_hyp = SSdata_rep(vars_ind,:);

% Plot
f2  = figure('DefaultAxesFontSize',14);
s2 = gobjects(vars_hyp_num);
for i = 1:vars_hyp_num
    s2(i) = subplot(3,3,i);
    histogram(s2(i),vars_hyp(i,:))
    xlabel(s2(i), vars_names(i), 'Interpreter','latex', 'FontSize',16)
end
hist_title = sprintf('%s',sex{sex_ind});
sgtitle(hist_title, 'FontSize',16)

%% Save data.

% save_data_name = sprintf('%s_%s_ss_data_scenario_Pri_Hyp_bs_rep10.mat', ...
%                          species{spe_ind},sex{sex_ind});
% save_data_name = strcat('Data/', save_data_name);
% save(save_data_name, 'SSdata_rep', 'num_sample', 'bs_rep_solve_time')

end % sex

% -------------------------------------------------------------------------
% Steady state solver function
% -------------------------------------------------------------------------

function SSdata = solve_ss_scenario(pars)

%% Parameters

varargin_input = {scenario{1},true};

%% Variables initial guess

% Load data for steady state initial guess. 
% Set name for data file to be loaded based upon sex.

load_data_name_SSdata = sprintf('%s_%s_ss_data_scenario_%s.mat', ...
                                species{spe_ind},sex{sex_ind},scenario{1});
load(load_data_name_SSdata, 'SSdata');
SSdataIG     = SSdata;
clear SSdata

% Initial guess for the variables.
% Find the steady state solution, so the derivative is 0.
% Arbitrary value for time to input, greater than tchange + deltat.
x0 = SSdataIG; x_p0 = zeros(num_vars,1); t = 30;

% Time at which to change place holder.
tchange = 0;

%% Find steady state solution

options = optimset('Display','off');
[SSdata  , ~     , ...
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
for vari = 1:length(SSdata)
    if abs(SSdata(vari)) < eps*100
        SSdata(vari) = 0;
    end
end

end % solve_ss_scenario

end % create_vp






























































