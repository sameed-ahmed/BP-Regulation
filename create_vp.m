% This script loads the bootstrap replicate parameter sets created for the
% virtual population. It then solves the corresponding steady state
% solution of the system for each parameter set and saves the result.

% Input:  none
% Output: saves steady state variables values corresponding to parameter
% bootstrap replicate set.

function create_vp

close all

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

% Parameters to perturb
% K_bar            - par 13; - 0%, +200%
% R_bv             - par 14; - 0%, +200%
% R_aass           - par 4 ; - 0%, +100%
% N_rs             - par 21; - 0%, +100%
% N_als_eq         - par 18; - 0%, +100%
% N_rsna           - par 3 ; - 0%, +100%
% N_adhs_eq        - par 15; - 0%, +100%
% sigmamyo_b       - par 41; - 0%, +900%
% Indices
pars_ind = [13;14;4;21;18;3;15;41];

% Scenario
scenario = {'Normal'};
% Index of scenario to fix.
fixed_ss = 1;

% Species
spe_ind = 2;

% Number of variables.
num_vars = 93;

species = {'human', 'rat'   };
sex     = {'male' , 'female'};

for sex_ind = 1:2 % sex

%% Load bootstrap replicate parameters created by create_par_bs_rep.m.

load_data_name_pars = sprintf('%s_%s_pars_scenario_Pri_Hyp_bs_rep1000.mat', ...
                              species{spe_ind},sex{sex_ind});
load(load_data_name_pars, 'pars_rep');
pars_hyp = pars_rep(pars_ind,:);

% pars0_est = [24.5664; 11.9122; 3.2875; 1.9027; 1.9448; 1.4909; 1.4893; 4.8474]; % diverges for j = 076
% pars_rep(pars_ind,sample_num) = pars0_est;

% Number of bootstrap replicates.
num_sample = size(pars_rep,2);

%% Load data for steady state initial guess. 

% Set name for data file to be loaded based upon sex.
load_data_name_SSdata = sprintf('%s_%s_ss_data_scenario_%s.mat', ...
                                species{spe_ind},sex{sex_ind},scenario{fixed_ss});
load(load_data_name_SSdata, 'SSdata');
SSdataIG     = SSdata;
clear SSdata

%% Create virtual population corresponding to each parameter set.

tic
SSdata_rep = zeros(num_vars, num_sample);
for j = 1:num_sample
    SSdata_rep(:,j) = solve_ss_scenario(pars_rep(:,j));
    fprintf('%s iteration = %s out of %s \n', ...
            sex{sex_ind},num2str(j),num2str(num_sample))
end
bs_rep_solve_time = toc

%% Save data.

save_data_name = sprintf('%s_%s_ss_data_scenario_Pri_Hyp_bs_rep1000.mat', ...
                         species{spe_ind},sex{sex_ind});
save_data_name = strcat('Data/', save_data_name);
save(save_data_name, 'SSdata_rep', 'num_sample', 'bs_rep_solve_time')

end % sex

% -------------------------------------------------------------------------
% Steady state solver function
% -------------------------------------------------------------------------

function SSdata = solve_ss_scenario(pars)

%% Parameters

varargin_input = {scenario{fixed_ss},true};

%% Variables initial guess

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



























