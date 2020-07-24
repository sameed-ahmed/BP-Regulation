
function compute_err_hyp_fit

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

fprintf('------------------------------------------ \n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Begin user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scenarios
% Normal - Normal conditions
% m_RAS  - male RAS pars
% m_Reab - male fractional sodium and water reabsorption
scenario = {'Normal', 'm_RAS', 'm_Reab', 'm_RAS_m_Reab'};
% Index of scenario to fix.
fixed_ss = 1;

% Dataset
% dataset = 'old';
% dataset = 'new';
dataset = 'newnew';
% dataset = 'particular';
% dataset = 'particular range';

% Bootstrap replicate sample number
% sample_num = random('Discrete Uniform',1000)
sample_num = 974
% sample_num_range = 010;
% ---
% female bad fit for new1
% sample_num = 239 
% sample_num = 723 
% sample_num = 261 
% sample_num = 823 
% sample_num = 098 
%   male bad fit for new1 ~
% sample_num = 975 
% sample_num = 080 
%   both bad fit for new1 ~
% sample_num = 003 
%  both good fit for new1
% sample_num = 318 
% sample_num = 951 
% sample_num = 035 
% ---
% female bad fit for old
% sample_num = 212 
% sample_num = 463 
% sample_num = 167 
%  both good fit for old
% sample_num = 206 
% sample_num = 573 
% sample_num = 742 
% sample_num = 077 
% sample_num = 207 
% sample_num = 862
% sample_num = 695
% ---
% new?
% sample_num = 001
% sample_num = 002
% sample_num = 003 % *      male
% sample_num = 004 % ***    female
% sample_num = 005 % *      female

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Species
spe_ind = 2;

% Number of variables; number of parameters; 
num_vars = 93; num_pars = 47; % + SF + fixed_var_pars + SSdata

species = {'human', 'rat'   };
sex     = {'male' , 'female'};

varargin_input = {scenario{fixed_ss},true};

for sex_ind = 1:2 % sex

%% Load bootstrap replicate parameters & variables created by create_par_bs_rep.m.

% Parameters
if     strcmp(dataset, 'old')
    load_data_name_pars = sprintf('%s_%s_pars_scenario_Pri_Hyp_bs_rep1000OLD.mat', ...
                                  species{spe_ind},sex{sex_ind});
elseif strcmp(dataset, 'new')
    load_data_name_pars = sprintf('%s_%s_pars_scenario_Pri_Hyp_bs_rep1000NEW.mat', ...
                                  species{spe_ind},sex{sex_ind});
elseif strcmp(dataset, 'newnew')
    load_data_name_pars = sprintf('%s_%s_pars_scenario_Pri_Hyp_bs_rep1000NEWNEW.mat', ...
                                  species{spe_ind},sex{sex_ind});
elseif strcmp(dataset, 'particular')
    load_data_name_pars = sprintf('%s_%s_pars_scenario_Pri_Hyp_bs_rep%s.mat', ...
                                  species{spe_ind},sex{sex_ind},num2str(sample_num));
elseif strcmp(dataset, 'particular range')
    load_data_name_pars = sprintf('%s_%s_pars_scenario_Pri_Hyp_bs_rep%s.mat', ...
                                  species{spe_ind},sex{sex_ind},num2str(sample_num_range));
end
load(load_data_name_pars, 'pars_rep');

% Variables
if     strcmp(dataset, 'old')
    load_data_name_vars = sprintf('%s_%s_ss_data_scenario_Pri_Hyp_bs_rep1000OLD.mat', ...
                                  species{spe_ind},sex{sex_ind});
elseif strcmp(dataset, 'new')
    load_data_name_vars = sprintf('%s_%s_ss_data_scenario_Pri_Hyp_bs_rep1000NEW.mat', ...
                                  species{spe_ind},sex{sex_ind});
elseif strcmp(dataset, 'newnew')
    load_data_name_vars = sprintf('%s_%s_ss_data_scenario_Pri_Hyp_bs_rep1000NEWNEW.mat', ...
                                  species{spe_ind},sex{sex_ind});
elseif strcmp(dataset, 'particular')
load_data_name_vars = sprintf('%s_%s_ss_data_scenario_Pri_Hyp_bs_rep%s.mat', ...
                              species{spe_ind},sex{sex_ind},num2str(sample_num));
elseif strcmp(dataset, 'particular range')
load_data_name_vars = sprintf('%s_%s_ss_data_scenario_Pri_Hyp_bs_rep%s.mat', ...
                              species{spe_ind},sex{sex_ind},num2str(sample_num_range));
end
load(load_data_name_vars, 'SSdata_rep');

%% Load bootstrap replicate data created by create_data_bs_rep.m.

load_data_name_data = sprintf('%s_%s_AngII_data_bs_rep.mat', ...
                              species{spe_ind},sex{sex_ind});
load(load_data_name_data, 'AngII_data_rep');
AngII_MAP_data = AngII_data_rep(sample_num,:);

%% Variables initial value

SSdataIG = SSdata_rep(:,sample_num);

% Initial guess for the variables.
% Find the steady state solution, so the derivative is 0.
% Arbitrary value for time to input, greater than tchange + deltat.
x0 = SSdataIG; x_p0 = zeros(num_vars,1); t = 30;

% Time at which to change place holder.
tchange = 0;

%% Compute error

fprintf('%s AngII error \n', sex{sex_ind})
AngII_err(pars_rep(:,sample_num))
fprintf('%s Sodin error \n', sex{sex_ind})
Sodin_err(pars_rep(:,sample_num))

end

%% Sub functions

% -------------------------------------------------------------------------
% Ang II inf error
% -------------------------------------------------------------------------

function err = AngII_err(pars)

%% Run Ang II infusion simulation. ----------------------------------------

% Initialization, etc.

% Number of days to run simulation after change; Day at which to induce change;
days = 14; day_change = 1;
% Number of points for plotting resolution
N = (days+1)*1 + 1;

% % Number of variables
% num_vars = 93;

% Initialize variables.
% X = (variables, points)
X = zeros(num_vars,N);

%% Input drugs.

% Ang II inf rate fmol/(ml min)
if     strcmp(sex{sex_ind}, 'male')
%     kappa_AngII = 2022; % Sampson 2008
    kappa_AngII = 910; % Sullivan 2010
%     kappa_AngII = 0; % none
%     kappa_AngII = 707; % Sullivan 2010
elseif strcmp(sex{sex_ind}, 'female')
%     kappa_AngII = 2060; % Sampson 2008
    kappa_AngII = 505; % Sullivan 2010
%     kappa_AngII = 0; % none
%     kappa_AngII = 707; % Sullivan 2010
end

varargin_input_angII = [varargin_input, 'AngII',kappa_AngII];

%% Solve DAE.

% Time at which to keep steady state, change a parameter, etc.
tchange_angII = day_change*1440;

% Initial time (min); Final time (min);
t0 = 0*1440; tend = tchange_angII + days*1440;

% Time vector
tspan = linspace(t0,tend,N);

% ode options
options_dae = odeset('MaxStep',1000); 
% default MaxStep is 0.1*abs(t0-tf)

% tic
% Solve dae
[~,x] = ode15i(@(t,x,x_p) ...
               bp_reg_mod(t,x,x_p,pars,tchange_angII,varargin_input_angII{:}), ...
               tspan, x0, x_p0, options_dae);

% X = (variables, points)
X(:,:) = x';

%% Compute error.

% Data from Sullivan 2010. MAP is in difference from baseline.
tdata       = [0+1 , 1+1 , 2+1 , 3+1 , 4+1 , 5+1 , 6+1 ,...
               7+1 , 8+1 , 9+1 , 10+1, 11+1, 12+1, 13+1, 14+1]*1 + 1;

% tic
% Substract MAP by baseline.
% X = (variable, points)
MAP = X(42,tdata) - X(42,1);

% Ang II inf error
AngII_MAP_err        = (MAP - AngII_MAP_data).^2;
AngII_MAP_err(2:end) = AngII_MAP_err(2:end) ./ AngII_MAP_data(2:end).^2;
% AngII_MAP_err        = sqrt(sum(AngII_MAP_err(8:end)) / (num_points-7));
% AngII_MAP_err        = sqrt(mean(AngII_MAP_err(8:end)));
AngII_MAP_err        = sqrt(mean(AngII_MAP_err(12:end)));
% toc3 = toc

% Error
err = AngII_MAP_err;

end % Ang II err

% -------------------------------------------------------------------------
% Sodium intake error
% -------------------------------------------------------------------------

function err = Sodin_err(pars)

% 4-fold increase in sodium intake.
pars(17) = 4 * pars(17);

%% Find steady state solution ---------------------------------------------

% tic
% Compute SSdata for increased sodium intake.
options_ss = optimset('Display','off');
[SSdata_sodin, ~, exitflag, ~] = ...
    fsolve(@(x) ...
           bp_reg_mod(t,x,x_p0,pars,tchange,varargin_input{:}), ...
           x0, options_ss);

% Check for solver convergence.
if exitflag == 0
    err = 5;
    return
end

% Check for imaginary solution.
if not (isreal(SSdata_sodin))
    err = 5;
    return
end
% toc2 = toc

%% Compute error.

% Data from several sources. See 'change in MAP.xlsx'.
if     strcmp(sex{sex_ind}, 'male'  )
    sodin_MAP_data = [15,25];
%     MAPdata = [0,100];
%     MAPdata = 20;
elseif strcmp(sex{sex_ind}, 'female')
    sodin_MAP_data = [ 5,10];
%     MAPdata = [0,100];
%     MAPdata = 7;
end

% tic
% Substract MAP by baseline.
% X = (variable, points)
MAP = SSdata_sodin(42) - x0(42);

% Sodin error
err = max( ( (MAP - (sodin_MAP_data(1) + sodin_MAP_data(2))/2)^2 -  ...
             (      (sodin_MAP_data(2) - sodin_MAP_data(1))/2)^2 ), ...
         0 );
err = err / mean(sodin_MAP_data)^2;
err = sqrt(err);
% err = abs(MAP - MAPdata) / MAPdata;
% toc3 = toc

end % Sodin err

end % solve_ss_hyp_fit





























