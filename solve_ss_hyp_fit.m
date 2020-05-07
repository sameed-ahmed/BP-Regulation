% This script calculates the steady state values of the blood pressure 
% regulation model bp_reg_solve_baseline.m by using fsolve.
% It is adopted with modifications from Karaaslan 2005 and Leete 2018.
% 
% Steady state data for the intial guess is inputted by solver_initial_guess_data.m.

function solve_ss_hyp_fit

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
par_ind = [13;14;4;21;18;3];
par_num = length(par_ind);
% Range for parameters
par_range_lower = [0  ;0  ;0  ;0  ;0  ;0  ]/100;
par_range_upper = [200;200;100;100;100;100]/100;

% Variables to check
% P_ma      - var 42; +30,40, +40,50
% Phi_co    - var 33; -5%   , +5%
% Phi_rb    - var 6 ; -5%   , +5%
% Phi_gfilt - var 7 ; -5%   , +5%
% Phi_u     - var 92; -5%   , +5%
% Phi_usod  - var 27; -5%   , +5%
% C_sod     - var 52; -2%   , +2%
% Indices
var_ind = [42;33;6;7;92;27;52];
num_vars_check = length(var_ind);
% Range for variables for each sex.
var_range_lower_change_m = [40*100;5;5;5;5;5;2]/100;
var_range_upper_change_m = [50*100;5;5;5;5;5;2]/100;
var_range_lower_change_f = [30*100;5;5;5;5;5;2]/100;
var_range_upper_change_f = [40*100;5;5;5;5;5;2]/100;

% Scenarios
% Normal - Normal conditions
% m_RAS  - male RAS pars
% m_Reab - male fractional sodium and water reabsorption
scenario = {'Normal', 'm_RAS', 'm_Reab', 'm_RAS_m_Reab'};
% Index of scenario to fix.
fixed_ss = 1;

% Species
spe_ind = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of variables; number of parameters; 
num_vars = 93; num_pars = 47; % + SF + fixed_var_pars + SSdata

species = {'human', 'rat'   };
sex     = {'male' , 'female'};

parpool
for sex_ind = 1:2 % sex

%% Parameters

varargin_input = {scenario{fixed_ss},true};

% Parameter input
pars0 = get_pars(species{spe_ind}, sex{sex_ind}, varargin_input{:});

% Set interval bounds for parameters.
lower = pars0(par_ind) - par_range_lower .* pars0(par_ind);
upper = pars0(par_ind) + par_range_upper .* pars0(par_ind);

%% Drugs

% drugs = [Ang II inf rate fmol/(ml min), ACEi target level, ARB target level]
if     strcmp(scenario{fixed_ss}, 'AngII')
    if     strcmp(sex{sex_ind}, 'male'  )
%         varargin_input = {'AngII',2022}; % Sampson 2008
        varargin_input = {'AngII',910 }; % Sullivan 2010
    elseif strcmp(sex{sex_ind}, 'female')
%         varargin_input = {'AngII',2060}; % Sampson 2008
        varargin_input = {'AngII',505 }; % Sullivan 2010
    end
elseif strcmp(scenario{fixed_ss}, 'ACEi' )
        varargin_input = {'ACEi' ,0.78 }; % Leete 2018
elseif strcmp(scenario{fixed_ss}, 'ARB'  )
        varargin_input = {'ARB'  ,0.67 }; % Leete 2018
end

%% Variables initial guess

% Load data for steady state initial guess. 
% Set name for data file to be loaded based upon sex.    
load_data_name = sprintf('%s_%s_ss_data_scenario_%s.mat', ...
                         species{spe_ind},sex{sex_ind},scenario{fixed_ss});
load(load_data_name, 'SSdata');
SSdataIG = SSdata;

% Set acceptable range for certain variables.
% MAP needs special treatment because it is an additive range, 
% not a relative range.
var_bl_value    = SSdata(var_ind(2:end));
if     strcmp(sex{sex_ind},'male')
    var_range_lower_change = var_range_lower_change_m;
    var_range_upper_change = var_range_upper_change_m;
elseif strcmp(sex{sex_ind},'female')
    var_range_lower_change = var_range_lower_change_f;
    var_range_upper_change = var_range_upper_change_f;
end
var_range_lower = var_bl_value - var_range_lower_change(2:end) .* var_bl_value;
var_range_upper = var_bl_value + var_range_upper_change(2:end) .* var_bl_value;
MAP_lower       = SSdata(var_ind(1)) + var_range_lower_change(1);
MAP_upper       = SSdata(var_ind(1)) + var_range_upper_change(1);
var_range_lower = [MAP_lower; var_range_lower];
var_range_upper = [MAP_upper; var_range_upper];

clear SSdata

% Order
% x  = [rsna; alpha_map; alpha_rap; R_r; beta_rsna; Phi_rb; Phi_gfilt; ...
%       P_f; P_gh; Sigma_tgf; Phi_filsod; Phi_ptsodreab; eta_ptsodreab; ...
%       gamma_filsod; gamma_at; gamma_rsna; Phi_mdsod; Phi_dtsodreab; ...
%       eta_dtsodreab; psi_al; Phi_dtsod; Phi_cdsodreab; eta_cdsodreab; ...
%       lambda_dt; lambda_anp; lambda_al; Phi_usod; Phi_sodin; V_ecf; ...
%       V_b; P_mf; Phi_vr; Phi_co; P_ra; vas; vas_f; vas_d; R_a; R_ba; ...
%       R_vr; R_tp; P_ma; epsilon_aum; a_auto; a_chemo; a_baro; C_adh; ...
%       N_adh; N_adhs; delta_ra; ...
%       M_sod; C_sod; nu_mdsod; nu_rsna; C_al; N_al; N_als; xi_ksod; ...
%       xi_map; xi_at; hatC_anp; AGT; nu_AT1; R_sec; PRC; PRA; AngI; ...
%       AngII; AT1R; AT2R; Ang17; AngIV; R_aa; R_ea; Sigma_myo; ...
%       Psi_AT1RAA; Psi_AT1REA; Psi_AT2RAA; Psi_AT2REA; ...
%       Phi_ptwreab; eta_ptwreab; mu_ptsodreab; ...
%       Phi_mdu; Phi_dtwreab; eta_dtwreab; mu_dtsodreab; Phi_dtu; ...
%       Phi_cdwreab; eta_cdwreab; mu_cdsodreab; mu_adh; Phi_u; Phi_win];

% Initial guess for the variables.
% Find the steady state solution, so the derivative is 0.
% Arbitrary value for time to input, greater than tchange + deltat.
x0 = SSdataIG; x_p0 = zeros(num_vars,1); t = 30;

% Time at which to change place holder.
tchange = 0;

%% Uniformly randomly sample perturbed parameters from range.

% Sample between 0 and 1.
ran_vec = random('unif',0,1,par_num,1);
% Sample within interval.
ran_vec = lower + ran_vec .* (upper - lower);
% Replace input parameters with newly sampled parameter.
pars0(par_ind) = ran_vec;
pars0_est = ran_vec;

% Initialize last pars to avoid unnecessary computation for SSdata.
pars_est_last = [];
% Initialize SSdata to avoid unnecessary computation for SSdata.
SSdata_iter = [];

%% Optimize.

% pars0
% spe_par = pars0(1);
% sex_par = pars0(2);

% Place holders for fmincon.
A = []; b = []; Aeq = []; beq = []; 
% Nonlinear constraints.
nonlcon = @mycon;
% Lower and upper bounds for parameters in fmincon.
lb = lower;
ub = upper;

tic
% % Edit options for optimizer. - fmincon
% opt_name = 'fm';
% options = optimoptions('fmincon', 'Display','iter', 'UseParallel',true);
% % [pars_est_min, residual_pars, exitflag_pars, output_pars] = ...
% %     fmincon(@(pars_est) ...
% %             cost_fun(t,x0,x_p0,pars0,pars_est,par_ind,tchange,varargin_input, ...
% %                      var_ind,var_range_lower,var_range_upper), ...
% %             pars0_est,A,b,Aeq,beq,lb,ub,nonlcon,options); % %#ok<ASGLU>
% [pars_est_min, residual_pars, exitflag_pars, output_pars] = ...
%     fmincon(@cost_fun,pars0_est,A,b,Aeq,beq,lb,ub,nonlcon,options); % %#ok<ASGLU>

% % Edit options for optimizer. - MultiStart
% opt_name = 'ms';
% N_ms = 5;
% options = optimoptions('fmincon', 'Display','iter', 'UseParallel',true);
% ms = MultiStart;
% % options = optimoptions('fmincon', 'Display','iter');
% % ms = MultiStart('UseParallel',true);
% % problem = ...
% %     createOptimProblem('fmincon','x0',pars0_est,...
% %                        'objective',@(pars_est) ...
% %                                    cost_fun(t,x0,x_p0,pars0,pars_est,par_ind,tchange,varargin_input, ...
% %                                             var_ind,var_range_lower,var_range_upper), ...
% %                        'lb',lb,'ub',ub,'nonlcon',nonlcon,'options',options);
% problem = ...
%     createOptimProblem('fmincon', 'x0',pars0_est, 'objective',@cost_fun, ...
%                        'lb',lb,'ub',ub,'nonlcon',nonlcon,'options',options);
% [pars_est_min, residual_pars, exitflag_pars, output_pars, solutions] = run(ms,problem,N_ms);

% % Edit options for optimizer. - GlobalSerach
% opt_name = 'gs';
% options = optimoptions('fmincon', 'Display','iter', 'UseParallel',true);
% gs = GlobalSearch('StartPointsToRun','bounds-ineqs', ...
%                   'MaxWaitCycle',3,'BasinRadiusFactor',0.3);
% % problem = ...
% %     createOptimProblem('fmincon','x0',pars0_est,...
% %                        'objective',@(pars_est) ...
% %                                    cost_fun(t,x0,x_p0,pars0,pars_est,par_ind,tchange,varargin_input, ...
% %                                             var_ind,var_range_lower,var_range_upper), ...
% %                        'lb',lb,'ub',ub,'nonlcon',nonlcon,'options',options);
% problem = ...
%     createOptimProblem('fmincon', 'x0',pars0_est, 'objective',@cost_fun, ...
%                        'lb',lb,'ub',ub,'nonlcon',nonlcon,'options',options);
% [pars_est_min, residual_pars, exitflag_pars, output_pars, solutions] = run(gs,problem);

% Edit options for optimizer. - pattersearch
opt_name = 'ps';
% options = optimoptions('patternsearch', 'Display','iter', 'UseCompletePoll',true);
options = optimoptions('patternsearch', 'Display','iter', 'UseCompletePoll',true, ...
                       'UseParallel',true, 'UseVectorized',false);
% [pars_est_min, residual_pars, exitflag_pars, output_pars] = ...
%     patternsearch(@(pars_est) ...
%                   cost_fun(t,x0,x_p0,pars0,pars_est,par_ind,tchange,varargin_input, ...
%                            var_ind,var_range_lower,var_range_upper), ...
%                   pars0_est,A,b,Aeq,beq,lb,ub,nonlcon,options);
[pars_est_min, residual_pars, exitflag_pars, output_pars] = ...
    patternsearch(@cost_fun,pars0_est,A,b,Aeq,beq,lb,ub,nonlcon,options);

% % Edit options for optimizer. - ga
% opt_name = 'ga';
% options = optimoptions('ga', 'Display','iter', 'UseParallel',true, 'UseVectorized',false);
% % [pars_est_min, residual_pars, exitflag_pars, output_pars] = ...
% %     ga(@(pars_est) ...
% %        cost_fun(t,x0,x_p0,pars0,pars_est,par_ind,tchange,varargin_input, ...
% %                 var_ind,var_range_lower,var_range_upper), ...
% %        length(pars0_est),A,b,Aeq,beq,lb,ub,nonlcon,options);
% [pars_est_min, residual_pars, exitflag_pars, output_pars] = ...
%     ga(@cost_fun,par_num,A,b,Aeq,beq,lb,ub,nonlcon,options);

% % Edit options for optimizer. - simulannealbnd
% opt_name = 'sa';
% options = optimoptions('simulannealbnd', 'Display','iter');
% % [pars_est_min, residual_pars, exitflag_pars, output_pars] = ...
% %     simulannealbnd(@(pars_est) ...
% %                    cost_fun(t,x0,x_p0,pars0,pars_est,par_ind,tchange,varargin_input, ...
% %                             var_ind,var_range_lower,var_range_upper), ...
% %                    pars0_est,lb,ub,options);
% [pars_est_min, residual_pars, exitflag_pars, output_pars] = ...
%     simulannealbnd(@cost_fun,pars0_est,lb,ub,options);
opt_time = toc

% % tic
% test1 = cost_fun(pars0_est);
% % toc1 = toc
% % tic
% test2 = mycon   (pars0_est);
% % toc1 = toc

% Place estimated pars in proper location.
pars = pars0;
pars(par_ind) = pars_est_min;

% pars_min
% spe_par = pars_min(1);
% sex_par = pars_min(2);

% Solve system with found pars.
options2 = optimset();
[SSdata, residual_ss, exitflag_ss, output_ss] = ...
    fsolve(@(x) ...
           bp_reg_mod(t,x,x_p0,pars,tchange,varargin_input{:}), ...
           x0, options2); % %#ok<ASGLU>

%% Save values.

% % Steady state data
% % save_data_name = sprintf('%s_%s_ss_data_scenario_Pri_Hyp_%s.mat', ...
% %                          species{spe_ind},sex{sex_ind},opt_name);
% save_data_name = sprintf('%s_%s_ss_data_scenario_Pri_Hyp.mat', ...
%                          species{spe_ind},sex{sex_ind});
% save_data_name = strcat('Data/', save_data_name);
% save(save_data_name, 'SSdata', 'residual_ss', 'exitflag_ss', 'output_ss')
% % Parameters
% % save_data_name = sprintf('%s_%s_pars_scenario_Pri_Hyp_%s.mat', ...
% %                          species{spe_ind},sex{sex_ind},opt_name);
% save_data_name = sprintf('%s_%s_pars_scenario_Pri_Hyp.mat', ...
%                          species{spe_ind},sex{sex_ind});
% save_data_name = strcat('Data/', save_data_name);
% if strcmp(opt_name, 'ms') || strcmp(opt_name, 'gs')
%     save(save_data_name, 'pars', 'solutions', 'residual_pars', 'exitflag_pars', 'output_pars', 'opt_time') %#ok<USESWNS>
% else
%     save(save_data_name, 'pars',              'residual_pars', 'exitflag_pars', 'output_pars', 'opt_time')
% end

end % sex
delete(gcp)

% -------------------------------------------------------------------------
% Cost function
% -------------------------------------------------------------------------

function tot_err = cost_fun(pars_est)

% Total error
% alpha = 0.0; % beta  = 2.0 - alpha;
% tot_err = (1.0*range_err + 1.0*AngII_MAP_err) / 2;
% tot_err = (0.0*range_err + 2.0*AngII_MAP_err) / 2;
% tot_err = AngII_err(pars_est);
% tot_err = range_err;
tot_err = (AngII_err(pars_est) + Sodin_err(pars_est)) / 2;

end % tot err

% -------------------------------------------------------------------------
% Ang II inf error
% -------------------------------------------------------------------------

function err = AngII_err(pars_est)

% Place estimated pars in proper location.
pars0(par_ind) = pars_est;

%% Find steady state solution ---------------------------------------------

% tic
% Check if computation is necessary for SSdata.
if ~isequal(pars_est,pars_est_last)
    options_ss = optimset('Display','off');
    [SSdata_iter, ~, ~, ~] = ...
        fsolve(@(x) ...
               bp_reg_mod(t,x,x_p0,pars0,tchange,varargin_input{:}), ...
               x0, options_ss);
    pars_est_last = pars_est;
end
% toc1 = toc

% % Check for solver convergence.
% if exitflag == 0
%     tot_err = 1;
%     return
% end
% 
% % Check for imaginary solution.
% if not (isreal(SSdata))
%     tot_err = 1;
%     return
% end
% 
% % Set any values that are within machine precision of 0 equal to 0.
% for i = 1:length(SSdata)
%     if abs(SSdata(i)) < eps*100
%         SSdata(i) = 0;
%     end
% end
% 
% % Compute error within range of specified variables.
% num_vars = length(var_ind);
% range_err = zeros(num_vars,1);
% for i = 1:num_vars
%     range_err(i) = max( ( (SSdata(var_ind(i)) - (var_range_lower(i) + var_range_upper(i))/2)^2 - ...
%                     (                     (var_range_upper(i) - var_range_lower(i))/2)^2 ) ...
%                   / SSdata(var_ind(i))^2, 0 );
% end
% 
% % Range error
% range_err = sqrt(sum(range_err)) / num_vars;

%% Run Ang II infusion simulation. ----------------------------------------

% % Retrieve species and sex identifier. 
% spe_par = pars0(1);
% sex_par = pars0(2);
% if     spe_par == 1
%     species_iter = 'human';
% elseif spe_par == 0
%     species_iter = 'rat';
% end
% if     sex_par == 1
%     sex_iter = 'male';
% elseif sex_par == 0
%     sex_iter = 'female';
% end

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
%     kappa_AngII = 707; % Sullivan 2010
elseif strcmp(sex{sex_ind}, 'female')
%     kappa_AngII = 2060; % Sampson 2008
    kappa_AngII = 505; % Sullivan 2010
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

% Initial value is steady state solution with given pars.
x0_angII = SSdata_iter;

% ode options
options_dae = odeset('MaxStep',1000); % default is 0.1*abs(t0-tf)
% tic
% Solve dae
[~,x] = ode15i(@(t,x,x_p) ...
               bp_reg_mod(t,x,x_p,pars0,tchange_angII,varargin_input_angII{:}), ...
               tspan, x0_angII, x_p0, options_dae);
% toc2 = toc

% Return if simulation crashed.
if size(x,1) < N
    err = 5;
    return
end
% X = (variables, points)
X(:,:) = x';

%% Compute error.

% Data from Sullivan 2010. MAP is in difference from baseline.
tdata       = [0+1 , 1+1 , 2+1 , 3+1 , 4+1 , 5+1 , 6+1 ,...
               7+1 , 8+1 , 9+1 , 10+1, 11+1, 12+1, 13+1, 14+1]*1 + 1;
if     strcmp(sex{sex_ind}, 'male'  )
    MAPdata = [0   , 1.1 , 2.3 , 8.9 , 15.5, 18.3, 22.7, 22.6, ...
               28.6, 31.2, 30.9, 32.8, 37.4, 41.4, 40.3];
elseif strcmp(sex{sex_ind}, 'female')
    MAPdata = [0   , 5.2 ,  5.3,  3.9,  3.6,  5.9,    8,   13, ...
               15.7, 17.4, 19.8, 23.7, 25.8,  23.5,  24];
end
num_points = length(tdata);

% tic
% Substract MAP by baseline.
% X = (variable, points)
MAP = X(42,tdata) - X(42,1);

% Ang II inf error
AngII_MAP_err        = (MAP - MAPdata).^2;
AngII_MAP_err(2:end) = AngII_MAP_err(2:end) ./ MAPdata(2:end).^2;
% AngII_MAP_err        = sqrt(sum(AngII_MAP_err(8:end)) / (num_points-7));
AngII_MAP_err        = sqrt(mean(AngII_MAP_err(8:end)));
% toc3 = toc

% Error
% alpha = 0.0; % beta  = 2.0 - alpha;
% tot_err = (1.0*range_err + 1.0*AngII_MAP_err) / 2;
% tot_err = (0.0*range_err + 2.0*AngII_MAP_err) / 2;
err = AngII_MAP_err;
% tot_err = range_err;

end % Ang II err

% -------------------------------------------------------------------------
% Sodiun intake error
% -------------------------------------------------------------------------

function err = Sodin_err(pars_est)

% Place estimated pars in proper location.
pars0(par_ind) = pars_est;

% 4-fold increase in sodium intake.
pars_sodin = pars0;
pars_sodin(17) = 4 * pars_sodin(17);

%% Find steady state solution ---------------------------------------------

% tic
% Check if computation is necessary for baseine SSdata.
if ~isequal(pars_est,pars_est_last)
    options_ss = optimset('Display','off');
    [SSdata_iter, ~, ~, ~] = ...
        fsolve(@(x) ...
               bp_reg_mod(t,x,x_p0,pars0,tchange,varargin_input{:}), ...
               x0, options_ss);
    pars_est_last = pars_est;
end
% toc1 = toc

% tic
% Compute SSdata for increased sodium intake.
options_ss = optimset('Display','off');
[SSdata_sodin, ~, exitflag, ~] = ...
    fsolve(@(x) ...
           bp_reg_mod(t,x,x_p0,pars_sodin,tchange,varargin_input{:}), ...
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
    MAPdata = [15,25];
%     MAPdata = [0,100];
%     MAPdata = 20;
elseif strcmp(sex{sex_ind}, 'female')
    MAPdata = [ 5,10];
%     MAPdata = [0,100];
%     MAPdata = 7;
end

% tic
% Substract MAP by baseline.
% X = (variable, points)
MAP = SSdata_sodin(42) - SSdata_iter(42);

% Sodin error
err = max( ( (MAP - (MAPdata(1) + MAPdata(2))/2)^2 -  ...
             (      (MAPdata(2) - MAPdata(1))/2)^2 ), ...
         0 );
err = err / mean(MAPdata)^2;
err = sqrt(err);
% err = abs(MAP - MAPdata) / MAPdata;
% toc3 = toc

end % Sodin err

% -------------------------------------------------------------------------
% Nonlinear constraints
% -------------------------------------------------------------------------

function [c,ceq] = mycon(pars_est)

% Place estimated pars in proper location.
pars0(par_ind) = pars_est;

%% Find steady state solution ---------------------------------------------

% tic
% Check if computation is necessary for SSdata.
if ~isequal(pars_est,pars_est_last)
    options_ss = optimset('Display','off');
    [SSdata_iter, ~, exitflag, ~] = ...
        fsolve(@(x) ...
               bp_reg_mod(t,x,x_p0,pars0,tchange,varargin_input{:}), ...
               x0, options_ss);
    pars_est_last = pars_est;
           
    % Check for solver convergence.
    if exitflag == 0
        c   = 1;
        ceq = [];
        return
    end

    % Check for imaginary solution.
    if not (isreal(SSdata_iter))
        c   = 1;
        ceq = [];
        return
    end
end
% toc1 = toc

% % Set any values that are within machine precision of 0 equal to 0.
% for i = 1:length(SSdata)
%     if abs(SSdata(i)) < eps*100
%         SSdata(i) = 0;
%     end
% end

% tic
% Nonlinear inequalities.
% num_vars_check = length(var_ind);
c = zeros(2*num_vars_check,1);
for i = 1:num_vars_check
    c(2*i-1) =    SSdata_iter(var_ind(i)) - var_range_upper(i)  ;
    c(2*i)   = -( SSdata_iter(var_ind(i)) - var_range_lower(i) );
end
% toc2 = toc

% Nonlinear equalities.
ceq = [];

end % mycon

end % solve_ss_hyp_fit





























