% This simulates the blood pressure regulation model blood_press_reg.m.
% 
% Parameters are given by:
% "Long-Term Mathematical Model Involving Renal Sympathetic Nerve Activity,
% Arterial Pressure, and Sodium Excretion" - 2005 - Karaaslan, et. al.
% "Sex-specific Long-term Blood Pressure Regulation: Modeling and Analysis"
% - 2018 - Leete, Layton.
% 
% Steady state data is calculated by solve_ss_numerical.m.

function run_sim_treatments(human,gg,IC, varargin)
species = {'rat','human'};
gender     = {'male', 'female'};

%% Define default imputs
AA = 1;
ACEI = 0;
furosemide = 0;
NSAID = 0;
myo_ind = 0;
water_ind = 0;

%% Read and assign optional variables
for i = 1:2:length(varargin)
    if strcmp(varargin{i},'ACEI')
        ACEI = varargin{i + 1}; %ACEI indicator
    elseif strcmp(varargin{i},'furosemide')
        furosemide = varargin{i + 1}; %furosemide indicator
    elseif strcmp(varargin{i},'NSAID')
        NSAID = varargin{i+1}; %indicator
   elseif strcmp(varargin{i},'Myogenic Response')
        myo_ind = varargin{i+1}; %indicator 0 for normal, 1 for impaired
   elseif strcmp(varargin{i},'Water Intake')
        water_ind = varargin{i+1};%indicator 0 for normal, 1 for low
    elseif strcmp(varargin{i},'RSNA')
        AA = varargin{i+1}; %multiply N_rsna by to simulate hypertension
    end
end

pars = get_pars(species{human+1},gender{gg},'',AA);


kappa_ACEI = 0;
kappa_f = 0;
kappa_f_md = 0;

   if ACEI == 1
       kappa_ACEI = 0.76;
   elseif ACEI > 1
       kappa_ACEI = 0.90;
   end
   if furosemide == 1
       kappa_f = 0.15;
       kappa_f_md = 0.4;
   elseif furosemide == 2
       kappa_f = 0.3;
       kappa_f_md = 0.5;
   end  

%% Solve DAE

% Initial value
% This initial condition is the steady state data value taken from
% experiments (CITE). Therefore, the initial condition of the derivative is
% 0.


% Load data for steady state initial value. 
load(IC,'SSdata');
SS_data_IC = SSdata;

% Initial condition for the variables and their derivatives. 
% System is initially at steady state, so the derivative is 0.
x0 = SS_data_IC; x_p0 = zeros(84,1);

% Time at which to keep steady state, change a parameter, etc.
tchange = 100;

% Initial time (min); Final time (min);\
num_days =2;
t0 = 0*1440; tf = tchange + num_days*1440;
% Time vector; %if any simulation fails, try adjusting this
%tspan = [t0, tf];
tspan = t0:10:tf;
% Solve dae
[t,x] = ode15i(@(t,x,x_p) bp_reg_mod(t,x,x_p,pars,tchange,...
                                     'ACEI',kappa_ACEI,'furosemide',[kappa_f,kappa_f_md],'NSAID',NSAID',...
                                     'Myogenic Response',myo_ind,'Water Intake',water_ind), ...
                tspan, x0, x_p0);

if not (isreal(x))
    disp('Imaginary number returned.')
end

% %Save Data
save_name_ending = '';
if water_ind
    save_name_ending = strcat(save_name_ending,'_lowwaterintake');
end
if myo_ind
    save_name_ending = strcat(save_name_ending,'_impairedmyo');
end


filename = '%s%s_%s_%s_%s_%s_%s%s.mat';


filename = sprintf( filename,'Data/',species{human+1},gender{gg},num2str(ACEI),num2str(furosemide),num2str(NSAID),num2str(AA),save_name_ending);
%filename
save(filename, 't','x');

end

































