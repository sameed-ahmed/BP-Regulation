% This solves the time course simulations. It saves output to Human_Data/.
% 
% Required input parameters:
% human - 1 to run a human simulation, 0 for the rat simulation
% gg - 1 for male, 2 for female
% IC - string containing the filename of the initial condition for the solver

% Optional input parameters: (each should be a string followed by a value)
% 'ACEi', # - value between 0 and 1 for the variable kappa_ACEi. Default is 0.
% 'furosemide', # - array of length 2 containing the values between 0 and 1 for the variables kappa_f and kappa_f_md.
% 'NSAID', # - 0 for no treatment 1 for normal dose, 2 for high dose. Default is [0 0].
% 'Impaired Myogenic Response', # - 0 for normal, 1 for impaired. Default is 0.
% 'Low Water Intake', # - 0 for normal, 1 for low water intake. Default is 0.
% 'RSNA', # - value by which N_rsna is multiplied by to induce hypertension. In my simulations it is 1 for normotensive simulations, 2.5 for hypertensive. Default is 1.

function run_sim_treatments(human,gg,IC, varargin)
species = {'rat','human'};
gender     = {'male', 'female'};

%% Define default imputs
AA = 1;
ACEi = 0;
furosemide = 0;
NSAID = 0;
myo_ind = 0;
water_ind = 0;

%% Read and assign optional variables
for i = 1:2:length(varargin)
    if strcmp(varargin{i},'ACEi')
        ACEi = varargin{i + 1}; %ACEi indicator
    elseif strcmp(varargin{i},'furosemide')
        furosemide = varargin{i + 1}; %furosemide indicator
    elseif strcmp(varargin{i},'NSAID')
        NSAID = varargin{i+1}; %indicator
   elseif strcmp(varargin{i},'Impaired Myogenic Response')
        myo_ind = varargin{i+1}; %indicator 0 for normal, 1 for impaired
   elseif strcmp(varargin{i},'Low Water Intake')
        water_ind = varargin{i+1};%indicator 0 for normal, 1 for low
    elseif strcmp(varargin{i},'RSNA')
        AA = varargin{i+1}; %multiply N_rsna by to simulate hypertension
    end
end

pars = get_pars(species{human+1},gender{gg},'RSNA',AA);


kappa_ACEi = 0;
kappa_f = 0;
kappa_f_md = 0;

   if ACEi == 1
       kappa_ACEi = 0.76;
   elseif ACEi > 1
       kappa_ACEi = 0.90;
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
                                     'ACEi',kappa_ACEi,'furosemide',[kappa_f,kappa_f_md],'NSAID',NSAID',...
                                     'Impaired Myogenic Response',myo_ind,'Low Water Intake',water_ind), ...
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


filename = sprintf( filename,'Human_Data/',species{human+1},gender{gg},num2str(ACEi),num2str(furosemide),num2str(NSAID),num2str(AA),save_name_ending);
%filename
save(filename, 't','x');

end

































