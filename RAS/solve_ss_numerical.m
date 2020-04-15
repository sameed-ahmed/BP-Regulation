% This script calculates the steady state solution to the RAS system using
% fsolve.

function solve_ss_numerical

% num_vars = 9;
num_vars = 8;

species  = {'human', 'rat'};
gender   = {'male' , 'female'};
SS_DATA  = zeros(num_vars,2);
RESIDUAL = zeros(num_vars,2);
EXITFLAG = zeros(1 ,2);
OUTPUT   = cell (1 ,2);

for ss = 2:2 % species
for gg = 1:2 % gender

%% Parameters

if     strcmp(gender{gg}, 'male')
    gen = 1;
elseif strcmp(gender{gg}, 'female')
    gen = 0;
end

% RAS
h_renin  = 12;      % min
h_AGT    = 10*60;   % min
h_AngI   = 0.5;     % min
h_AngII  = 0.66;    % min
h_Ang17  = 30;      % min
h_AngIV  = 0.5;     % min
h_AT1R   = 12;      % min
h_AT2R   = 12;      % min

% Male and female different parameters for RAS
if     strcmp(species{ss}, 'human')
    X_PRCPRA = 18.02/17.72;
    if     strcmp(gender{gg}, 'male')
        k_AGT   = 577.04;
        c_ACE   = 0.88492;
        c_Chym  = 0.09315;
        c_NEP   = 0.038189;
        c_ACE2  = 0.0078009;
        c_IIIV  = 0.25056;
        c_AT1R  = 0.17008;
        c_AT2R  = 0.065667;
        AT1R_eq = 13.99;
        AT2R_eq = 5.0854;
    elseif strcmp(gender{gg}, 'female')
        k_AGT   = 610.39;
        c_ACE   = 1.4079;
        c_Chym  = 0.1482;
        c_NEP   = 0.060759;
        c_ACE2  = 0.0037603;
        c_IIIV  = 0.038644;
        c_AT1R  = 0.027089;
        c_AT2R  = 0.038699;
        AT1R_eq = 3.78;
        AT2R_eq = 5.0854;
    end
elseif strcmp(species{ss}, 'rat')
    if     strcmp(gender{gg}, 'male')
        X_PRCPRA = 135.59/17.312;
        k_AGT    = 801.02;
        c_ACE    = 0.096833;
        c_Chym   = 0.010833;
        c_NEP    = 0.012667;
        c_ACE2   = 0.0026667;
        c_IIIV   = 0.29800;
        c_AT1R   = 0.19700;
        c_AT2R   = 0.065667;
        AT1R_eq  = 20.46;
        AT2R_eq  = 6.82;
    elseif strcmp(gender{gg}, 'female')
        X_PRCPRA = 114.22/17.312;
        k_AGT    = 779.63;
        c_ACE    = 0.11600;
        c_Chym   = 0.012833;
        c_NEP    = 0.0076667;
        c_ACE2   = 0.00043333;
        c_IIIV   = 0.29800;
        c_AT1R   = 0.19700;
        c_AT2R   = 0.065667;
        AT1R_eq  = 20.46;
        AT2R_eq  = 6.82;
    end
end

% Ang II infustion rate fmol / ml min
k_AngII = 3000;
% ACEi blocking percentage
k_ACEi  = 000;
% ARB  blocking percentage
k_ARB   = 000;
% Constant renin secretion rate
R_sec = 1;

pars = [X_PRCPRA; h_renin; h_AGT; h_AngI; h_AngII; h_Ang17; h_AngIV; ...
        h_AT1R; h_AT2R; k_AGT; c_ACE; c_Chym; c_NEP; c_ACE2; c_IIIV; ...
        c_AT1R; c_AT2R; AT1R_eq; AT2R_eq; k_AngII;  k_ACEi; k_ARB; gen; R_sec];

%% Variables

% Order
% x = [R_sec; PRC; AGT; AngI; AngII; AT1R; AT2R; Ang17; AngIV];

% Initial guess for the variables; Arbitrary value for time to input;
% x0 = random('unif', 1,100, num_vars-1,1); x0 = [1;x0]; t = 0;
x0 = random('unif', 1,100, num_vars,1); t = 0;

%% Find steady state solution

options = optimset(); %options = optimset('MaxFunEvals',900+100000, 'MaxIter',400+10000);
[ss_data, residual, ...
 exitflag, output] = fsolve(@(x) RAS(t,x,pars), x0, options);

% Check for solver convergence.
if exitflag == 0
    disp('Solver did not converge.')
    disp(output)
end

% Check for imaginary solution.
if not (isreal(ss_data))
    disp('Imaginary number returned.')
end

% Set any values that are within machine precision of 0 equal to 0.
for i = 1:length(SS_DATA(:,gg))
    if abs(ss_data(i)) < eps*100
        SS_DATA(i) = 0;
    end
end

SS_DATA(:,gg)  = ss_data;
RESIDUAL(:,gg) = residual;
EXITFLAG(gg)   = exitflag;
OUTPUT{gg}     = output;

end % gender

% save_data_name = sprintf('%s_ss_data.mat', species{ss});
% % save_data_name = strcat('Data/', save_data_name);
% save(save_data_name, 'SS_DATA', 'RESIDUAL', 'EXITFLAG', 'OUTPUT')

save_data_name = sprintf('%s_ss_AngII_inf_data.mat', species{ss});
% save_data_name = strcat('Data/', save_data_name);
save(save_data_name, 'SS_DATA', 'RESIDUAL', 'EXITFLAG', 'OUTPUT')

end % species

end






























