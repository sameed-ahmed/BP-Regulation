% RAS Submodel

function run_sim

close all

species = {'human', 'rat'   };
sex     = {'male' , 'female'};
X       = cell(1,2);
T       = cell(1,2);

load('rat_ss_data.mat', 'SS_DATA')

for spe_ind = 2:2 % species
for sex_ind = 1:2 % gender

%% Parameters

if     strcmp(sex{sex_ind}, 'male')
    gen = 1;
elseif strcmp(sex{sex_ind}, 'female')
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
if     strcmp(species{spe_ind}, 'human')
    X_PRCPRA = 18.02/17.72;
    if     strcmp(sex{sex_ind}, 'male')
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
    elseif strcmp(sex{sex_ind}, 'female')
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
elseif strcmp(species{spe_ind}, 'rat')
    if     strcmp(sex{sex_ind}, 'male')
        X_PRCPRA = 135.59/17.312;
        k_AGT   = 801.02;
        c_ACE   = 0.096833;
        c_Chym  = 0.010833;
        c_NEP   = 0.012667;
        c_ACE2  = 0.0026667;
        c_IIIV  = 0.29800;
        c_AT1R  = 0.19700;
        c_AT2R  = 0.065667;
        AT1R_eq = 20.46;
        AT2R_eq = 6.82;
    elseif strcmp(sex{sex_ind}, 'female')
        X_PRCPRA = 114.22/17.312;
        k_AGT   = 779.63;
        c_ACE   = 0.11600;
        c_Chym  = 0.012833;
        c_NEP   = 0.0076667;
        c_ACE2  = 0.00043333;
        c_IIIV  = 0.29800;
        c_AT1R  = 0.19700;
        c_AT2R  = 0.065667;
        AT1R_eq = 20.46;
        AT2R_eq = 6.82;
    end
end

% Ang II infusion rate fmol / ml min
k_AngII = 000;
% ACEi blocking percentage
% k_ACEi  = 0.95;
% whole_title = sprintf('ACEi %s%%',num2str(k_ACEi*100));
k_ACEi  = 0.0;
% ARB  blocking percentage
k_ARB   = 0.86;
whole_title = sprintf('ARB %s%%',num2str(k_ARB*100));
% k_ARB   = 0.0;

pars = [X_PRCPRA; h_renin; h_AGT; h_AngI; h_AngII; h_Ang17; h_AngIV; ...
        h_AT1R; h_AT2R; k_AGT; c_ACE; c_Chym; c_NEP; c_ACE2; c_IIIV; ...
        c_AT1R; c_AT2R; AT1R_eq; AT2R_eq; k_AngII; k_ACEi; k_ARB; gen];

%% Variables

names = {'$R_{sec}$'; '$PRC$'; '$AGT$'; '$AngI$'; '$AngII$'; '$AT1R$'; ...
         '$AT2R$'; '$Ang(1-7)$'; '$AngIV$'};

% Initialization value
R_sec = 1; PRA = 1; PRC = 1; AGT = 1; 
AngI = 1; AngII = 1; 
AT1R = 10; AT2R = 1; 
Ang17 = 1; AngIV = 1; 

% x0 = [R_sec; PRC; AGT; AngI; AngII; AT1R; AT2R; Ang17; AngIV];
x0 = SS_DATA(:,sex_ind); x0 = [R_sec; x0];

%% Solve ode

% Initial time; Final time;
t0 = 0; tf = 200;
% Points per minute; Number of points;
ppm = 1; N = (tf-t0)*ppm+1;
% Time vector;
tspan = [t0, tf]; %linspace(t0,tf,N);

[t,x] = ode15s(@RAS,tspan,x0,[],pars);
T{sex_ind} = t;
X{sex_ind} = x;

end % gender

%% Plot results

% Retrieve male and female.
t1 = T{1}; t2 = T{2};
x1 = X{1}; x2 = X{2};

f = gobjects(1,1);
f(1) = figure;
s = gobjects(1,9);
% Loop through each subplot within a set of subplots.
for j = 1:9
    s(j) = subplot(3,3,j);
    plot(s(j), t1,x1(:,j), t2,x2(:,j));

    xlabel('Time (min)')
%     legend('Male', 'Female')
    title(names(j), 'Interpreter','latex', 'FontSize',15)
end
sgtitle(whole_title, 'FontSize',16)

% if     strcmp(species{ss}, 'human')
%     whole_title = suptitle('Human');
%     set(whole_title,'FontSize',20)
%     savefig(f, 'all_vars_const_Rsec_human.fig')
% elseif strcmp(species{ss}, 'rat')
%     whole_title = suptitle('Rat');
%     set(whole_title,'FontSize',20)
%     savefig(f, 'all_vars_const_Rsec_rat.fig')
% end

end % species

end




































