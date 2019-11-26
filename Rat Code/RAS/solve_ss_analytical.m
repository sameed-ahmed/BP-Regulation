% This script finds the steady state solution to the system of differential
% equations that model the RAS pathway. This is done by setting the 
% derivative = 0, setting up the symbolic variables and parameters, and 
% using MATLAB's 'solve' function to solve this system of equations.

% The analytic expression is saved. This can be loaded into another
% function. A given set of parameters can then be plugged into this
% analytic expression for the numerical value of the steady state solution.

function ss_val_matr = solve_ss_analytical

num_vars = 8;
num_pars = 21;

eqns = sym('eqns', [num_vars,1]);
vars = sym('vars', [num_vars,1]); 
pars = sym('pars', [num_pars,1]); 

%% Retreive parameters by name.

% pars = [X_PRCPRA; h_renin; h_AGT; h_AngI; h_AngII; h_Ang17; h_AngIV; ...
%         h_AT1R; h_AT2R; k_AGT; c_ACE; c_Chym; c_NEP; c_ACE2; c_IIIV; ...
%         c_AT1R; c_AT2R; AT1R_eq; AT2R_eq; k_AngII; R_sec];
pars(1 ) = 'X_PRCPRA'; pars(2 ) = 'h_renin'; pars(3 ) = 'h_AGT';      
pars(4 ) = 'h_AngI';   pars(5 ) = 'h_AngII'; pars(6 ) = 'h_Ang17';
pars(7 ) = 'h_AngIV';  pars(8 ) = 'h_AT1R';  pars(9 ) = 'h_AT2R'; 
pars(10) = 'k_AGT';    pars(11) = 'c_ACE';   pars(12) = 'c_Chym'; 
pars(13) = 'c_NEP';    pars(14) = 'c_ACE2';  pars(15) = 'c_IIIV';         
pars(16) = 'c_AT1R';   pars(17) = 'c_AT2R';  pars(18) = 'AT1R_eq'; 
pars(19) = 'AT2R_eq';  pars(20) = 'k_AngII'; pars(21) = 'R_Sec';

X_PRCPRA = pars(1 );
h_renin  = pars(2 );
h_AGT    = pars(3 );
h_AngI   = pars(4 );
h_AngII  = pars(5 );
h_Ang17  = pars(6 );
h_AngIV  = pars(7 );
h_AT1R   = pars(8 );
h_AT2R   = pars(9 );
k_AGT    = pars(10);
c_ACE    = pars(11);
c_Chym   = pars(12);
c_NEP    = pars(13);
c_ACE2   = pars(14);
c_IIIV   = pars(15);
c_AT1R   = pars(16);
c_AT2R   = pars(17);
AT1R_eq  = pars(18);
AT2R_eq  = pars(19);
k_AngII  = pars(20);
R_sec    = pars(21);

%% Retrieve variables by name.

% vars = [PRC; AGT; AngI; AngII; AT1R; AT2R; Ang17; AngIV];
% R_sec = 1

vars(1 ) = 'PRC';   vars(2 ) = 'AGT';   vars(3 ) = 'AngI'; 
vars(4 ) = 'AngII'; vars(5 ) = 'AT1R';  vars(6 ) = 'AT2R'; 
vars(7 ) = 'Ang17'; vars(8 ) = 'AngIV'; 

PRC   = vars(1);
AGT   = vars(2);
AngI  = vars(3);
AngII = vars(4);
AT1R  = vars(5);
AT2R  = vars(6);
Ang17 = vars(7);
AngIV = vars(8);

%% Equations

% PRC
eqns(1) = R_sec - log(2)/h_renin * PRC;
% AGT
eqns(2) = k_AGT - X_PRCPRA * PRC - log(2)/h_AGT * AGT;
% AngI
eqns(3) = X_PRCPRA * PRC - (c_ACE + c_Chym + c_NEP) * AngI - log(2)/h_AngI * AngI;
% AngII
eqns(4) = k_AngII + (c_ACE + c_Chym) * AngI - (c_ACE2 + c_IIIV + c_AT1R + c_AT2R) * AngII - log(2)/h_AngII * AngII;
% AT1R
eqns(5) = c_AT1R * AngII - log(2)/h_AT1R * AT1R;
% AT2R
eqns(6) = c_AT2R * AngII - log(2)/h_AT2R * AT2R;
% Ang17
eqns(7) = c_NEP * AngI + c_ACE2 * AngII - log(2)/h_Ang17 * Ang17;
% AngIV
eqns(8) = c_IIIV * AngII - log(2)/h_AngIV * AngIV;

%% Find analytic steady state solution.

tic
ss_sol = solve(eqns, vars);
ss_sol = struct2array(ss_sol)';
CPUtime = toc;

save('ss_sol.mat', 'ss_sol', 'CPUtime')

%% Substitute parameters values into analytic expression.

species     = {'human', 'rat'};
gender      = {'male' , 'female'};
ss_val_anal = cell(2, 2); % (species, gender)
ss_val_matr = cell(2, 2); % (species, gender)

for ss = 1:2 % species
for gg = 1:2 % gender

% RAS
h_renin  = 12;      % min
h_AGT    = 10*60;   % min
h_AngI   = 0.5;     % min
h_AngII  = 0.66;    % min
h_Ang17  = 30;      % min
h_AngIV  = 0.5;     % min
h_AT1R   = 12;      % min
h_AT2R   = 12;      % min
R_sec    = 1;
% R_sec    = 1.01902709680815; % male
% R_sec    = 0.998789033750504; % female

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

pars_num = [X_PRCPRA; h_renin; h_AGT; h_AngI; h_AngII; h_Ang17; ...
            h_AngIV; h_AT1R; h_AT2R; k_AGT; c_ACE; c_Chym; c_NEP; ...
            c_ACE2; c_IIIV; c_AT1R; c_AT2R; AT1R_eq; AT2R_eq; k_AngII; R_sec];

if     strcmp(species{ss}, 'human')
    if     strcmp(gender{gg}, 'male')
        key = 'human male';
    elseif strcmp(gender{gg}, 'female')
        key = 'human female';
    end
elseif strcmp(species{ss}, 'rat')
    if     strcmp(gender{gg}, 'male')
        key = 'rat male';
    elseif strcmp(gender{gg}, 'female')
        key = 'rat female';
    end
end

ss_val_anal{ss, gg} = double(subs(ss_sol, pars, pars_num));
ss_val_anal{ss, gg} = {ss_val_anal{ss, gg}, key};

%% Using x = A^(-1)b

b = zeros(8,1);
A = zeros(8,8);

b(1) = -R_sec; b(2) = -k_AGT; b(4) = -k_AngII;

A(1,1) = -log(2)/h_renin; 
A(2,1) = -X_PRCPRA; A(2,2) = -log(2)/h_AGT; 
A(3,1) = X_PRCPRA; A(3,3) = -(c_ACE + c_Chym + c_NEP + log(2)/h_AngI); 
A(4,3) = c_ACE + c_Chym; A(4,4) = -(c_ACE2 + c_IIIV + c_AT1R + c_AT2R + log(2)/h_AngII); 
A(5,4) = c_AT1R; A(5,5) = -log(2)/h_AT1R; 
A(6,4) = c_AT2R; A(6,6) = -log(2)/h_AT2R; 
A(7,3) = c_NEP; A(7,4) = c_ACE2; A(7,7) = -log(2)/h_Ang17; 
A(8,4) = c_IIIV; A(8,8) = -log(2)/h_AngIV;

% vars = [PRC; AGT; AngI; AngII; AT1R; AT2R; Ang17; AngIV];
ss_val_matr{ss, gg} = A\b;
ss_val_matr{ss, gg} = {ss_val_matr{ss, gg}, key};


end % gender
end % species

% save('ss_val.mat', 'ss_val_anal', 'ss_val_matr')

% err = cell(2, 2);
% for i = 1:2
%     for j = 1:2
%         err{i,j} = max(abs(ss_val_anal{i,j}{1} - ss_val_matr{i,j}{1}));
%     end
% end

end


































