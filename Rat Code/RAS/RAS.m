function dx = RAS(t,x,pars)

%% Retrieve parameters by name.

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

gen      = pars(21);
if     gen == 1
    gender = 'male';
elseif gen == 0
    gender = 'female';
end

%% Retrieve variables by name.

% R_sec = x(1);  
% PRC   = x(2); 
% AGT   = x(3); 
% AngI  = x(4); 
% AngII = x(5); 
% AT1R  = x(6); 
% AT2R  = x(7); 
% Ang17 = x(8); 
% AngIV = x(9); 

R_sec = pars(end);
PRC   = x(1); 
AGT   = x(2); 
AngI  = x(3); 
AngII = x(4); 
AT1R  = x(5); 
AT2R  = x(6); 
Ang17 = x(7); 
AngIV = x(8); 

%% Differential equation system dxdt = A(t,x).

dx = zeros(length(x),1);

% % R_sec
% dx(1) = 0;
% % dx(1) = 10^(0.0102) * (-0.95) * (AT1R / AT1R_eq)^(-0.95-1) * ( c_AT1R * AngII - log(2)/h_AT1R * AT1R ) / AT1R_eq;
% % PRC
% dx(2) = R_sec - log(2)/h_renin * PRC;
% % AGT
% dx(3) = k_AGT - X_PRCPRA * PRC - log(2)/h_AGT * AGT;
% % AngI
% dx(4) = X_PRCPRA * PRC - (c_ACE + c_Chym + c_NEP) * AngI - log(2)/h_AngI * AngI;
% % AngII
% dx(5) = (c_ACE + c_Chym) * AngI - (c_ACE2 + c_IIIV + c_AT1R + c_AT2R) * AngII - log(2)/h_AngII * AngII;
% % AT1R
% dx(6) = c_AT1R * AngII - log(2)/h_AT1R * AT1R;
% % AT2R
% dx(7) = c_AT2R * AngII - log(2)/h_AT2R * AT2R;
% % Ang17
% dx(8) = c_NEP * AngI + c_ACE2 * AngII - log(2)/h_Ang17 * Ang17;
% % AngIV
% dx(9) = c_IIIV * AngII - log(2)/h_AngIV * AngIV;

% PRC
dx(1) = R_sec - log(2)/h_renin * PRC;
% AGT
dx(2) = k_AGT - X_PRCPRA * PRC - log(2)/h_AGT * AGT;
% AngI
dx(3) = X_PRCPRA * PRC - (c_ACE + c_Chym + c_NEP) * AngI - log(2)/h_AngI * AngI;
% AngII
dx(4) = k_AngII + (c_ACE + c_Chym) * AngI - (c_ACE2 + c_IIIV + c_AT1R + c_AT2R) * AngII - log(2)/h_AngII * AngII;
% AT1R
dx(5) = c_AT1R * AngII - log(2)/h_AT1R * AT1R;
% AT2R
dx(6) = c_AT2R * AngII - log(2)/h_AT2R * AT2R;
% Ang17
dx(7) = c_NEP * AngI + c_ACE2 * AngII - log(2)/h_Ang17 * Ang17;
% AngIV
dx(8) = c_IIIV * AngII - log(2)/h_AngIV * AngIV;

% b = zeros(8,1);
% A = zeros(8,8);
% 
% b(1) = -R_sec; b(2) = -k_AGT;
% 
% A(1,1) = -log(2)/h_renin; 
% A(2,1) = -X_PRCPRA; A(2,2) = -log(2)/h_AGT; 
% A(3,1) = X_PRCPRA; A(3,3) = -(c_ACE + c_Chym + c_NEP + log(2)/h_AngI); 
% A(4,3) = c_ACE + c_Chym; A(4,4) = -(c_ACE2 + c_IIIV + c_AT1R + c_AT2R + log(2)/h_AngII); 
% A(5,4) = c_AT1R; A(5,5) = -log(2)/h_AT1R; 
% A(6,4) = c_AT2R; A(6,6) = -log(2)/h_AT2R; 
% A(7,3) = c_NEP; A(7,4) = c_ACE2; A(7,7) = -log(2)/h_Ang17; 
% A(8,4) = c_IIIV; A(8,8) = -log(2)/h_AngIV;

% dx = A*x - b;

end

































