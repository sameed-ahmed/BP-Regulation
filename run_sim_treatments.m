% This simulates the blood pressure regulation model blood_press_reg.m.
% 
% Parameters are given by:
% "Long-Term Mathematical Model Involving Renal Sympathetic Nerve Activity,
% Arterial Pressure, and Sodium Excretion" - 2005 - Karaaslan, et. al.
% "Sex-specific Long-term Blood Pressure Regulation: Modeling and Analysis"
% - 2018 - Leete, Layton.
% 
% Steady state data is calculated by solve_ss_numerical.m.

function run_sim_treatments(treatments_to_run)
close all

%AA = 8;
human = 1;% 1 for human, 0 for rat
species = {'rat','human'};
gender   = {'male', 'female'};
disp([species{human+1}])

%SSDATA   = zeros(82,2);
%residual = zeros(82,2);
X        = cell(1,2);
T        = cell(1,2);
%columns are ACEI,diuretic,NSAID
models = [0,0,0;... %normal dose
    1,0,0;...
    0,1,0;...
    0,0,1;...
    1,1,0;...
    0,1,1;...
    1,0,1;...
    1,1,1;...
    2,0,0;...%High ACEI 9
    2,1,0;...
    2,0,1;...
    2,1,1;...
    0,2,0;...%high furosemide 13
    1,2,0;...
    0,2,1;...
    1,2,1;...
    0,0,2;...%high NSAID 17
    0,1,2;...
    1,0,2;...
    1,1,2;...
    2,2,0;...% All high combos 21
    0,2,2;...
    2,0,2;...
    2,2,2];
sex_to_run = 1:2;
hyp = [2.5];
for gg = 1%sex_to_run %sex
    for mm= treatments_to_run 
    for AA = 1:length(hyp)
         disp([gg,mm,hyp(AA)])
         
pars = get_params(species{human+1},gender{gg},hyp(AA));

ACEI = models(mm,1);
diuretic = models(mm,2);
NSAID = models(mm,3);

kappa_ACEI = 0;
kappa_d = 0;
kappa_d_tgf = 0;
kappa_d_renin = 0;

   if ACEI == 1
        kappa_ACEI = 0.76;
   elseif ACEI > 1
       kappa_ACEI = 1;
   end
   if diuretic == 1
       kappa_d = 0.15*diuretic;
       kappa_d_tgf = 0.4*diuretic;
       kappa_d_renin = 0.4*diuretic;
   elseif diuretic == 2
       kappa_d = 0.15*diuretic;
       kappa_d_tgf = 0.5;
       kappa_d_renin = 0.5;
   end  


drugs = [kappa_ACEI,kappa_d,kappa_d_tgf,kappa_d_renin,NSAID];

%disp(drugs)
%% Solve DAE

% Initial value
% This initial condition is the steady state data value taken from
% experiments (CITE). Therefore, the initial condition of the derivative is
% 0.


% Load data for steady state initial value. 
name = sprintf('Data_testsexdiff/%s_%s_ss_%s_%s_%s_rsna%s_femalegammaRSNA.mat',...
    species{human+1},gender{gg},num2str(0),num2str(0),num2str(0),num2str(hyp(AA)));
 if (mm == 8) || (mm == 12) || (mm == 16)|| (mm == 20) || (mm==24) %triple treatment
     name = sprintf('Data_testsexdiff/%s_%s_ss_%s_%s_%s_rsna%s_femalegammaRSNA.mat',...
         species{human+1},gender{gg},num2str(0),num2str(0),num2str(models(mm,3)),num2str(hyp(AA)));
 elseif (mm ==6) || (mm==15) || (mm==18) || (mm==22) % F+N
      name = sprintf('Data_testsexdiff/%s_%s_ss_%s_%s_%s_rsna%s_femalegammaRSNA.mat', ...
          species{human+1},gender{gg},num2str(0),num2str(0),num2str(models(mm,3)),num2str(hyp(AA)));
 elseif (mm ==7) || (mm==11) || (mm==19) || (mm==23) %A+N
      name = sprintf('Data_testsexdiff/%s_%s_ss_%s_%s_%s_rsna%s_femalegammaRSNA.mat',...
          species{human+1},gender{gg},num2str(0),num2str(0),num2str(models(mm,3)),num2str(hyp(AA)));
 elseif (mm==4) || (mm==17)
     name = sprintf('Data_testsexdiff/%s_%s_ss_%s_%s_%s_rsna%s_femalegammaRSNA.mat',...
         species{human+1},gender{gg},num2str(0),num2str(0),num2str(models(mm,3)),num2str(hyp(AA)));
 end
   load(name,'SSdata');
   SS_data_IG = SSdata;
   
   if length(SS_data_IG) < 84
    SS_data_IG(84) = 0.126;
   end
   
names  = {'$rsna$'; '$\alpha_{map}$'; '$\alpha_{rap}$'; '$R_{r}$'; ...
          '$\beta_{rsna}$'; '$\Phi_{rb}$'; '$\Phi_{gfilt}$'; '$P_{f}$'; ...
          '$P_{gh}$'; '$\Sigma_{tgf}$'; '$\Phi_{filsod}$'; ...
          '$\Phi_{pt-sodreab}$'; '$\eta_{pt-sodreab}$'; ...
          '$\gamma_{filsod}$'; '$\gamma_{at}$'; '$\gamma_{rsna}$'; ...
          '$\Phi_{md-sod}$'; '$\Phi_{dt-sodreab}$'; ...
          '$\eta_{dt-sodreab}$'; '$\psi_{al}$'; '$\Phi_{dt-sod}$'; ...
          '$\Phi_{cd-sodreab}$'; '$\eta_{cd-sodreab}$'; ...
          '$\lambda_{dt}$'; '$\lambda_{anp}$'; '$\Phi_{u-sod}$'; ...
          '$\Phi_{win}$'; '$V_{ecf}$'; '$V_{b}$'; '$P_{mf}$'; ...
          '$\Phi_{vr}$'; '$\Phi_{co}$'; '$P_{ra}$'; '$vas$'; ...
          '$vas_{f}$'; '$vas_{d}$'; '$R_{a}$'; '$R_{ba}$'; '$R_{vr}$'; ...
          '$R_{tp}$'; '$P_{ma}$'; '$\epsilon_{aum}$'; '$a_{auto}$'; ...
          '$a_{chemo}$'; '$a_{baro}$'; '$C_{adh}$'; '$N_{adh}$'; ...
          '$N_{adhs}$'; '$\delta_{ra}$'; '$\Phi_{t-wreab}$'; ...
          '$\mu_{al}$'; '$\mu_{adh}$'; '$\mu_{Na}$'; '$\Phi_{u}$'; '$M_{sod}$'; ...
          '$C_{sod}$'; '$\nu_{md-sod}$'; '$\nu_{rsna}$'; '$C_{al}$'; ...
          '$N_{al}$'; '$N_{als}$'; '$\xi_{k/sod}$'; '$\xi_{map}$'; ...
          '$\xi_{at}$'; '$\hat{C}_{anp}$'; '$AGT$'; '$\nu_{AT1}$'; ...
          '$R_{sec}$'; '$PRC$'; '$PRA$'; '$Ang I$'; '$Ang II$'; ...
          '$Ang II_{AT1R-bound}$'; '$Ang II_{AT2R-bound}$'; ...
          '$Ang (1-7)$'; '$Ang IV$'; '$R_{aa}$'; '$R_{ea}$'; ...
          '$\Sigma_{myo}$'; '$\Psi_{AT1R-AA}$'; '$\Psi_{AT1R-EA}$'; ...
          '$\Psi_{AT2R-AA}$'; '$\Psi_{AT2R-EA}$';'$\Phi_{sod-in}$'};

% Initial condition for the variables and their derivatives. 
% System is initially at steady state, so the derivative is 0.
x0 = SS_data_IG; x_p0 = zeros(84,1);

% Time at which to keep steady state, change a parameter, etc.
tchange = 100;

% Initial time (min); Final time (min);\
num_days =2;
t0 = 0*1440; tf = tchange + num_days*1440;
% Time vector;
tspan = t0:10:tf;%[t0, tf];

% Solve dae
[t,x] = ode15i(@(t,x,x_p) blood_press_reg_sim(t,x,x_p,pars,tchange,drugs), tspan, x0, x_p0);
T{gg,mm,AA} = t;
X{gg,mm,AA} = x;

%%Save Data
if     strcmp(gender{gg}, 'male')
    filename = '%s%s_male_%s_%s_%s_%s_femalegammaRSNA.mat';
elseif strcmp(gender{gg}, 'female')
    filename = '%s%s_female_%s_%s_%s_%s_femalegammaRSNA.mat';
end

if not (isreal(x))
    disp('Imaginary number returned.')
end
filename = sprintf( filename,'Data_testsexdiff/',species{human+1},num2str(ACEI),num2str(diuretic),num2str(NSAID),num2str(hyp(AA)));

save(filename, 't','x');

    end
end
end
































