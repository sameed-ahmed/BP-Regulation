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
plot_it = 1;
plot_these = [7,41,54,84];

%AA = 8;
human = 1;% 1 for human, 0 for rat
species = {'rat','human'};
gender   = {'male', 'female'};
disp([species{human+1}])

%SSDATA   = zeros(82,2);
%residual = zeros(82,2);
X        = cell(1,2);
T        = cell(1,2);
models = [0,0,0;...
    1,0,0;...
    0,1,0;...
    0,0,1;...
    1,1,0;...
    0,1,1;...
    1,0,1;...
    1,1,1];%control;ACEI;diuretic;NSAID

sex_to_run = 1:2;
hyp = [1,2.5];
for gg = sex_to_run %sex
    for mm= treatments_to_run 
    for AA = 1:length(hyp)
         disp([gg,mm,hyp(AA)])
         
%% Parameters

if     strcmp(gender{gg}, 'male')
    gen = 1;
elseif strcmp(gender{gg}, 'female')
    gen = 0;
end

% Scaling factor
% Rat flow = Human flow x SF
if human
    SF = 1.0;
else
    SF = 4.5*10^(-3);
end

N_rsna      = 1*hyp(AA);
R_aass      = 31.67 / SF;   % mmHg min / l
R_eass      = 51.66 / SF;   % mmHg min / l
P_B         = 18;           % mmHg
P_go        = 28;           % mmHg
C_gcf       = 0.00781 * SF;
eta_etapt   = 0.8; 
eta_epsdt   = 0.5; 
eta_etacd   = 0.93; 
K_vd        = 0.00001;
K_bar       = 16.6 / SF;    % mmHg min / l
R_bv        = 3.4 / SF;     % mmHg min / l
T_adh       = 6;            % min
Phi_sodin   = 0.126 * SF;   % mEq / min
%if     strcmp(gender{gg}, 'male')
%   C_K      = 6;            % mEq / l 
%elseif strcmp(gender{gg}, 'female')
   C_K      = 5;            % mEq / l 
%end
T_al        = 30;           % min LISTED AS 30 IN TABLE %listed as 60 in text will only change dN_al
N_rs        = 1;            % ng / ml / min

% RAS
h_renin     = 12;      % min
h_AGT       = 10*60;   % min
h_AngI      = 0.5;     % min
h_AngII     = 0.66;    % min
h_Ang17     = 30;      % min
h_AngIV     = 0.5;     % min
h_AT1R      = 12;      % min
h_AT2R      = 12;      % min

c_GPautoreg = 5;
P_ghnom     = 62;      % mmHg

% Male and female different parameters for RAS
if human
    X_PRCPRA = 61/60.0; %fmol/min/pg
    if strcmp(gender{gg}, 'male')     
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
        ALD_eq = 85;
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
        ALD_eq = 69.1775;


    end
    pt_sod_reab_EQ = 13.909;%15.285;%;14.4;%
    dt_sod_reab_EQ = 1.5859;%2.0714;%1.8;
    cd_sod_reab_EQ = 1.6909;%1.6934;%1.674;
else %rat
    if     strcmp(gender{gg}, 'male')
        C_K      = 6;            % mEq / l 
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
        ALD_eq = 395.3;
    elseif strcmp(gender{gg}, 'female')
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
        ALD_eq = 379.4;
    end
    pt_sod_reab_EQ = 0.068294039337572;
    dt_sod_reab_EQ = 0.008094477703862;
    cd_sod_reab_EQ = 0.007616239742696;
end
amount = 0.08;
slope = 8.5;
if     strcmp(gender{gg}, 'male')
    if     AA>1
        baseline = 0.0183;
    else
        baseline = 0.0199;
    end
else
    if     AA>1
        baseline = 0.01843;
    else
        baseline = 0.0222;
    end
end

pars = [N_rsna; R_aass; R_eass; P_B; P_go; C_gcf; eta_etapt; eta_epsdt; ...
        eta_etacd; K_vd; K_bar; R_bv; T_adh; Phi_sodin; C_K; T_al; ...
        N_rs; X_PRCPRA; h_renin; h_AGT; h_AngI; h_AngII; h_Ang17; ...
        h_AngIV; h_AT1R; h_AT2R; c_GPautoreg; P_ghnom; k_AGT; c_ACE; ...
        c_Chym; c_NEP; c_ACE2; c_IIIV; c_AT1R; c_AT2R; AT1R_eq; ...
        AT2R_eq; ALD_eq; gen;...
        pt_sod_reab_EQ; dt_sod_reab_EQ; cd_sod_reab_EQ;amount; slope; baseline; human;SF];

ACEI = models(mm,1);
diuretic = models(mm,2);
NSAID = models(mm,3);

kappa_ACEI = 0;
kappa_d = 0;
kappa_d_tgf = 0;
kappa_d_renin = 0;
if human
   if ACEI
        kappa_ACEI = 0.55;
   end
   if diuretic
       kappa_d = 0.15;
       kappa_d_tgf = 0.4;
       kappa_d_renin = 0.4;
   end  

else
    if ACEI
        kappa_ACEI = 0.646;
    end
    if diuretic
        kappa_d = 0.08;
        kappa_d_tgf = 0.04;
        kappa_d_renin = 0.23;
    end
end

drugs = [kappa_ACEI,kappa_d,kappa_d_tgf,kappa_d_renin,NSAID];


%% Solve DAE

% Initial value
% This initial condition is the steady state data value taken from
% experiments (CITE). Therefore, the initial condition of the derivative is
% 0.

% Add directory containing data.
%mypath = pwd;
%mypath = strcat(mypath, '/Data');
%addpath(genpath(mypath))

% Load data for steady state initial value. 
name = sprintf('Data/%s_%s_ss_%s_%s_%s_rsna%s_fixedparams.mat', species{human+1},gender{gg},num2str(0),num2str(0),num2str(0),num2str(hyp(AA)));
 if mm == 8
     name = sprintf('Data/%s_%s_ss_%s_%s_%s_rsna%s_fixedparams.mat', species{human+1},gender{gg},num2str(1),num2str(1),num2str(0),num2str(hyp(AA)));
 elseif mm ==6
      name = sprintf('Data/%s_%s_ss_%s_%s_%s_rsna%s_fixedparams.mat', species{human+1},gender{gg},num2str(0),num2str(1),num2str(0),num2str(hyp(AA)));
 elseif mm ==7
      name = sprintf('Data/%s_%s_ss_%s_%s_%s_rsna%s_fixedparams.mat', species{human+1},gender{gg},num2str(1),num2str(0),num2str(0),num2str(hyp(AA)));
 elseif (mm==4)
     name = sprintf('Data/%s_%s_ss_%s_%s_%s_rsna%s_fixedparams.mat', species{human+1},gender{gg},num2str(0),num2str(0),num2str(1),num2str(hyp(AA)));
 
 end
   load(name,'SSdata');
   SS_data_IG = SSdata;
   SS_data_IG(84) = 0.126;

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
          '$\Psi_{AT2R-AA}$'; '$\Psi_{AT2R-EA}$'};

% Initial condition for the variables and their derivatives. 
% System is initially at steady state, so the derivative is 0.
x0 = SS_data_IG; x_p0 = zeros(84,1);

% Time at which to keep steady state, change a parameter, etc.
tchange = 100;

% Initial time (min); Final time (min);\
num_days = 3;
t0 = 0*1440; tf = tchange + num_days*1440;
% Time vector;
tspan = [t0, tf];

% Solve dae
[t,x] = ode15i(@(t,x,x_p) blood_press_reg_sim(t,x,x_p,pars,tchange,drugs), tspan, x0, x_p0);
T{gg,mm,AA} = t;
X{gg,mm,AA} = x;

%%Save Data
if     strcmp(gender{gg}, 'male')
    filename = '%smale_%s_%s_%s_%s_fixedparams.mat';
elseif strcmp(gender{gg}, 'female')
    filename = '%sfemale_%s_%s_%s_%s_fixedparams.mat';
end

if not (isreal(x))
    disp('Imaginary number returned.')
    plot_it = 0;
end
filename = sprintf( filename,'Data/',num2str(ACEI),num2str(diuretic),num2str(NSAID),num2str(hyp(AA)));

save(filename, 't','x');

    end
end
end

% Plot
if 0
% Retrieve data
% t1 = T{1}; t2 = T{2}; t3 = T{3}; t4 = T{4};%t5 = T{5};
% x1 = X{1}; x2 = X{2}; x3 = X{3}; x4 = X{4};%x5 = X{5};

linestyle = {'-','--','--','--',':'};
f = gobjects(6,1);
s = gobjects(6,15);
% Loop through each set of subplots.
for i = 1:6
    f(i) = figure;
    % This is to avoid the empty plots in the last subplot set.
    if i == 6
        last_plot = 8;
    else
        last_plot = 15;
    end
    
    % Loop through each subplot within a set of subplots.
    for j = 1:last_plot
        s(i,j) = subplot(3,5,j);
        for mm=treatments_to_run
            for gg=sex_to_run%1:2
                for AA = [1,2]
                x1=X{gg,mm,AA};
                plot(s(i,j), T{gg,mm,AA},x1(:,(i-1)*15 + j), linestyle{mm});%, t2,x2(:,(i-1)*15 + j),t3,x3(:,(i-1)*15 + j),t4,x4(:,(i-1)*15 + j));%,t5,x5(:,(i-1)*15 + j));
                hold on
                end
            end
        end
        % Set the y-axis limits to be within 10% of the steady state value.
%         ss_value_max = max(x1(1,(i-1)*15 + j));%, x2(1,(i-1)*15 + j));
%         ss_value_min = min(x1(1,(i-1)*15 + j));%, x2(1,(i-1)*15 + j));
%         lower_lim = ss_value_min - 0.1*abs(ss_value_min);
%         upper_lim = ss_value_max + 0.1*abs(ss_value_max);
%         if ss_value_max == 0
%            lower_lim = -10^(-5); %-0.01; 
%            upper_lim =  10^(-5); % 0.01; 
%         end
%         ylim([lower_lim, upper_lim])
         %xlim([-1000,tf+1000])
         xlim([-10,tf+10])
        xlabel('Time (min)')
        if j==1
                lnames = {'Control','ACEI', 'Diuretic', 'NSAID','AD','DN','AN','Triple'};
            legend(lnames{treatments_to_run})
        end
        title(names((i-1)*15 + j), 'Interpreter','latex', 'FontSize',15)
    end
end
savefig(f, sprintf('Figures/%s_%s.fig',species{human+1},num2str(treatments_to_run)))

end


if 0
    ylim_labels = {'l/min','mmHg','l/min','ng/l','fmol/ml/min','mEq/min'};
    ylimits = [0.06,0.16;80,150; 0,0.008];
    linestyle = {'-',':',':',':','--','--','--',':'};
    num_columns = 3;
    f = gobjects(length(plot_these),1);
    s = gobjects(length(plot_these),length(sex_to_run)*length(hyp)*num_columns);
    
    for i = 1:length(plot_these)
        f(i) = figure;
        for j = 1:length(sex_to_run)*length(hyp)*num_columns
            s(i,j) = subplot(length(sex_to_run)*length(hyp),num_columns,j);
            if j < (num_columns +1) %male
                gg = 1;
                AA =1;
            elseif j < (2*num_columns +1)
                gg = 2;
                AA = 1;
            elseif j < (3*num_columns +1)
                gg = 1;
                AA = 2;
            else 
                gg = 2;
                AA = 2; 
            end
            if mod(j,num_columns) == 1 %single treatments
                treat_to_plot = 1:4;
            elseif mod(j,num_columns) == 2
                treat_to_plot = [1,5:7];
            elseif num_columns > 1
                treat_to_plot = [1,8];
            else
                treat_to_plot = treatments_to_run;
            end
            for mm=1:length(treat_to_plot)
                %disp([gg,treat_to_plot(mm),AA,i,plot_these(i)])
                x1=X{gg,treat_to_plot(mm),AA};  
                plot(s(i,j), T{gg,treat_to_plot(mm),AA},x1(:,plot_these(i)), linestyle{mm},'Linewidth',3);
                hold on
                if i==3
                    hline(0.0003,'k')
                end
            end
            xlim([-10,tf+10])
            %ylim(ylimits(i,:))
            set(gca,'YGrid','on','FontSize',20);
            if j == 1
                ylabel({'Male';ylim_labels{i}})
                title('Single Treatments')
                %legend({'Control','ACEI','Furosemide','NSAID'})
            elseif j ==num_columns + 1
                ylabel({'Female';ylim_labels{i}})
            elseif j==2*num_columns + 1
                ylabel({'M HTN';ylim_labels{i}})
            elseif j==3*num_columns + 1
                ylabel({'F HTN';ylim_labels{i}})
                xlabel('Time (min)')
            elseif j==2
                title('Double Treatments')
                %legend({'Control','A+F','F+N','A+N'})
            elseif j==3
                title('Triple Treatment')
                %legend({'Control','A+F+N'})
            
            elseif j > 3*num_columns + 1
                xlabel('Time (min)')
            end
        end
        
        %suptitle(names(plot_these(i)), 'Interpreter','latex', 'FontSize',15)
    end
end

if 0
    treat_to_plot = [7,41,54];
    names  = {'$rsna$'; '$\alpha_{map}$'; '$\alpha_{rap}$'; '$R_{r}$'; ...
          '$\beta_{rsna}$'; '$\Phi_{rb}$'; 'GFR'; '$P_{f}$'; ...
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
          '$\mu_{al}$'; '$\mu_{adh}$'; '$\mu_{Na}$';'$\Phi_{u}$'; '$M_{sod}$'; ...
          '$C_{sod}$'; '$\nu_{md-sod}$'; '$\nu_{rsna}$'; '$C_{al}$'; ...
          '$N_{al}$'; '$N_{als}$'; '$\xi_{k/sod}$'; '$\xi_{map}$'; ...
          '$\xi_{at}$'; '$\hat{C}_{anp}$'; '$AGT$'; '$\nu_{AT1}$'; ...
          '$R_{sec}$'; '$PRC$'; '$PRA$'; '$Ang I$'; '$Ang II$'; ...
          '$Ang II_{AT1R-bound}$'; '$Ang II_{AT2R-bound}$'; ...
          '$Ang (1-7)$'; '$Ang IV$'; '$R_{aa}$'; '$R_{ea}$'; ...
          '$\Sigma_{myo}$'; '$\Psi_{AT1R-AA}$'; '$\Psi_{AT1R-EA}$'; ...
          '$\Psi_{AT2R-AA}$'; '$\Psi_{AT2R-EA}$'};
   ylim_labels = {'l/min','mmHg','l/min'}; 
    ylimits = [0.06,0.16;80,150; 0,0.008];
   f = gobjects(length(plot_these),1);
        s = gobjects(length(plot_these),num_days);
   for i = 1:length(treat_to_plot)
        f(i) = figure;
    
        data = zeros(8,2*2,num_days);
        for gg= sex_to_run
            for AA = 1:length(hyp)
                for mm = 1:8
                    x1=X{gg,mm,AA};  
                    t = T{gg,mm,AA};
                    if t(end) > tchange
                        for dd = 1:num_days
                            mask1 = (t>tchange+(dd-1)*1440);
                            mask2 = (t < tchange + dd*1440);
                            mask = mask1&mask2;
                            day = t(mask);
                            dt = day(2:end) - day(1:end-1);
                            x_day = x1(mask,treat_to_plot(i));
                            ave = sum(dt .* x_day(2:end))/(day(end) - day(1));
                            data(mm,2*(AA-1) + gg,dd) = ave;
                        end
                    else
                        data(mm,2*(AA-1) + gg,:) = zeros(1,num_days);
                    end
                end
            end
        end
        for dd = 1:num_days
            s(i,dd) = subplot(num_days,1,dd);
            bar(data(:,:,dd))
            
            labels = {'Control','ACEI','Furosemide','NSAID','Furosemide\newline+ ACEI','Furosemide\newline+ NSAID','ACEI\newline+ NSAID','Furosemide\newline+ ACEI\newline+ NSAID'};
            if dd==1
                %title(names{treat_to_plot(i)}, 'Interpreter','latex')%, 'FontSize',24)
                legend({'Male','Female','M HTN','F HTN'});
                set(gca,'xtick',[])
            elseif dd < num_days
                set(gca,'xtick',[])
            elseif dd==num_days
                set(gca,'xticklabel',labels)%,'FontSize',24)
            end
            if i==3
                hline([0.0003,0.001])
            end
            title(sprintf('Day %s',num2str(dd)))
            ylabel(ylim_labels{i})
            xlim([0.4,8.6])
            ylim(ylimits(i,:))
            set(gca,'YGrid','on','FontSize',24);
        end
   end
end

if 1

    names  = {'$rsna$'; '$\alpha_{map}$'; '$\alpha_{rap}$'; '$R_{r}$'; ...
          '$\beta_{rsna}$'; '$\Phi_{rb}$'; 'GFR'; '$P_{f}$'; ...
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
          '$\mu_{al}$'; '$\mu_{adh}$'; '$\mu_{Na}$';'$\Phi_{u}$'; '$M_{sod}$'; ...
          '$C_{sod}$'; '$\nu_{md-sod}$'; '$\nu_{rsna}$'; '$C_{al}$'; ...
          '$N_{al}$'; '$N_{als}$'; '$\xi_{k/sod}$'; '$\xi_{map}$'; ...
          '$\xi_{at}$'; '$\hat{C}_{anp}$'; '$AGT$'; '$\nu_{AT1}$'; ...
          '$R_{sec}$'; '$PRC$'; '$PRA$'; '$Ang I$'; '$Ang II$'; ...
          '$Ang II_{AT1R-bound}$'; '$Ang II_{AT2R-bound}$'; ...
          '$Ang (1-7)$'; '$Ang IV$'; '$R_{aa}$'; '$R_{ea}$'; ...
          '$\Sigma_{myo}$'; '$\Psi_{AT1R-AA}$'; '$\Psi_{AT1R-EA}$'; ...
          '$\Psi_{AT2R-AA}$'; '$\Psi_{AT2R-EA}$';'$\Phi_{sodin}$'};
   ylim_labels = {'l/min','mmHg','l/min','mEq/min'}; 
   linecolor = {'g';'b';'r'};
   markerstyle = {'o';'+';'*';'X';'s';'d';'^';'p'};
   linestyle = {'-','--','--','--',':',':',':','-'};
   ylimits = [0.06,0.16;...
       80,150;...
       0,0.008;...
       0.05,0.14];
   f = gobjects(length(plot_these),1);
   titles = {'Male','Female','M HTN','F HTN'};
   for i = 1:length(plot_these)
        f(i) = figure;
    
        data = zeros(8,2*2,num_days);
        for gg= sex_to_run
            for AA = 1:length(hyp)
                s(AA,gg) = subplot(1,4,2*(AA-1)+gg);
                disp([AA,gg,2*(AA-1)+gg])
                for mm = treatments_to_run
                    x1=X{gg,mm,AA};  
                    t = T{gg,mm,AA};
                    if t(end) > tchange
                        for dd = 1:(num_days)
                            mask1 = (t>tchange+(dd-1)*1440);
                            mask2 = (t < tchange + dd*1440);
                            mask = mask1&mask2;
                            day = t(mask);
                            dt = day(2:end) - day(1:end-1);
                            x_day = x1(mask,plot_these(i));
                            ave = sum(dt .* x_day(2:end))/(day(end) - day(1));
                            data(mm,2*(AA-1) + gg,dd) = ave;
                        end
                    else
                        data(mm,2*(AA-1) + gg,:) = zeros(1,num_days);
                    end
                    plot(1:num_days,reshape(data(mm,2*(AA-1)+gg,:),[1,num_days]),...
                         sprintf('%s%s',linestyle{mm},markerstyle{mm}),...
                        'Color',linecolor{mod(mm,3)+1},...
                        'LineWidth',2,...
                        'MarkerSize',10,...
                         'MarkerEdgeColor','black',...
                         'MarkerFaceColor',linecolor{mod(mm,3)+1})
                    hold on
                end
                xlim([0.7,num_days+0.3])
                xticks(1:num_days)                       
                ylim(ylimits(i,:))
                title(titles{2*(AA-1)+gg})%names{plot_these(i)}, 'Interpreter','latex')%, 'FontSize',24)
                set(gca,'YGrid','on','Fontsize',18)
                if i==3
                 hline(0.0003,'k','AKI')
                    %hline(0.001)
                end
                    xlabel('Day')
                    xticklabels({'1','2','3'})
                if (2*(AA-1)+gg) ==1
                     ylabel(ylim_labels{i})
                else 
                    yticklabels([])
                end
            end
        end
        legend({'Control','ACEI','Furosemide','NSAID',...
            'A+F','F+N','A+N','A+F+N'})
   end
end
































