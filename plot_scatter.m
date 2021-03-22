
%% read in data
human = 1;% 1 for human, 0 for rat
species = {'rat','human'};
gender   = {'male', 'female'};
num_days = 2;
tchange = 100;

models = [0,0,0;...
    1,0,0;...
    0,1,0;...
    0,0,1;...
    1,1,0;...
    0,1,1;...
    1,0,1;...
    1,1,1];
%models = 2*models; %for high dose simulations
treatments_to_run = 1:8;
sex_to_run = 1:2;
hyp = [1,2.5];
T = {};
X = {};
for gg = sex_to_run %sex
    for mm= treatments_to_run 
        for AA = 1:length(hyp)
            filename = '%s%s_%s_%s_%s_%s_%s.mat';
            filename = sprintf( filename,'Human_Data/',species{human+1},gender{gg},num2str(models(mm,1)),num2str(models(mm,2)),num2str(models(mm,3)),num2str(hyp(AA)));
            load(filename, 't','x');
            T{gg,mm,AA} = t;
            X{gg,mm,AA} = x;
            tf = t(end);
        end
    end
end




names  = {'$rsna$'; '$\alpha_{map}$'; '$\alpha_{rap}$'; '$R_{r}$'; ...
      '$\beta_{rsna}$'; '$\Phi_{rb}$'; 'GFR'; '$P_{f}$'; ...
      '$P_{gh}$'; '$\Sigma_{tgf}$'; '$\Phi_{filsod}$'; ...
      '$\Phi_{pt-sodreab}$'; '$\eta_{pt-sodreab}$'; ...
      '$\gamma_{filsod}$'; '$\gamma_{at}$'; '$\gamma_{rsna}$'; ...
      '$\Phi_{md-sod}$'; '$\Phi_{dt-sodreab}$'; ...
      '$\eta_{dt-sodreab}$'; '$\psi_{al}$'; '$\Phi_{dt-sod}$'; ...
      '$\Phi_{cd-sodreab}$'; '$\eta_{cd-sodreab}$'; ...
      '$\lambda_{al}$';'$\lambda_{dt}$'; '$\lambda_{anp}$'; '$\Phi_{u-sod}$'; ...
      '$\Phi_{win}$'; '$V_{ecf}$'; '$V_{b}$'; '$P_{mf}$'; ...
      '$\Phi_{vr}$'; '$\Phi_{co}$'; '$P_{ra}$'; '$vas$'; ...
      '$vas_{f}$'; '$vas_{d}$'; '$R_{a}$'; '$R_{ba}$'; '$R_{vr}$'; ...
      '$R_{tp}$'; '$P_{ma}$'; '$\epsilon_{aum}$'; '$a_{auto}$'; ...
      '$a_{chemo}$'; '$a_{baro}$'; '$C_{adh}$'; '$N_{adh}$'; ...
      '$N_{adhs}$'; '$\delta_{ra}$'; '$\Phi_{u}$'; '$M_{sod}$'; ...
      '$C_{sod}$'; '$\nu_{md-sod}$'; '$\nu_{rsna}$'; '$C_{al}$'; ...
      '$N_{al}$'; '$N_{als}$'; '$\xi_{k/sod}$'; '$\xi_{map}$'; ...
      '$\xi_{at}$'; '$\hat{C}_{anp}$'; '$AGT$'; '$\nu_{AT1}$'; ...
      '$R_{sec}$'; '$PRC$'; '$PRA$'; '$Ang I$'; '$Ang II$'; ...
      '$Ang II_{AT1R-bound}$'; '$Ang II_{AT2R-bound}$'; ...
      '$Ang (1-7)$'; '$Ang IV$'; '$R_{aa}$'; '$R_{ea}$'; ...
      '$\Sigma_{myo}$'; '$\Psi_{AT1R-AA}$'; '$\Psi_{AT1R-EA}$'; ...
      '$\Psi_{AT2R-AA}$'; '$\Psi_{AT2R-EA}$';'$\Phi_{sodin}$';...
      '$\Phi_{t-wreab}$'; ...
      '$\mu_{adh}$'; '$\mu_{Na}$';};
plot_these = [7,42,83,66,80,81,82,52,68,55,73];%,52,30,11,17,21,27,68,55];
ylim_labels = {'l/min','mmHg','l','fmol/ml/min','l/min','','','mEq/l','fmol/ml','mEq/L','l','mEq/min','mEq/min','mEq/min','mEq/min'}; 
label_names = {'GFR','MAP','Urine Volume','PRA','Tubular water reabsorption','Mu_{adh}','Mu_{na}','Plasma Na concentration','Ang II','ALD','Blood Volume','filtered Na','MD Na','DT Na','Urine Na Flow'};
linecolor = {'g';'b';'r'};
markerstyle = {'o';'+';'*';'X';'s';'d';'^';'p'};
linestyle = {'-','--','--','--',':',':',':','-'};
ylimits = [0.06,0.15;...
   80,155;...
   0,12;...
   15,55;140,160;0,20;12,22;2.5,7;1.3,2.6;0,0.6];
f = gobjects(length(plot_these),1);
titles = {'Male','Female','M HTN','F HTN'};
data = zeros(length(plot_these),8,2*2,num_days);
for i = 1:length(plot_these)
    if i <= 4
        f= figure;
    end
    for gg= sex_to_run
        for AA = 1:length(hyp)
            if i<=4
                s(AA,gg) = subplot(1,4,2*(AA-1)+gg);
            end
            %disp([AA,gg,2*(AA-1)+gg])
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
                        
                        if plot_these(i) == 83
                            %plot urine flow volume for each day
                            data(i,mm,2*(AA-1) + gg,dd) = sum(dt .* x_day(2:end));%ave;
                        else
                            %plot value daily average
                            data(i,mm,2*(AA-1) + gg,dd) = ave;
                        end
                        
                    end
                else
                    data(i,mm,2*(AA-1) + gg,:) = zeros(1,num_days);
                end
                if i <=4 %only plot first 4 variables
                    plot(1:num_days,reshape(data(i,mm,2*(AA-1)+gg,:),[1,num_days]),...
                         sprintf('%s%s',linestyle{mm},markerstyle{mm}),...
                        'Color',linecolor{mod(mm,3)+1},...
                        'LineWidth',2,...
                        'MarkerSize',10,...
                         'MarkerEdgeColor','black',...
                         'MarkerFaceColor',linecolor{mod(mm,3)+1})
                    hold on
                end
            end
              if i <=4 
                xlim([0.7,num_days+0.3])
                xticks(1:num_days)                       
                ylim(ylimits(i,:))
                title(titles{2*(AA-1)+gg})%names{plot_these(i)}, 'Interpreter','latex')%, 'FontSize',24)
                if (2*(AA-1)+gg) ==1
                     ylabel(sprintf('%s (%s)',label_names{i},ylim_labels{i}))
                else 
                    yticklabels([])
                end
                set(gca,'YGrid','on','Fontsize',16)

                if i==3
                 hline(0.0003*1440,'k','AKI')
                end
                    xlabel('Day')
                    %xticklabels(num2str(1:numdays))
              end

        end
    end
 if i <= 4
    %legend({'Control','ACEI','Furosemide','NSAID',...
    %    'A+F','F+N','A+N','A+F+N'})
     saveas(f, sprintf('Human_Figures/%s_scatter.png',label_names{i}))
     savefig(f, sprintf('Human_Figures/%s_scatter.fig',label_names{i}))
 end
end

% Data values used in the paper

% rise in MAP from hypertension (compare steady state controls)
disp("Male rise in MAP from hypertension")
data(2,1,3,1) - data(2,1,1,1)
disp("Fenale rise in MAP from hypertension")
data(2,1,4,1) - data(2,1,2,1)

disp("Male percent change in GFR from hypertension")
(data(1,1,3,1) - data(1,1,1,1))/data(1,1,1,1)*100
disp("Female percent change in GFR from hypertension")
(data(1,1,4,1) - data(1,1,1,1))/data(1,1,2,1)*100

disp("Male normotensive percent change from ACEI")
(data(:,2,1,2) - data(:,1,1,2))
(data(:,2,1,2) - data(:,1,1,2))./data(:,1,1,2)*100
disp("Male hypertensive percent change from ACEI")
(data(:,2,3,2) - data(:,1,3,2))
(data(:,2,3,2) - data(:,1,3,2))./data(:,1,3,2)*100
disp("Female hypertensive percent change from ACEI")
(data(:,2,4,1) - data(:,1,4,1))./data(:,1,4,1)*100


disp("Male normotensive values during furosemide treatment")
(data(:,3,1,1))
disp("Male normotensive change from furosemide")
(data(:,3,1,1) - data(:,1,1,1))
disp("Male normotensive change from furosemide")
(data(:,3,1,1) - data(:,1,1,1))./data(:,1,1,1)*100
disp("Female normotensive values during furosemide treatment")
(data(:,3,2,1) )
disp("Female normotensive change from furosemide")
(data(:,3,2,1) - data(:,1,2,1))
disp("Female normotensive perecnt change from furosemide")
(data(:,3,2,1) - data(:,1,2,1))./data(:,1,2,1)*100

disp("Male hypertensive percent change from furosemide")
(data(:,3,3,1) - data(:,1,3,1))./data(:,1,3,1)*100
disp("Female hypertensive percent change from furosemide")
(data(:,3,4,1) - data(:,1,4,1))./data(:,1,4,1)*100

disp("Male normotensive percent change from NSAID")
(data(:,4,1,1) - data(:,1,1,1))./data(:,1,1,1)*100
disp("Female normotensive percent change from NSAID")
(data(:,4,2,1) - data(:,1,2,1))./data(:,1,2,1)*100

disp("Male hypertensive percent change from triple treatment")
(data(:,8,3,1) - data(:,1,3,1))./data(:,1,3,1)*100
disp("Female hypertensive percent change from triple treatment")
(data(:,8,4,1) - data(:,1,4,1))./data(:,1,4,1)*100

disp("Male normotensive percent change from triple treatment")
(data(:,8,1,1) - data(:,1,1,1))./data(:,1,1,1)*100
disp("Female normotensive percent change from triple treatment")
(data(:,8,2,1) - data(:,1,2,1))./data(:,1,2,1)*100




