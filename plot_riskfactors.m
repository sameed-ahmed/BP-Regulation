
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
    1,1,1;...
    2,2,0;...
    0,2,2;...
    2,0,2;...
    2,2,2];
%models = 2*models;
treatments_to_run = 5:12;
sex_to_run = 1:2;
hyp = [1,2.5];


AA = 2;
endings = {'';'_impairedmyo';'_lowwaterintake';'_lowwaterintake_impairedmyo';};
filename = '%s%s_%s_%s_%s_%s_%s%s.mat';

T = zeros(8,4);
X = zeros(8,4,3);
plot_these = [7,42,83];
ylim_labels = {'GFR (l/min)','MAP (mmHg)','Urine Flow (l/min)'};
savenames = {'GFR','MAP','urineflow'};
y_limits = [0,0.15;0,160;0,0.01];
for gg = sex_to_run %sex
    for mm= treatments_to_run 
        for name_end = 1:length(endings)
            %disp({'Data',mm,name_end})
            filename2 = sprintf( filename,'Human_Data/',species{human+1},gender{gg},num2str(models(mm,1)),num2str(models(mm,2)),num2str(models(mm,3)),num2str(hyp(AA)),endings{name_end});
            %disp({mm,filename2})
            load(filename2, 't','x');
            
            mask = (t < tchange + 2*1440);
            day = t(mask);
            x_day = x(mask,:);
            T(mm-4,name_end) = day(end);
            X(mm-4,name_end,:) = x_day(end,plot_these);
        end
    end
    
    for i=1:length(plot_these)
        f=figure('Position', [100 100 650 500]);
        data = zeros(4,8);
        data(:,1:4) = X(1:4,:,i);
        data(:,5:8) = X(5:8,:,i);

        bar(data)
        hold on
        xticklabels({'ACEI\newline+ Furosemide';'Furosemide\newline+ NSAID';'ACEI\newline+ NSAID';'Triple'})
        %legend({'Normal Conditions';'Impaired Myo';'Low Water';'Impaired Myo + Low Water'; 'High Dose';'High Dose + Impaired Myo';'High Dose + Low Water';'High Dose + Impaired Myo + Low Water'})
        ylabel(ylim_labels{i})

        if i==1
            if gg==1
                hline(0.1471,'k','Control')
            else
                hline(0.1421,'k','Control')
            end
        elseif i==2
            if gg==1
                hline(155.3558,'k','Control')
            else
                hline(141.4468,'k','Control')
            end
        elseif i==3
            hline(0.001,'k','Control')
            hline(0.0003,'k','AKI')
        end
        ylim(y_limits(i,:))
        set(gca,'YGrid','on','Fontsize',18)
        savefig(f, sprintf('Human_Figures/%s_%s_riskfactors.fig',gender{gg},savenames{i}))
        saveas(f,sprintf('Human_Figures/%s_%s_riskfactors.png',gender{gg},savenames{i}))
    end
end
