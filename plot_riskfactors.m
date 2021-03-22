
%% read in data
human = 1;% 1 for human, 0 for rat
species = {'rat','human'};
gender   = {'male', 'female'};
num_days = 2;
tchange = 100;

models = [%0,0,0;...
%    1,0,0;...
%    0,1,0;...
%    0,0,1;...
    1,1,0;...
    0,1,1;...
    1,0,1;...
    1,1,1;...
    2,2,0;...
    0,2,2;...
    2,0,2;...
    2,2,2];
%models = 2*models;
treatments_to_run = 1:8;%[1,5:12];
sex_to_run = 1:2;
hyp = [1,2.5];


AA = 2;
endings = {'';'_impairedmyo';'_lowwaterintake';'_lowwaterintake_impairedmyo';};
filename = '%s%s_%s_%s_%s_%s_%s%s.mat';
%water84,80,81,82];
%GFR factors
%6,73,74,75,76:79,5,10];
plot_these = [7,42,83];%,80];%,84,81,82];
T = zeros(8,4);
X = zeros(8,4,length(plot_these));
contol = zeros(2,84);
%cumulative_urine = zeros(8,4,3);

ylim_labels = {'GFR (l/min)','MAP (mmHg)','urine volume over two days(l)'};
savenames = {'GFR','MAP','urinevol'};
y_limits = [0,0.15;0,150;0,21];
for gg = sex_to_run %sex
    %read in control
    filename_ss = '%s%s_%s_ss_%s_%s_%s_rsna%s.mat';
    filename_ss2 = sprintf( filename_ss,'Human_Data/',species{human+1},gender{gg},num2str(0),num2str(0),num2str(0),num2str(rsna(AA)));
    %disp({mm,filename2})
    load(filename_ss2, 'SSdata');
    control(gg,:) = SSdata;
    for mm= treatments_to_run 
        for name_end = 1:length(endings)
            %disp({'Data',mm,name_end})
            filename2 = sprintf( filename,'Human_Data/',species{human+1},gender{gg},num2str(models(mm,1)),num2str(models(mm,2)),num2str(models(mm,3)),num2str(hyp(AA)),endings{name_end});
            %disp({mm,filename2})
            load(filename2, 't','x');
            
            mask = (t < tchange + 2*1440);
            day = t(mask);
            x_day = x(mask,:);
            T(mm,name_end) = day(end);
            X(mm,name_end,:) = x_day(end,plot_these);

            mask2 = (t > tchange);
            day2 = t(mask2);
            %disp((day2(end)-100)/1440)
            x_2days = x(mask2,:);
            dt = day2(2:end) - day2(1:end-1);
            %size(x_2days(2:end,83))
            %size(dt)
            X(mm,name_end,3) = sum(x_2days(2:end,83).*dt);
            %X(mm-4,:,4)
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
        %legend({'Normal Conditions';'Impaired Myo';'Low Water';'Impaired Myo + Low Water'; 'High Sensitivity';'High Sensitivity + Impaired Myo';'High Sensitivity + Low Water';'High Sensitivity + Impaired Myo + Low Water'})
        ylabel(ylim_labels{i})
        
        if plot_these(i)==83
            hline(0.0003*1440,'k','AKI')
            hline(control(gg,plot_these(i))*1440,'k','Control')
        else
            hline(control(gg,plot_these(i)),'k','Control')
        end
%         if i==1
%             hline(control(gg,7),'k','Control')
%             if gg==1
%                 hline(0.1471,'k','Control')
%             else
%                 hline(0.1421,'k','Control')
%             end
%         elseif i==2
%             if gg==1
%                 hline(155.3558,'k','Control')
%             else
%                 hline(141.4468,'k','Control')
%             end
%         elseif i==3
%             hline(0.001*1440,'k','Control')
%             hline(0.0003*1440,'k','AKI')
%         end
        if i==3
            hline(0.0003*1440,'k','AKI')
        end
        ylim(y_limits(i,:))
        set(gca,'YGrid','on','Fontsize',18)
        savefig(f, sprintf('Human_Figures/%s_%s_riskfactors.fig',gender{gg},savenames{i}))
        saveas(f,sprintf('Human_Figures/%s_%s_riskfactors.png',gender{gg},savenames{i}))
    end
    disp(' Reduction in GFR from ACEI + Furosemide')
    (X(1,1,1) - control(gg,7))/control(gg,7)*100
    disp(' Reduction in GFR from ACEI + Furosemide with impaired myogenic response')
    (X(1,2,1) - control(gg,7))/control(gg,7)*100
    disp(' Reduction in GFR from ACEI + Furosemide with impaired myogenic response and low water intake')
    (X(1,4,1) - control(gg,7))/control(gg,7)*100
    disp(' Reduction in GFR from triple treatment with high drug sensitivity compared to normal sensitivity')
    X(8,1,1)
    (X(8,1,1) - X(4,1,1))/X(4,1,1)*100
%     X(1,1,3)/X(1,1,1)
%     X(5,1,3)/X(5,1,1)
%     X(10,1,3)/X(10,1,1) 
%     disp([X(1,1,1),X(1,1,3),X(1,1,4),X(1,1,4)/X(1,1,1)])
%     disp([X(5,1,1),X(5,1,3),X(5,1,4),X(5,1,4)/X(5,1,1)])
%     disp([X(10,1,1),X(10,1,3),X(10,1,4),X(10,1,4)/X(10,1,1)]) 
end
