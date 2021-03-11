function man_figs_graph_abs

close all

%% Experiment figure

% Number of points in plot
N = 100;
% Number of individuals in population
M = 1000;

t = linspace(0,10,N)';

am = 7.43;
bm = 4.64;
% ---
af = 3.8;
bf = 3;

noisem = randn(2,M)/10;
apopm = am + am .* noisem(1,:);
bpopm = am + am .* noisem(2,:);
% ---
noisef = randn(2,M)/10;
apopf = af + af .* noisef(1,:);
bpopf = af + af .* noisef(2,:);

ypopm = zeros(N,M);
ypopf = zeros(N,M);
for i = 1:M
    ypopm(:,i) = apopm(i)*t ./ (bpopm(i) + t);
    ypopf(:,i) = apopf(i)*t ./ (bpopf(i) + t);
end
ypopm_mean  = mean(ypopm,  2);
ypopm_std   = std (ypopm,0,2);
ypopm_lower = ypopm_mean - 2*ypopm_std;
ypopm_upper = ypopm_mean + 2*ypopm_std;
% ---
ypopf_mean  = mean(ypopf,  2);
ypopf_std   = std (ypopf,0,2);
ypopf_lower = ypopf_mean - 2*ypopf_std;
ypopf_upper = ypopf_mean + 2*ypopf_std;

exp_fig = figure;
plot(t,ypopm_mean ,'-' , 'Color',[0.203, 0.592, 0.835], 'LineWidth',3)
hold on
plot(t,ypopf_mean ,'-' , 'Color',[0.835, 0.203, 0.576], 'LineWidth',3)
plot(t,ypopm_lower,'--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',1)
plot(t,ypopm_upper,'--', 'Color',[0.203, 0.592, 0.835], 'LineWidth',1)
plot(t,ypopf_lower,'--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',1)
plot(t,ypopf_upper,'--', 'Color',[0.835, 0.203, 0.576], 'LineWidth',1)
hold off
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
grid on

%% Distribution figure

distm = random('beta',2,5,[M,1]);
distf = random('beta',2,5,[M,1])+0.25;

dist_fig = figure;
hm = histogram(distm);
hm.FaceColor = [0.203, 0.592, 0.835];
hold on
hf = histogram(distf);
hf.FaceColor = [0.835, 0.203, 0.576];
hold off
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])

%% Save figures

save_data_name = sprintf('exp.png');
exportgraphics(exp_fig, save_data_name)
% 
save_data_name = sprintf('dist.png');
exportgraphics(dist_fig, save_data_name)

end


























