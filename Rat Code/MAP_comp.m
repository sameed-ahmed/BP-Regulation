% This script is a post processing script. It loads the data from some
% different scenarios to compute the change in mean arterial pressure.

function MAP_comp

close all

gender   = {'male', 'female'};

% MAP = [drug, gender];
%      Male    Female
% ACEi
% ARB
deltaMAP = zeros(2,2);

for gg = 1:2 % gender

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

if     strcmp(gender{gg}, 'male')
    load(  'male_ss_data_scenario_Normal.mat', 'SSdata');
    SSdata_Normal = SSdata;
    clear SSdata;
elseif strcmp(gender{gg}, 'female')
    load('female_ss_data_scenario_Normal.mat', 'SSdata');
    SSdata_Normal = SSdata;
    clear SSdata;
end

if     strcmp(gender{gg}, 'male')
    load(  'male_ss_data_scenario_ACEi.mat', 'SSdata');
    SSdata_ACEi = SSdata;
    clear SSdata;
elseif strcmp(gender{gg}, 'female')
    load('female_ss_data_scenario_ACEi.mat', 'SSdata');
    SSdata_ACEi = SSdata;
    clear SSdata;
end

if     strcmp(gender{gg}, 'male')
    load(  'male_ss_data_scenario_ARB.mat', 'SSdata');
    SSdata_ARB = SSdata;
    clear SSdata;
elseif strcmp(gender{gg}, 'female')
    load('female_ss_data_scenario_ARB.mat', 'SSdata');
    SSdata_ARB = SSdata;
    clear SSdata;
end

% MAP = [drug, gender];
%      Male    Female
% ACEi
% ARB
deltaMAP(1,gg) = SSdata_ACEi(41) - SSdata_Normal(41);
deltaMAP(2,gg) = SSdata_ARB (41) - SSdata_Normal(41);

end % gender

deltaMAP;

% male ACEi; female ACEi; male ARB; female ARB;
f = figure('DefaultAxesFontSize',30, 'pos',[100 100 650 450]);
c = categorical({'ACEi','ARB'});
b = bar(c,deltaMAP);
ylabel('Change in MAP (mmHg)')
legend('Male','Female', 'Location','Southeast')

b(1).FaceColor = [0.203, 0.592, 0.835];
b(2).FaceColor = [0.835, 0.203, 0.576];


% Save figure.

savefig(f, 'Pma_change_ACEi_ARB.fig')

end



































