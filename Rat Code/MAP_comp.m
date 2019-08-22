% This script is a post processing script. It loads the data from some
% different drug administration scenarios to compute the change in mean 
% arterial pressure and plot.

function MAP_comp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Begin user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scenarios
% Normal - Normal conditions
% m_RSNA - male RSNA
% m_AT2R - male AT2R
% m_RAS  - male RAS pars
% m_Reab - male fractional sodium and water reabsorption
% ACEi   - Angiotensin convernting enzyme inhibitor
% AngII  - Ang II infusion
scenario = {'Normal', 'ACEi', 'ARB'};
num_scen = length(scenario);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End user input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Data');
addpath(genpath(mypath))

gender   = {'male', 'female'};

% MAP = [scenario, gender];
deltaMAP = zeros(num_scen,2);

for ss = 1:num_scen % scenario
for gg = 1:2 % gender

% Set name for data file to be loaded based upon gender and scenario.    
load_data_name = sprintf('%s_ss_data_scenario_%s.mat', gender{gg},scenario{ss});

% Load and retrieve data.
load(load_data_name, 'SSdata');
if     strcmp(scenario{ss}, 'Normal')
    SSdata_Normal = SSdata;
end
deltaMAP(ss,gg) = SSdata(42) - SSdata_Normal(42);
clear SSdata;

end % gender
end % scenario

% Plot 
% male ACEi; female ACEi; male ARB; female ARB;
f1 = figure('DefaultAxesFontSize',30, 'pos',[100 100 650 450]);
c1 = categorical({'Normal','ACEi','ARB'});
b1 = bar(c1,deltaMAP);
ylabel('Change in MAP (mmHg)')
legend('Male','Female', 'Location','Southeast')
b1(1).FaceColor = [0.203, 0.592, 0.835];
b1(2).FaceColor = [0.835, 0.203, 0.576];

f2 = figure('DefaultAxesFontSize',30, 'pos',[100 100 650 450]);
c2 = categorical({'ACEi','ARB'});
b2 = bar(c2,deltaMAP(2:end,:));
ylabel('Change in MAP (mmHg)')
legend('Male','Female', 'Location','Southeast')
b2(1).FaceColor = [0.203, 0.592, 0.835];
b2(2).FaceColor = [0.835, 0.203, 0.576];

% % Save figure.
% savefig(f, 'Pma_change_ACEi_ARB.fig')

end



































