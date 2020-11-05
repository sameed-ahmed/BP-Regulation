% This script creates bootstrap replicates of the dataset.
% The dataset is the time course response of MAP to Ang II infusion. 
% It is from Sullivan 2010 Hypertension.

% Input:  none
% Output: saves dataset of 1000 bootstrap replicates.

function create_data_bs_rep

close all
clear

% Species index
spe_ind = 2;

species = {'human', 'rat'   };
sex     = {'male' , 'female'};

for sex_ind = 1:2 % sex

% Load data.
% Data is from Sullivan 2010 Hypertension paper, Figure 1a.
load_data_name = sprintf('%s_AngII_data.xlsx', sex{sex_ind});
if     strcmp(sex{sex_ind}, 'male')
    cell_range = 'B2:P7';
elseif strcmp(sex{sex_ind}, 'female')
    cell_range = 'B2:P8';
end
% data = [individual rat, time MAP]
data = readmatrix(load_data_name,'Range',cell_range);
% data = [11 12; 21 22; 31 32; 41 42; 51 52];
data_points = size(data,1);
time_points = size(data,2);

% Number of bootstrap samples.
num_sample = 1000;

% Create bootstrap replicates of sample sets of rats.
% tic
bs_sample = zeros(data_points, time_points, num_sample);
for j = 1:num_sample
    bs_sample(:,:,j) = datasample(data,data_points);
end
% bs_time = toc

% Compute average over all rats.
% Make all data values as change in MAP from baseline.
% tic
AngII_data_rep = zeros(num_sample,time_points);
for j = 1:num_sample
    AngII_data_rep(j,:) = mean(bs_sample(:,:,j));
    AngII_data_rep(j,:) = AngII_data_rep(j,:) - AngII_data_rep(j,1);
end
% mean_time = toc

% mean(data)
% mean(time_course_data)
% min(time_course_data)
% max(time_course_data)

% Save data.
save_data_name = sprintf('%s_%s_AngII_data_bs_rep.mat', ...
                         species{spe_ind},sex{sex_ind});
save_data_name = strcat('Data/', save_data_name);
save(save_data_name, 'AngII_data_rep', 'num_sample')

end % sex

end % create_data_bs_rep













