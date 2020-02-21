close all
clear

data_size = 50;
num_sample = 1000;

% data = [1;2;3;4;5;6;7;8;9;10];
data = random('Normal',3,1,data_size,1);
% data = random('Weibull',1,0.5,data_size,1);
x_bar = mean(data);

tic
bs_sample1 = zeros(data_size, num_sample);
for j = 1:num_sample
for i = 1:data_size
    ind = random('Discrete Uniform',data_size);
    bs_sample1(i,j) = data(ind);
end
end
toc

% bs_sample1 = sort(bs_sample1,1);

x_bar_star_data1 = mean(bs_sample1);
x_bar_star1      = mean(x_bar_star_data1);

% figure
% histogram(data)
% hold on
% plot([x_bar,x_bar],[0,data_size/5])
% hold off
% 
% figure
% histogram(x_bar_star_data)
% hold on
% plot([x_bar_star,x_bar_star],[0,num_sample/100])
% % plot([x_bar,x_bar],[0,num_sample/10])
% % legend('', 'bs mean', 'sam mean')
% hold off

tic
[~,bs_sample2] = bootstrp(num_sample,@mean,data);
for i = 1:num_sample
    bs_sample2(:,i) = data(bs_sample2(:,i));
end
toc

x_bar_star_data2 = mean(bs_sample2);
x_bar_star2      = mean(x_bar_star_data2);

% bs_sample2 = sort(bs_sample2,1);

% % % Fastest % % %
tic
bs_sample3 = zeros(data_size, num_sample);
for j = 1:num_sample
    bs_sample3(:,j) = datasample(data,data_size);
end
toc

x_bar_star_data3 = mean(bs_sample3);
x_bar_star3      = mean(x_bar_star_data3);










