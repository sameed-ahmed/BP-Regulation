x = [2,4];
y = [1,2,3,4;...
     5,6,7,8];

figure;
barh(x,y,'grouped'); 
title('Grouped Style')

figure;
bar(x,y,'grouped'); 
title('Grouped Style')

% close all
% u = [1 2 3]; v = [1 2; 
%                   1 2; 
%                   1 2];
% bar(u,v(:,1),v(:,2))









