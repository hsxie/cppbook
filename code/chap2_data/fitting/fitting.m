close all;clear;clc;
load carsmall; x1 = Weight;
x2 = Horsepower; % Contains NaN data
y = MPG; X = [ones(size(x1)) x1 x2 x1.*x2];
b = regress(y,X); % Removes NaN data
figure; set (gcf,'DefaultAxesFontSize',15);
scatter3(x1,x2,y,'filled'); hold on;
x1fit = min(x1):100:max(x1); x2fit = min(x2):10:max(x2);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.*X2FIT;
mesh(X1FIT,X2FIT,YFIT); view(50,10); 
title('linear regress fitting, example');
xlabel('weight(x_1)'); ylabel('horsepower(x_2)'); zlabel('MPG(y)');
text(0,50,45,['y=',num2str(b(1)),'+(',num2str(b(2)),')x_1+(',...
    num2str(b(3)),')x_2+(',num2str(b(4)),')x_1x_2'],'FontSize',15);

