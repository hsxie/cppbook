clear;clc;
t = 0:.05:500;
x0 = sin(t)+1.1.*sin(1.1*t); % signal data

y=x0; x=t;

interpMethod='spline'; %'linear', 'spline', 'cubic'
extrMaxValue = y(find(diff(sign(diff(y)))==-2)+1);
extrMaxIndex = find(diff(sign(diff(y)))==-2)+1;
extrMinValue = y(find(diff(sign(diff(y)))==+2)+1);
extrMinIndex = find(diff(sign(diff(y)))==+2)+1;
up = extrMaxValue;
up_x = x(extrMaxIndex);
down = extrMinValue;
down_x = x(extrMinIndex);
up = interp1(up_x,up,x,interpMethod); 
down = interp1(down_x,down,x,interpMethod);

x1u=up; x1d=down;

plot(t,x0,'g',t,x1u,'r',t,x1d,'r');
title('direct evenlope analysis');
legend('original data','evenlope up','evenlope down');
