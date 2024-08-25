function dy=fun_newton_rand(x,y,w,s,a)
Veff=a*(cos(x)+(s*x-a*sin(x))*sin(x))/(1+(s*x-a*sin(x))^2)+w^2;
dy=[ y(2)
    -(2.0*(s*x-a*sin(x))*(s-a*cos(x))/(1+(s*x-a*sin(x))^2)*y(2)+Veff*y(1))];