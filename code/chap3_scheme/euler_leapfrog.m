close all; clear;clc;
x0=1; v0=0; dt=0.05; nt=1000; omega2=1;
x_euler(1)=x0; % x_euler(n)=x(n)
v_euler(1)=v0; % v_euler(n)=x(n)
x_leapfrog(1)=x0-0.5*dt*v0; % x_leapfrog(n)=x(n-0.5)
x_leapfrog(2)=x0+0.5*dt*v0;
v_leapfrog(1)=v0; % v_leapfrog(n)=v(n)
t(1)=0;
for it=1:nt
    t(it+1)=t(it)+dt;
    v_euler(it+1)=v_euler(it)-omega2*x_euler(it)*dt;
    x_euler(it+1)=x_euler(it)+v_euler(it)*dt;
    v_leapfrog(it+1)=v_leapfrog(it)-omega2*x_leapfrog(it+1)*dt;
    x_leapfrog(it+2)=x_leapfrog(it+1)+v_leapfrog(it+1)*dt;
end
figure; set(gcf,'DefaultAxesFontSize',15);
plot(t,x_euler,'g--',t,x_leapfrog(2:end),'r','LineWidth',2);
xlabel('t');ylabel('x');title('x-t');
legend('Euler','Leapfrog',2); legend('boxoff');