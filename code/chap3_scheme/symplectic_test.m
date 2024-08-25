% Hua-sheng XIE, 2015-04-29 19:14
close all; clear; clc;
nt=1000; dt=0.05; p1(1)=1; q1(1)=0; p2=p1; q2=q1; t=(1:nt)*dt;
figure('DefaultAxesFontSize',15);
for it=1:nt-1
    p1(it+1)=p1(it)-q1(it)*dt; q1(it+1)=q1(it)+p1(it)*dt;
    p2(it+1)=p2(it)-q2(it)*dt; q2(it+1)=q2(it)+p2(it+1)*dt;
end
subplot(1,2,1); plot(t,q1,t,q2,'r--','Linewidth',2);
legend('Euler','Symplectic',2); legend('boxoff');
xlabel('t'); ylabel('q'); box on;
subplot(1,2,2); plot(p1,q1,p2,q2,'r--','Linewidth',2); 
xlabel('p'); ylabel('q'); box on;
