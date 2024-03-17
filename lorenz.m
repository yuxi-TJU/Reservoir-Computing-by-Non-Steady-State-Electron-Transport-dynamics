clc;
clear;
time =linspace(0,150,50001)
ini=[0.01,0.01,0.01];
% [t x]=ode45('chaos',[0,150],ini);
[t x]=ode45('chaos',time,ini);
% plot(t,x(:,1));
% grid on
plot3(x(:,1),x(:,2),x(:,3));