function dx=chaos(t,x)
a=10;
b=28;
c=8/3;
dx=[a*(x(2)-x(1));x(1)*(b-x(3))-x(2);x(1)*x(2)-c*x(3)];