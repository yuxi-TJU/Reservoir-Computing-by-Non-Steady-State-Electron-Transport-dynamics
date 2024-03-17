clear;clc;

eg1=0.30;
lambda1=0.7;
gr=10^13;
gl=10^13;
eg2=0.20;
lambda2=0.8;
gl2=10^12;
gr2=10^11;
NH=[0.1,1,5,10,20,30,40,50,100];
DT=[0.0001,0.001,0.005,0.01,0.05,0.1,0.2,0.5,1,10];
% scanr=0.5;
for m=1:length(NH)
    for n=1:length(DT)
clear AQ1 HAQ1 kp kn Input_ex
        
nH=NH(m);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kp=10*nH;
kn =0.05; 

AQ1=0.6;
HAQ1=0.4;

ML = 4;
N = 10;
Vmax = 1.5;
Vmin = -0.5;
dt=DT(n);%programming speed%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sample = 8;
step = 2000;
Data = zeros(1, 2*step);
p1 = sin(pi*2*(0:sample-1)/sample);
p2(1:sample/2) = 1;
p2(sample/2+1:sample) = -1;
for i = 1:2*step/sample
    q = unidrnd(2);
    if q == 1
        Data(sample*(i-1)+1:sample*i) = p1;
        Label(sample*(i-1)+1:sample*i) = 0;
    else
        Data(sample*(i-1)+1:sample*i) = p2;
        Label(sample*(i-1)+1:sample*i) = 1;
    end
end

Input = Data(1:step);

Target = Label(1:step);

Mask = 2*unidrnd(2, N, ML)-3;
Input_ex = [];
for j = 1:N
    for i = 1:step
        Input_ex(j, (i-1)*ML+1:ML*i) = Input(i)*Mask(j, :);
    end
end
UL = max(max(Input_ex));
DL = min(min(Input_ex));
Input_ex = (Input_ex-DL)/(UL-DL)*(Vmax - Vmin)+Vmin;

for i = 1:length(Input_ex(1, :))

[amoAQ(:,i),amoh2AQ(:,i),CurrentsAQ(:,i),Currentsh2AQ(:,i),Currents(:,i)]=It3pulse(eg1,lambda1,gr,gl,eg2,lambda2,gr2,gl2,dt,Input_ex(:,i),kn,kp,AQ1,HAQ1);

AQ1=amoAQ(:,i);
HAQ1=amoh2AQ(:,i);end
states = [];
for i = 1:step
    a = Currents(:, ML*(i-1)+1:ML*i);
    states(:, i) = a(:);
end
X = [ones(1,step); states];
Wout = Target*pinv(X);% 
Input = Data(step+1:end);


Target = Label(step+1:end);


Input_ex = [];
for j = 1:N
    for i = 1:step
        Input_ex(j, (i-1)*ML+1:ML*i) = Input(i)*Mask(j, :);
    end
end
UL = max(max(Input_ex));
DL = min(min(Input_ex));
Input_ex = (Input_ex-DL)/(UL-DL)*(Vmax - Vmin)+Vmin;


Currents = [];
states = [];

for i = 1:length(Input_ex(1, :))


[amoAQ(:,i),amoh2AQ(:,i),CurrentsAQ(:,i),Currentsh2AQ(:,i),Currents(:,i)]=It3pulse(eg1,lambda1,gr,gl,eg2,lambda2,gr2,gl2,dt,Input_ex(:,i),kn,kp,AQ1,HAQ1);

AQ1=amoAQ(:,i);
HAQ1=amoh2AQ(:,i);
end
for i = 1:step
    a = Currents(:, ML*(i-1)+1:ML*i);
    states(:,i) = a(:);
end
X = [ones(1,step);states];


Out = Wout*X;
NRMSE = sqrt(mean((Out(10:end)-Target(10:end)).^2)./var(Target(10:end)));
sprintf('%s',['NRMSE:',num2str(NRMSE),'nh=',num2str(nH),'dt=',num2str(dt)])

nrmse(m,n)=NRMSE;

end
end
for i=1:length(NH)
    for j=length(DT)
sprintf('%s',['NRMSE:',num2str(nrmse(i,j))])
    end
end
