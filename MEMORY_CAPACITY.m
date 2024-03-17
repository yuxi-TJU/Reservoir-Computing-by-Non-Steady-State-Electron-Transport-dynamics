%% WAVE PREDICT
clear;clc;
load('NARMA10data');
eg1=0.30;
lambda1=0.7;
gr=10^13;
gl=10^13;
eg2=0.20;
lambda2=0.8;
gl2=10^12;
gr2=10^11;

nH=20;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kp=10*nH;
kn =0.05; 

AQ1=0.6;
HAQ1=0.4;
% timestep2=1;

amplifier=10^10;
ML = 5;
N = 8;
Vmax = 1.5;
Vmin = -0.5;
dt=0.5;%programming speed%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
step = 2000;

Data=data(101:end);



tuo = 2;

MC_C = zeros(tuo,1);
for m = 1:tuo


Input=Data(1:(step))';
Target = Data(m+1-1:(step+m-1))';


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
HAQ1=amoh2AQ(:,i);
    sprintf('%s', ['train:', num2str(i), ', Vmax:', num2str(Vmax), ', ML:', num2str(ML)])
end

states = [];
for i = 1:step
    a = Currents(:, ML*(i-1)+1:ML*i);
    states(:, i) = a(:);
end
states=states*amplifier;
X = [ones(1,step); states];

Wout = Target*pinv(X);


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
HAQ1=amoh2AQ(:,i);
    sprintf('%s', ['train:', num2str(i), ', Vmax:', num2str(Vmax), ', ML:', num2str(ML)])
end

states = [];
for i = 1:step
    a = Currents(:, ML*(i-1)+1:ML*i);
    states(:, i) = a(:);
end
states=states*amplifier;
X = [ones(1,step); states];
Out=Wout*X;

COV=cov(Target,Out);
COV1=cov(Target,Target);
MR_C(m)=COV(1,2)*COV(2,1)/COV1(1,1)/COV(2,2);


fprintf('%f\n',MR_C(m))
end


figure(8)
plot(MR_C)
fprintf('%f\n',sum(MR_C(MR_C>0.1)))
ylabel('memory capacity')
xlabel('tuo')


figure(9);
subplot(2, 1, 1);
plot(Input, 'b', 'linewidth', 1);


ylabel('Input')
set(gca,'FontName', 'Arial', 'FontSize', 20);
subplot(2, 1, 2);
plot(Target, 'k', 'linewidth', 2);
hold on;
plot(Out, 'r', 'linewidth',1);

str1 = '\color{black}Target';
str2 = '\color{red}Output';
lg = legend(str1, str2);
set(lg, 'Orientation', 'horizon');
ylabel('Prediction')
xlabel('Time (\tau)')
set(gca,'FontName', 'Arial', 'FontSize', 20);
set(gcf, 'unit', 'normalized', 'position', [0.2, 0.2, 0.6, 0.35]);

%% 可以用的MC
clear;clc;

eg1=0.30;
lambda1=0.7;
gr=10^13;
gl=10^13;
eg2=0.2;
lambda2=0.8;
gl2=10^12;
gr2=10^11;

nH=20;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kp=10*nH;
kn =0.05; 

AQ1=0.6;
HAQ1=0.4;
amplifier=10^10;
ML = 4;
N = 25;
Vmax = 1.5;
Vmin = -0.5;

dt=0.5;%programming speed%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
step = 1000;
dataset = HenonMap(2*step+1);

tuo = 40;

MC_C = zeros(tuo,1);
for m = 1:tuo
Input = dataset(1:step);
Target = dataset(1+m:step+m);
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
HAQ1=amoh2AQ(:,i);
end

states = [];
for i = 1:step
    a = Currents(:, ML*(i-1)+1:ML*i);
    states(:, i) = a(:);
end
states=states*amplifier;
X = [ones(1,step); states];

Wout = Target*pinv(X);
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
states=states*amplifier;
X = [ones(1,step);states];

Out = Wout*X;
Out=Out;
NRMSE = sqrt(mean((Out(10:end)-Target(10:end)).^2)./var(Target(10:end)));
sprintf('%s',['NRMSE:',num2str(NRMSE)])
COV=cov(Target,Out);
COV1=cov(Target,Target);
MR_C(m)=COV(1,2)*COV(2,1)/COV1(1,1)/COV(2,2);


fprintf('%f\n',MR_C(m))
end
% ----------------------PLOT---------------------
figure(8)
plot(MR_C)
fprintf('%f\n',sum(MR_C(MR_C>0.05)))
ylabel('memory capacity')
xlabel('tuo')


figure(9);
plot(Target, 'k', 'linewidth', 2);
hold on;
plot(Out, 'r', 'linewidth',1);

str1 = '\color{black}Target';
str2 = '\color{red}Output';
lg = legend(str1, str2);
set(lg, 'Orientation', 'horizon');
ylabel('Prediction')
xlabel('Time (\tau)')
set(gca,'FontName', 'Arial', 'FontSize', 20);
set(gcf, 'unit', 'normalized', 'position', [0.2, 0.2, 0.6, 0.35]);

%% 可以用的MC narma10来测试 循环一直到τ=50
clear;clc;%计算：目标函数为u(t-k),输入为u(t)

eg1=0.30;

lambda1=0.7;
gr=10^13;
gl=10^13;
eg2=0.2;
lambda2=0.8;
gl2=10^12;
gr2=10^11;

nH=1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kp=10*nH;
kn =0.05; 
AQ1=0.6;
HAQ1=0.4;
amplifier=10^10;
ML = 10;
N = 4;
Vmax = 1.5;
Vmin = -0.5;
step = 2000;
load('NARMA10data');
dataset=data(50:end)';

tuo = 30;
xaix=linspace(0,tuo,tuo+1);

number1=2;
freq=linspace(log10(0.01),log10(1000),51);
DT=1./(ML*10.^freq);
for number=1:length(freq)%行是频率
    dt=DT(number);
for m = 0:tuo
Target = dataset(1:step);
Input = dataset(1+m:step+m);
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
HAQ1=amoh2AQ(:,i);
end

states = [];
for i = 1:step
    a = Currents(:, ML*(i-1)+1:ML*i);
    states(:, i) = a(:);
end
states=states*amplifier;
X = [ones(1,step); states];

Wout = Target*pinv(X);

% ----------------------TEST----------------------
Input=[];
Target=[];
Target = dataset(1+step+m+5000:step+step+m+5000);

Input = dataset(1+m+step+m+5000:step+m+step+m+5000);

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
states=states*amplifier;
X = [ones(1,step);states];

Out = Wout*X;
Out=Out;
NRMSE(m+1,number) = sqrt(mean((Out(10:end)-Target(10:end)).^2)./var(Target(10:end)));
sprintf('%s',['NRMSE:',num2str(NRMSE(m+1,number))])
COV=cov(Target,Out);
COV1=cov(Input,Target);
MR_C(m+1,number)=COV(1,2)*COV(2,1)/COV1(1,1)/COV(2,2);
MR_C1(m+1,number)=COV(1,2)*COV(2,1)/COV(1,1)/COV(2,2);


sprintf('%s',['number:',num2str(number),'  tuo:',num2str(m),'  MCi:',num2str(MR_C(m+1,number))])
end


end
% 
for i=1:number
MC(i)=sum(MR_C(:,i));
MC1(i)=sum(MR_C1(:,i));
end
figure(8)
plot(xaix,MR_C)
ylabel('memory capacity')
xlabel('tuo')


figure(9)
plot(freq,MC)
ylabel('total memory capacity')
xlabel('frequency')

y=ones(length(xaix));
figure(10)
plot3(xaix,y*freq(1),MR_C(:,1),xaix,y*freq(11),MR_C(:,11),xaix,y*freq(21),MR_C(:,21),xaix,y*freq(31),MR_C(:,31),xaix,y*freq(41),MR_C(:,41),'Color','r','LineWidth',2);
xlabel('x');
ylabel('y');
zlabel('z');

