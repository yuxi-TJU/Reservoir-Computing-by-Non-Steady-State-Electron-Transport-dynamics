% 20240105 M10N4
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
ML = 10;
N = 4;
Vmax = 1.5;
Vmin = -0.5;
freq=linspace(log10(0.01),log10(1000),51);
DT=1./(ML*10.^freq);
for avenumber=1:10
for numberfreq=1:length(freq)
dt=DT(numberfreq)
time=linspace(0,dt*1999,2000);
time2=linspace(0,dt*1999,2000*10);
% ----------------------DATASET----------------------
sample = 8;
sample2=sample*10;
step = 2000;
step2=step/sample*sample2;
Data = zeros(1, 2*step);
Data2 = zeros(1, 2*step2);
p1 = sin(pi*2*(0:sample-1)/sample);
p11=sin(pi*2*(0:sample2-1)/sample2);
p2(1:sample/2) = 1;
p2(sample/2+1:sample) = -1;
p22(1:sample2/2) = 1;
p22(sample2/2+1:sample2) = -1;
number=[];
numberx=1;
for i = 1:2*step/sample
    q = unidrnd(2);
    if q == 1
        Data(sample*(i-1)+1:sample*i) = p1;
        Data2(sample2*(i-1)+1:sample2*i) = p11;
        Label(sample*(i-1)+1:sample*i) = 0;
            if i>step/sample
            number(numberx)=i;
            numberx=numberx+1;
            end 
        else
        Data(sample*(i-1)+1:sample*i) = p2;
        Data2(sample2*(i-1)+1:sample2*i) = p22;
        Label(sample*(i-1)+1:sample*i) = 1;
    end
end

Input = Data(1:step);
Input2 = Data2(1:step2);
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

Input = Data(step+1:end);
Input2 = Data2(step2+1:end);
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
states=states*amplifier;
X = [ones(1,step);states];

Out = Wout*X;
Out=Out;
NRMSE = sqrt(mean((Out(10:end)-Target(10:end)).^2)./var(Target(10:end)));




for m=1:length(Mask(1,:))
    if Mask(1,m)>0
        mm=m*4-3;
for i=1:step/sample-1
deltaP(i)=states(mm,8*i+3)-states(mm,8*i-5);
end
[Delta,o]=sort(deltaP);
DletaPave(numberfreq,avenumber)=sum(Delta(end-4:end))/5;

for i=1:length(number)
   thetaP(i)=states(mm,(number(i)-step/sample)*8-5); 
end
ThetaP(numberfreq,avenumber)=std(thetaP);
    break
    end
end
sprintf('%s', ['freq:', num2str(freq(numberfreq)), ', average number:', num2str(avenumber), ', ML:', num2str(ML)])
end
end
for i=1:length(ThetaP(:,1))
ThetaPave(i)=sum(ThetaP(i,:))/avenumber;
DletaPaveave(i)=sum(DletaPave(i,:))/avenumber;
end
 figure(5)
yyaxis right;
plot(freq,ThetaPave/max(ThetaPave),'LineWidth',3.7);
xlabel('log10(frequency/Hz)','fontsize',18)
ylabel('State richness','fontsize',18)
hold on
yyaxis left;
plot(freq,DletaPaveave/max(DletaPaveave),'LineWidth',3.7);
ylabel('Feedback strength','fontsize',18)
set(gcf,'color','white');
set(gca,'FontSize',18);
hold off

%% 原始
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
ML = 10;
N = 4;
Vmax = 1.5;
Vmin = -0.5;
freq=linspace(log10(0.01),log10(1000),51);
DT=1./(ML*10.^freq);

for numberfreq=1:length(freq)
dt=DT(numberfreq);
time=linspace(0,dt*1999,2000);
time2=linspace(0,dt*1999,2000*10);
sample = 8;
sample2=sample*10;
step = 2000;
step2=step/sample*sample2;
Data = zeros(1, 2*step);
Data2 = zeros(1, 2*step2);
p1 = sin(pi*2*(0:sample-1)/sample);
p11=sin(pi*2*(0:sample2-1)/sample2);
p2(1:sample/2) = 1;
p2(sample/2+1:sample) = -1;
p22(1:sample2/2) = 1;
p22(sample2/2+1:sample2) = -1;
number=[];
numberx=1;
for i = 1:2*step/sample
    q = unidrnd(2);
    if q == 1
        Data(sample*(i-1)+1:sample*i) = p1;
        Data2(sample2*(i-1)+1:sample2*i) = p11;
        Label(sample*(i-1)+1:sample*i) = 0;
            if i>step/sample
            number(numberx)=i;
            numberx=numberx+1;
            end 
        else
        Data(sample*(i-1)+1:sample*i) = p2;
        Data2(sample2*(i-1)+1:sample2*i) = p22;
        Label(sample*(i-1)+1:sample*i) = 1;
    end
end

Input = Data(1:step);
Input2 = Data2(1:step2);
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
Input = Data(step+1:end);
Input2 = Data2(step2+1:end);
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
states=states*amplifier;
X = [ones(1,step);states];

Out = Wout*X;
Out=Out;
NRMSE = sqrt(mean((Out(10:end)-Target(10:end)).^2)./var(Target(10:end)));



for n=1:length(Mask(:,1))
for m=1:length(Mask(1,:))
        mm=m*4-4+n;
    if Mask(n,m)>0

for i=1:step/sample-1
deltaP(i)=states(mm,8*i+3)-states(mm,8*i-5);
end
[Delta,o]=sort(deltaP);
DletaPave(numberfreq,m*4-4+n)=sum(Delta(end-4:end))/5;%feedback strength

for i=1:length(number)
   thetaP(i)=states(mm,(number(i)-step/sample)*8-5); 
end
ThetaP(numberfreq,m*4-4+n)=std(thetaP);%state richness
    end
           if Mask(n,m)<0
for i=1:step/sample-1
deltaP(i)=states(mm,8*i+7)-states(mm,8*i-1);
end
[Delta,o]=sort(deltaP);
DletaPave(numberfreq,m*4-4+n)=sum(Delta(end-4:end))/5;%feedback strength

for i=1:length(number)
   thetaP(i)=states(mm,(number(i)-step/sample)*8-1); 
end
ThetaP(numberfreq,m*4-4+n)=std(thetaP);%state richness
            end
end
end
sprintf('%s', ['freq:', num2str(freq(numberfreq))])
end

for i=1:length(ThetaP(:,1))
ThetaPave(i)=sum(ThetaP(i,:))/40;
DletaPaveave(i)=sum(DletaPave(i,:))/40;
end
figure(5)
yyaxis right; % 激活左边的轴
plot(freq,ThetaPave/max(ThetaPave),'LineWidth',3.7);
xlabel('log10(freq/Hz)','fontsize',18)
ylabel('\sigma(Vp1)V','fontsize',18)
hold on
yyaxis left; % 激活左边的轴
plot(freq,DletaPaveave/max(DletaPaveave),'LineWidth',3.7);
ylabel('\delta Vp1bar','fontsize',18)
set(gcf,'color','white');
set(gca,'FontSize',18);
hold off;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%画图
 figure(1)
plot(time2,Input2, 'b', 'linewidth', 1);
hold on;
plot(time,Input, '.r');
axis([0, 0.025*dt*2000, -1.2, 1.2])
ylabel('Input')
set(gca,'FontName', 'Arial', 'FontSize', 20);

 figure(2)%one state
plot(time,states(1,:), 'b', 'linewidth', 1);
axis([0, 0.025*dt*2000, -1, 40])
ylabel('One state')
set(gca,'FontName', 'Arial', 'FontSize', 20);

figure(3)%all state
plot(time,states(1,:),time,states(5,:),time,states(9,:),time,states(13,:),time,states(17,:),time,states(21,:),time,states(25,:),time,states(29,:),time,states(33,:),time,states(37,:));
axis([0, 0.025*dt*2000, -1, 40])
ylabel('All state')
set(gca,'FontName', 'Arial', 'FontSize', 20);

figure(4)
plot(time,Target, 'k', 'linewidth', 2);
hold on;
plot(time,Out, 'r', 'linewidth',1);
axis([0, dt*2000, -0.2, 1.2])
str1 = '\color{black}Target';
str2 = '\color{red}Output';
lg = legend(str1, str2);
set(lg, 'Orientation', 'horizon');
ylabel('Prediction')
xlabel('Time (\tau)')
set(gca,'FontName', 'Arial', 'FontSize', 20);
set(gcf, 'unit', 'normalized', 'position', [0.2, 0.2, 0.6, 0.35]);
%% 所有态
for n=1:length(Mask(:,1))
for m=1:length(Mask(1,:))
        mm=m*4-4+n;
    if Mask(n,m)>0

for i=1:step/sample-1
deltaP(i)=states(mm,8*i+3)-states(mm,8*i-5);
end
[Delta,o]=sort(deltaP);
DletaPave(numberfreq,m*4-4+n)=sum(Delta(end-4:end))/5;%feedback strength

for i=1:length(number)
   thetaP(i)=states(mm,(number(i)-step/sample)*8-5); 
end
ThetaP(numberfreq,m*4-4+n)=std(thetaP);%state richness
    end
           if Mask(n,m)<0
for i=1:step/sample-1
deltaP(i)=states(mm,8*i+5)-states(mm,8*i-3);
end
[Delta,o]=sort(deltaP);
DletaPave(numberfreq,m*4-4+n)=sum(Delta(end-4:end))/5;%feedback strength

for i=1:length(number)
   thetaP(i)=states(mm,(number(i)-step/sample)*8-3); 
end
ThetaP(numberfreq,m*4-4+n)=std(thetaP);%state richness
            end
end
end
for i=1:length(ThetaP(:,1))
ThetaPave(i)=sum(ThetaP(i,:))/40;
DletaPaveave(i)=sum(DletaPave(i,:))/40;
end

 figure(5)
yyaxis right; % 激活左边的轴
plot(freq,ThetaPave,'LineWidth',3.7);
xlabel('log10(freq/Hz)','fontsize',18)
ylabel('\sigma(Vp1)V','fontsize',18)
hold on
yyaxis left; % 激活左边的轴
plot(freq,DletaPaveave,'LineWidth',3.7);
ylabel('\delta Vp1bar','fontsize',18)
set(gcf,'color','white');
set(gca,'FontSize',18);
hold off
%% 两通道电流差很大情况 回滞大
clear;clc;

eg1=0.30;
lambda1=0.7;
gr=10^13;
gl=10^13;
eg2=0.2;
lambda2=0.8;
gl2=10^10;%%%%%%%%%%%%%%%%
gr2=10^10;%%%%%%%%%%%%%%%%%%

% scanr=0.5;
nH=20;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kp=10*nH;
kn =0.05; 

AQ1=0.6;
HAQ1=0.4;
amplifier=10^10;
ML = 5;
N = 8;
Vmax = 1.5;
Vmin = -0.5;
freq=linspace(log10(0.01),log10(1000),51);
DT=1./(ML*10.^freq);
for avenumber=1:10
for numberfreq=1:length(freq)
dt=DT(numberfreq);%programming speed%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time=linspace(0,dt*1999,2000);
time2=linspace(0,dt*1999,2000*10);
% ----------------------DATASET----------------------
sample = 8;
sample2=sample*10;
step = 2000;
step2=step/sample*sample2;
Data = zeros(1, 2*step);
Data2 = zeros(1, 2*step2);
p1 = sin(pi*2*(0:sample-1)/sample);
p11=sin(pi*2*(0:sample2-1)/sample2);
p2(1:sample/2) = 1;
p2(sample/2+1:sample) = -1;
p22(1:sample2/2) = 1;
p22(sample2/2+1:sample2) = -1;
number=[];
numberx=1;
for i = 1:2*step/sample
    q = unidrnd(2);
    if q == 1
        Data(sample*(i-1)+1:sample*i) = p1;
        Data2(sample2*(i-1)+1:sample2*i) = p11;
        Label(sample*(i-1)+1:sample*i) = 0;
            if i>step/sample
            number(numberx)=i;
            numberx=numberx+1;
            end 
        else
        Data(sample*(i-1)+1:sample*i) = p2;
        Data2(sample2*(i-1)+1:sample2*i) = p22;
        Label(sample*(i-1)+1:sample*i) = 1;
    end
end

% ----------------------TRAIN----------------------
Input = Data(1:step);
Input2 = Data2(1:step2);
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

Input = Data(step+1:end);
Input2 = Data2(step2+1:end);
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
states=states*amplifier;
X = [ones(1,step);states];

Out = Wout*X;
Out=Out;
NRMSE = sqrt(mean((Out(10:end)-Target(10:end)).^2)./var(Target(10:end)));




for m=1:length(Mask(1,:))
    if Mask(1,m)>0
        mm=m*8-7;
for i=1:step/sample-1
deltaP(i)=states(mm,8*i+3)-states(mm,8*i-5);
end
[Delta,o]=sort(deltaP);
DletaPave(numberfreq,avenumber)=sum(Delta(end-4:end))/5;%feedback strength

for i=1:length(number)
   thetaP(i)=states(mm,(number(i)-step/sample)*8-5); 
end
ThetaP(numberfreq,avenumber)=std(thetaP);%state richness
    break
    end
end
sprintf('%s', ['freq:', num2str(freq(numberfreq)), ', average number:', num2str(avenumber), ', ML:', num2str(ML)])
end
end
for i=1:length(ThetaP(:,1))
ThetaPave(i)=sum(ThetaP(i,:))/avenumber;
DletaPaveave(i)=sum(DletaPave(i,:))/avenumber;
end
 figure(5)
yyaxis right; % 激活左边的轴
plot(freq,ThetaPave,'LineWidth',3.7);
xlabel('log10(freq/Hz)','fontsize',18)
ylabel('\sigma(Vp1)V','fontsize',18)
hold on
yyaxis left; % 激活左边的轴
plot(freq,DletaPaveave,'LineWidth',3.7);
ylabel('\delta Vp1bar','fontsize',18)
set(gcf,'color','white');
set(gca,'FontSize',18);
hold off