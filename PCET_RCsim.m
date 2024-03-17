%% Waveform classification
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
AQ1=1;
HAQ1=0;

amplifier=10^10;
ML = 10;
N = 4;
Vmax = 1.5;
Vmin = -0.5;
freq=linspace(log10(0.01),log10(1000),51);
DT=1./(ML*10.^freq);

dt=DT(31);%programming speed%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
for i = 1:2*step/sample
    q = unidrnd(2);
    if q == 1
        Data(sample*(i-1)+1:sample*i) = p1;
        Data2(sample2*(i-1)+1:sample2*i) = p11;
        Label(sample*(i-1)+1:sample*i) = 0;
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
    sprintf('%s', ['test:', num2str(i), ', Vmax:', num2str(Vmax), ', ML:', num2str(ML)])
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

%% Spoken-digit recognition
clc;
clear;
addpath('E:\研究生二年级工作安排\动态器件\Dynamic-memristor-based-reservoir-computing-v1.0.0\Stuka255-Dynamic-memristor-based-reservoir-computing-9193344\Auditory Toolbox\Auditory Toolbox\');


eg1=0.30;
lambda1=0.7;
gr=10^13;
gl=10^13;
eg2=0.20;
lambda2=0.8;
gl2=10^12;
gr2=10^11;

% scanr=0.5;
nH=20;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kp=10*nH;
kn =0.05; 

AQ1=0.6;
HAQ1=0.4;
amplifier=10^10;
Vmax = 1.5;
Vmin = -0.5;
dt=0.5;%programming speed%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ML = 10;
N = 40;
Mask = 2*randi([0,1],64,ML,N)-1;

for i = 1:5
    for j = 1:10
        for k = 1:10
            filename(k+(i-1)*10,j) = {['1',num2str(i),'\','0',num2str(j-1),'f',num2str(i),'set',num2str(k-1),'.wav']};
        end
    end
end

WRR = 0;
TF=zeros(10,10);
for u=1:10
S = [];
for i = 1:10
    r = randperm(size(filename,1));
    res = filename(:,i);
    res = res(r,:);
    S = [S,res];
end

words = 450;
VL = zeros(words);
Target = [];X = [];
q = 0;p = 1;
for j = 1:words   
    q = q+1;
    if q > 10
        q = 1;
        p = p+1;
    end
    
    a = audioread(S{p,q});
    a = resample(a,8000,12500);
    f = LyonPassiveEar(a,8000,250);
    L = zeros(10,length(f(1,:)));
    L(q,:) = ones(1,length(f(1,:)));
    VL(j) = length(f(1,:));
    Target(:,sum(VL(1:j))-VL(j)+1:sum(VL(1:j))) = L;
    
    Input = [];
    for k = 1:N
        for i = 1:VL(j)
            Input(k, ML*(i-1)+1:ML*i) = abs(f(:,i))'*Mask(:,:,k);            
        end
    end
    UL = max(max(Input));
    DL = min(min(Input));
    Input = (Input-DL)/(UL-DL)*(Vmax - Vmin)+Vmin;
    
    memout = [];
    for i = 1:length(Input(1, :))
[amoAQ(:,i),amoh2AQ(:,i),CurrentsAQ(:,i),Currentsh2AQ(:,i),memout(:,i)]=It3pulse(eg1,lambda1,gr,gl,eg2,lambda2,gr2,gl2,dt,Input(:,i),kn,kp,AQ1,HAQ1);

AQ1=amoAQ(:,i);
HAQ1=amoh2AQ(:,i);

    end
    
    for i = 1:VL(j)
        a = memout(:, ML*(i-1)+1:ML*i);
        X(:,sum(VL(1:j))-VL(j)+i) = a(:);
    end
    
    sprintf('%s',['loop:',num2str(u),',train:',num2str(j),',',num2str(u-1),'acc:',num2str(WRR)])
end

Wout = Target*X'*pinv(X*X');

clc;
VL = zeros(words);
Target = [];X=[];
words = 50;q = 0;p = 46;
for j=1:words
    q = q+1;
    if q > 10
        q = 1;
        p = p+1;
    end

    % data preprocess
    a = audioread(S{p,q});
    a = resample(a,8000,12500);
    f = LyonPassiveEar(a,8000,250);
    L = zeros(10,length(f(1,:)));
    L(q,:) = ones(1,length(f(1,:)));
    VL(j) = length(f(1,:));
    Target(:,sum(VL(1:j))-VL(j)+1:sum(VL(1:j))) = L;
    
    % mask process
    Input = [];
    for k=1:N
        for i=1:VL(j)
            Input(k, ML*(i-1)+1:ML*i) = abs(f(:,i))'*Mask(:,:,k);            
        end
    end
    UL = max(max(Input));
    DL = min(min(Input));
    Input = (Input-DL)/(UL-DL)*(Vmax - Vmin)+Vmin;
    
    memout = [];
    for i = 1:length(Input(1, :))
[amoAQ(:,i),amoh2AQ(:,i),CurrentsAQ(:,i),Currentsh2AQ(:,i),memout(:,i)]=It3pulse(eg1,lambda1,gr,gl,eg2,lambda2,gr2,gl2,dt,Input(:,i),kn,kp,AQ1,HAQ1);

AQ1=amoAQ(:,i);
HAQ1=amoh2AQ(:,i);

    end
    
    % states collection
    for i = 1:VL(j)
        a = memout(:, ML*(i-1)+1:ML*i);
        X(:,sum(VL(1:j))-VL(j)+i) = a(:);
    end

    sprintf('%s',['loop:',num2str(u),',test:',num2str(j)])
end

% system output
Y = Wout*X;

% accuracy calculation
Mout = [];
rl = zeros(10,10);
real = zeros(10,words);
for i=1:words
    Mout(:,i) = mean(Y(:,sum(VL(1:i))-VL(i)+1:sum(VL(1:i))),2);
    [~,id] = max(Mout(:,i));
    real(id,i) = 1;
    if mod(i,10) == 0
        rl = rl+real(:,(i/10-1)*10+1:i);
    end
end
TF= TF+rl;
end
WRR = 100*sum(sum(TF.*eye(10,10)))/(u*words);

% ----------------------PLOT----------------------
figure(1);
x = [0 9];y = [0 9];
imagesc(x, y, TF);
ylabel('Predicted output digit')
xlabel('Correct output digit')
title(['Acc: ',num2str(WRR),'%'])
colorbar;
colormap(flipud(hot)); 
set(gca,'FontName', 'Arial', 'FontSize', 15);
FF=linspace(1,length(f(1,:)),length(f(1,:)));
figure(3);
plot(FF,f(1,:),FF,f(5,:),FF,f(15,:),FF,f(25,:),FF,f(35,:))

figure(2);
subplot(2,1,1)
plot(Input(1, :));
ylabel('Input (V)')
axis([0,inf,-inf,inf]);
set(gca,'FontName', 'Arial', 'FontSize', 15);
subplot(2,1,2)
plot(memout(1, :), 'r');
xlabel('Time step')
ylabel('Output (uA)')
axis([0,inf,-inf,inf]);
set(gca,'FontName', 'Arial', 'FontSize', 15);

