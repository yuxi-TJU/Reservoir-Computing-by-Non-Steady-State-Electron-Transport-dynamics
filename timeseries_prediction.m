% henon预测
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
ML = 5;
N = 25;
Vmax = 1.5;
Vmin = -0.5;


dt=0.5;%或者0.01/programming speed%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
step = 500;
dataset = HenonMap(2*step+10);

Input = dataset(1:step+1);

Target = Input(2:end);

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

Input= dataset(step+2:2*step+3);
Target=Input(2:step+1);


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

NRMSE = sqrt(mean((Out(10:end)-Target(10:step)).^2)./var(Target(10:step)));
sprintf('%s',['NRMSE:',num2str(NRMSE)])

figure(6);
plot(Target(1:200), 'k', 'linewidth', 2);
hold on;
plot(Out(1:200), 'r', 'linewidth',1);
axis([0, 200, -2, 2])
str1 = '\color{black}Target';
str2 = '\color{red}Output';
lg = legend(str1, str2);
set(lg, 'Orientation', 'horizon', 'box', 'off');
ylabel('Prediction')
xlabel('Time (\tau)')
set(gca,'FontName', 'Arial', 'FontSize', 20);
set(gcf, 'unit', 'normalized', 'position', [0.2, 0.2, 0.6, 0.35]);

% 2D map
figure(5);
plot(Target(2:end), 0.3*Target(1:end-1), '.k', 'markersize', 12);
hold on;
plot(Out(2:end), 0.3*Out(1:end-1), '.r', 'markersize', 12);
str1 = '\color{black}Target';
str2 = '\color{red}Output';
lg = legend(str1,str2);
set(lg, 'box', 'off');
ylabel('{\ity} (n)');
xlabel('{\itx} (n)');
axis([-2, 2, -0.4, 0.4]);
set(gca, 'FontName', 'Arial', 'FontSize', 20);
set(gcf, 'unit', 'normalized', 'position', [0.2,0.2,0.3,0.45]);
%% henon预测循环
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
ML = 5;
N = 25;
Vmax = 1.5;
Vmin = -0.5;


dt=0.01;%programming speed%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
step = 4000;
dataset = HenonMap(2*step+500);

Input = dataset(1:step+1);

Target = Input(2:end);

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

Target= dataset(step+1:end);
Input=Target(1:step);
for k=1:200
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
    sprintf('%s', ['loop', num2str(k),'test:', num2str(i), ', Vmax:', num2str(Vmax), ', ML:', num2str(ML)])
end
for i = 1:step
    a = Currents(:, ML*(i-1)+1:ML*i);
    states(:,i) = a(:);
end
states=states*amplifier;
X = [ones(1,step);states];
Out = Wout*X;

NRMSE(k) = sqrt(mean((Out(10:end)-Target(10+k:step+k)).^2)./var(Target(10+k:step+k)));
sprintf('%s',['NRMSE:',num2str(NRMSE)])

for i=1:step-1
Input(i)=Input(i+1);
end
Input(step)=Out(step);
end
figure(1);
plot(Target(10+k:step+k), 'k', 'linewidth', 2);
hold on;
plot(Out(10:end), 'r', 'linewidth',1);
axis([step-500, step, -2, 2])
str1 = '\color{black}Target';
str2 = '\color{red}Output';
lg = legend(str1, str2);
set(lg, 'Orientation', 'horizon', 'box', 'off');
ylabel('Prediction')
xlabel('Time (\tau)')
set(gca,'FontName', 'Arial', 'FontSize', 20);
set(gcf, 'unit', 'normalized', 'position', [0.2, 0.2, 0.6, 0.35]);
hold off
figure(2);
plot(Target(1510+k:step+k), 0.3*Target(1510+k-1:step+k-1), '.k', 'markersize', 12);
hold on;
plot(Out(1510:end), 0.3*Out(1509:end-1), '.r', 'markersize', 12);
str1 = '\color{black}Target';
str2 = '\color{red}Output';
lg = legend(str1,str2);
set(lg, 'box', 'off');
ylabel('{\ity} (n)');
xlabel('{\itx} (n)');
axis([-2, 2, -0.4, 0.4]);
set(gca, 'FontName', 'Arial', 'FontSize', 20);
set(gcf, 'unit', 'normalized', 'position', [0.2,0.2,0.3,0.45]);
hold off

%% 洛伦兹
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
ML = 5;
N = 25;
Vmax = 1.5;
Vmin = -0.5;

dt=0.5;%programming speed%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
step=50001;
time =linspace(0,150,step)
ini=[0.01,0.01,0.01];
[t x]=ode45('chaos',time,ini);
x1=x(:,1)';
y1=x(:,2)';
z1=x(:,3)';
xx=x';
Input=[x1(1:(step-1)/2),y1(1:(step-1)/2),z1(1:(step-1)/2)];
Target=[x1(2:(step+1)/2),y1(2:(step+1)/2),z1(2:(step+1)/2)];


step=length(Input);
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
step=50001;
Input=[x1((step+1)/2:step-1),y1((step+1)/2:step-1),z1((step+1)/2:step-1)];
Target=[x1((step+3)/2:end),y1((step+3)/2:end),z1((step+3)/2:end)];
step=length(Input);
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

step=25000;
figure(1);
plot(Target(1:step), 'k', 'linewidth', 2);
hold on;
plot(Out(1:step), 'r', 'linewidth',1);
axis([0, 200, -2, 2])
str1 = '\color{black}Target';
str2 = '\color{red}Output';
lg = legend(str1, str2);
set(lg, 'Orientation', 'horizon', 'box', 'off');
ylabel('Prediction')
xlabel('Time (\tau)')
set(gca,'FontName', 'Arial', 'FontSize', 20);
set(gcf, 'unit', 'normalized', 'position', [0.2, 0.2, 0.6, 0.35]);


step=25000;
o(:,1)=Target(1:step);
o(:,2)=Target(step+1:step*2);
o(:,3)=Target(step*2+1:end);

p(:,1)=Out(1:step);
p(:,2)=Out(step+1:step*2);
p(:,3)=Out(step*2+1:end);
figure(3)
plot3(o(:,1),o(:,2),o(:,3),'.b', 'markersize', 12);
figure(4)
plot3(p(:,1),p(:,2),p(:,3), '.r', 'markersize', 12);

%% 洛伦兹三个坐标分别计算

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
step=50001;
time =linspace(0,150,step)
ini=[0.01,0.01,0.01];
[t x]=ode45('chaos',time,ini);
x1=x(:,1)';
y1=x(:,2)';
z1=x(:,3)';
xx=x';
Input1=x1(1:(step-1)/2);
Target1=x1(2:(step+1)/2);
Input=Input1;
Target=Target1;

step=length(Input);
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
step=50001;
Input1=x1((step+1)/2:step-1);
Target1=x1((step+3)/2:end);
Input=Input1;
Target=Target1;
step=length(Input);
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
Out1=Out;
NRMSE = sqrt(mean((Out1(10:end)-Target1(10:end)).^2)./var(Target1(10:end)));
sprintf('%s',['NRMSE:',num2str(NRMSE)])

Input2=y1(1:(step-1)/2);
Target2=y1(2:(step+1)/2);
Input=Input2;
Target=Target2;


step=length(Input);
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
step=50001;
Input2=y1((step+1)/2:step-1);
Target2=y1((step+3)/2:end);
Input=Input2;
Target=Target2;
step=length(Input);
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
Out2=Out;
NRMSE = sqrt(mean((Out2(10:end)-Target2(10:end)).^2)./var(Target2(10:end)));
sprintf('%s',['NRMSE:',num2str(NRMSE)])
Input3=z1(1:(step-1)/2);
Target3=z1(2:(step+1)/2);
Input=Input3;
Target=Target3;
step=length(Input);
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

step=50001;
Input3=z1((step+1)/2:step-1);
Target3=z1((step+3)/2:end);
Input=Input3;
Target=Target3;
step=length(Input);
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
Out3=Out;
NRMSE = sqrt(mean((Out3(10:end)-Target3(10:end)).^2)./var(Target3(10:end)));
sprintf('%s',['NRMSE:',num2str(NRMSE)])
step=25000;
figure(1);
plot(Target1(1:step), 'k', 'linewidth', 2);
hold on;
plot(Out1(1:step), 'r', 'linewidth',1);
str1 = '\color{black}Target';
str2 = '\color{red}Output';
lg = legend(str1, str2);
ylabel('Prediction x')
xlabel('Time (\tau)')
set(gca,'FontName', 'Arial', 'FontSize', 20);
set(gcf, 'unit', 'normalized', 'position', [0.2, 0.2, 0.6, 0.35]);

figure(11)
plot(Target2(1:step), 'k', 'linewidth', 2);
hold on;
plot(Out2(1:step), 'r', 'linewidth',1);
str1 = '\color{black}Target';
str2 = '\color{red}Output';
lg = legend(str1, str2);
ylabel('Prediction y')
xlabel('Time (\tau)')
set(gca,'FontName', 'Arial', 'FontSize', 20);
set(gcf, 'unit', 'normalized', 'position', [0.2, 0.2, 0.6, 0.35]);

figure(111)
plot(Target3(1:step), 'k', 'linewidth', 2);
hold on;
plot(Out3(1:step), 'r', 'linewidth',1);
str1 = '\color{black}Target';
str2 = '\color{red}Output';
lg = legend(str1, str2);
ylabel('Prediction z')
xlabel('Time (\tau)')
set(gca,'FontName', 'Arial', 'FontSize', 20);
set(gcf, 'unit', 'normalized', 'position', [0.2, 0.2, 0.6, 0.35]);


step=25000;
o(:,1)=Target1;
o(:,2)=Target2;
o(:,3)=Target3;

p(:,1)=Out1;
p(:,2)=Out2;
p(:,3)=Out3;
figure(3)
plot3(o(:,1),o(:,2),o(:,3),'.b', 'markersize', 12);
figure(4)
plot3(p(:,1),p(:,2),p(:,3), '.r', 'markersize', 12);