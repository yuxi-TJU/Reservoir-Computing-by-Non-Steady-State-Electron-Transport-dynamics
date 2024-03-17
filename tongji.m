clear;clc;


A = xlsread('1','聚集')% 数据表\sheet\行列
Bb=[];%zhuanhuan 
Cc=zeros(6,71);%verynice
Dd=zeros(6,71);%shaoliang maoci
Ee=zeros(6,71);%nengkanchulai
Ff=zeros(6,71);%buxingde
m=length(A(:,1));
n=length(A(1,:));

for i=1:m
            a=mod(i,6);
        if a==0
            a=6;
        end
    for j=1:n

        if A(i,j)<0.1
            Bb(i,j)=1;
            Cc(a,j)=Cc(a,j)+1;
        elseif A(i,j)<0.2
            Bb(i,j)=2;
            Dd(a,j)=Dd(a,j)+1;
        elseif A(i,j)<0.3
            Ee(a,j)=Ee(a,j)+1;
             Bb(i,j)=3;
        else 
            Ff(a,j)=Ff(a,j)+1;
            Bb(i,j)=4;
        end
    end
end

%% 其他情况
clear;clc;

A = xlsread('1','聚集')% 数据表\sheet\行列
Bb=[];%zhuanhuan 
Cc=zeros(1,61);%verynice
Dd=zeros(1,61);%shaoliang maoci
Ee=zeros(1,61);%nengkanchulai
Ff=zeros(1,61);%buxingde
m=length(A(:,1));
n=length(A(1,:));

for i=1:m
            a=mod(i,1);
        if a==0
            a=1;
        end
    for j=1:n

        if A(i,j)<0.1
            Bb(i,j)=1;
            Cc(a,j)=Cc(a,j)+1;
        elseif A(i,j)<0.2
            Bb(i,j)=2;
            Dd(a,j)=Dd(a,j)+1;
        elseif A(i,j)<0.3
            Ee(a,j)=Ee(a,j)+1;
             Bb(i,j)=3;
        else 
            Ff(a,j)=Ff(a,j)+1;
            Bb(i,j)=4;
        end
    end
end
%%
%si图用的
clear;clc;


A = xlsread('1','聚集')% 数据表\sheet\行列
Bb=[];%zhuanhuan 
Cc=zeros(5,10);%verynice
Dd=zeros(5,10);%shaoliang maoci
Ee=zeros(5,10);%nengkanchulai
Ff=zeros(5,10);%buxingde
m=length(A(:,1));
n=length(A(1,:));

for i=1:m
            a=mod(i,5);
        if a==0
            a=5;
        end
    for j=1:n

        if A(i,j)<0.12
            Bb(i,j)=1;
            Cc(a,j)=Cc(a,j)+1;
        elseif A(i,j)<0.2
            Bb(i,j)=2;
            Dd(a,j)=Dd(a,j)+1;
        elseif A(i,j)<0.3
            Ee(a,j)=Ee(a,j)+1;
             Bb(i,j)=3;
        else 
            Ff(a,j)=Ff(a,j)+1;
            Bb(i,j)=4;
        end
    end
end
% plot(Cc(1,:)/43)
%%
x=linspace(1,length(Cc(1,:)),length(Cc(1,:)));
n=37;
figure(2);
plot(x,Cc(1,:)/n,x,Cc(2,:)/n,x,Cc(3,:)/n,x,Cc(4,:)/n,x,Cc(5,:)/n, 'k', 'linewidth', 2);

% plot(x,Cc(1,:)/n, 'k', 'linewidth', 2);
lg = legend( '10^-^2 Hz','10^-^1 Hz','10^0 Hz','10^1 Hz','10^2 Hz');
xlabel('parallel devices')
ylabel('recognition rate')
set(gca,'FontName', 'Arial', 'FontSize', 20);
%%
x=linspace(-3,3,61);
n=15;
figure(2);
plot(x,Cc(1,:)/n,x,Cc(2,:)/n,x,Cc(3,:)/n,x,Cc(4,:)/n,x,Cc(5,:)/n,x,Cc(6,:)/n,x,Cc(7,:)/n, 'k', 'linewidth', 2);

% plot(x,Cc(1,:)/n, 'k', 'linewidth', 2);
lg = legend( '0.1','1','5','10','20','50','100');
xlabel('parallel devices')
ylabel('recognition rate')
set(gca,'FontName', 'Arial', 'FontSize', 20);
%%
x=linspace(-3,3,61);
n=1;
figure(4);
plot(x,Cc(1,:)/n, 'k', 'linewidth', 2);


% str1 = '\color{black}Target';
% str2 = '\color{red}Output';
% lg = legend( '10^-1 Hz','10^0 Hz','10^1 Hz','10^2 Hz','10^3 Hz');
% set(lg, 'Orientation', 'horizon');
xlabel('log10(frequency(Hz))')
ylabel('Normalized count')
set(gca,'FontName', 'Arial', 'FontSize', 20);