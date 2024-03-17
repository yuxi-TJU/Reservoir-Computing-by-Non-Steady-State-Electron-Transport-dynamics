function [amoAQ,amoh2AQ,CurrentsAQ,Currentsh2AQ,Currents] = It3pulse(eg1,lambda1,gr,gl,eg2,lambda2,gr2,gl2,dt,VV,kn,kp,AQ,H2AQ)
%输入的电压序列是每个并行节点的序列，在时间尺度上无关，所以这个函数对应于电压不相关序列，如果是一个连续的电压序列输入进来则无效果。

q=1.602e-19;
kB=0.000086;

result1=hoppingchannel(eg1,lambda1,gr,gl,VV );
result2=hoppingchannel(eg2,lambda2,gr2,gl2,VV );


% dt=1/timestep; %时间变量
% t=0;%初始时间
% tpoint = tall/dt; %循环数
% Timess=linspace(0,tall,tpoint+1);%时间数组

    m1 = result1(:,4);%氧化态
    m2 = result1(:,3);%还原态
    n1 = result2(:,4);%质子化氧化态
    n2 = result2(:,3);%质子化还原态
    
    kn1=kn*n1;
    kp1=kp*m2;
 AQ=(AQ-(kn1./(kp1+kn1))).*exp(-(kp1+kn1)*dt)+(kn1./(kp1+kn1));
%     AQ=(AQ-(kn1/(kp1+kn1))).*exp(-(kp1+kn1).*Timess)+(kn1/(kp1+kn1));
    H2AQ=1-AQ;
    amoAQ = AQ;
    amoh2AQ = H2AQ;
    CurrentsAQ = result1(:,2).*AQ;
    Currentsh2AQ =result2(:,2).*H2AQ;

Currents=CurrentsAQ+Currentsh2AQ;
% amoAQ1=amoAQ';
% amoh2AQ1=amoh2AQ';
% CurrentsAQ1=CurrentsAQ';
% Currentsh2AQ1=Currentsh2AQ';
% Currents1=Currents';

end
