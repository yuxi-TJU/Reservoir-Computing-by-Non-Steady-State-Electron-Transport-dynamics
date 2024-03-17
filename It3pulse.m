function [amoAQ,amoh2AQ,CurrentsAQ,Currentsh2AQ,Currents] = It3pulse(eg1,lambda1,gr,gl,eg2,lambda2,gr2,gl2,dt,VV,kn,kp,AQ,H2AQ)
%����ĵ�ѹ������ÿ�����нڵ�����У���ʱ��߶����޹أ��������������Ӧ�ڵ�ѹ��������У������һ�������ĵ�ѹ���������������Ч����

q=1.602e-19;
kB=0.000086;

result1=hoppingchannel(eg1,lambda1,gr,gl,VV );
result2=hoppingchannel(eg2,lambda2,gr2,gl2,VV );


% dt=1/timestep; %ʱ�����
% t=0;%��ʼʱ��
% tpoint = tall/dt; %ѭ����
% Timess=linspace(0,tall,tpoint+1);%ʱ������

    m1 = result1(:,4);%����̬
    m2 = result1(:,3);%��ԭ̬
    n1 = result2(:,4);%���ӻ�����̬
    n2 = result2(:,3);%���ӻ���ԭ̬
    
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
