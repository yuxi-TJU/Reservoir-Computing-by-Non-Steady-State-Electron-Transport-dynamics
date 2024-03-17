function [result1] = hoppingchannel(eg1,lambda1,gr,gl,VV )
q=1.602e-19;
kB=0.000086;

T=300; %�¶�
k=10^8; %�ڲ����ת������,���ù�
N=1;%λ����
alphal=0.5;%����ѹ
alphar=1-alphal;%�Ҳ��ѹ



NE=501;
E=linspace(-5,5,NE);
dE=E(2)-E(1);

for i= 1: length(VV)
%����������̬��ԭ̬�ֲ�
    FCrab=(1/(4*pi*lambda1*kB*T)^0.5).*exp(-((E+0.5*VV(i))-(eg1-(alphar)*VV(i))-lambda1).^2./(4*lambda1*kB*T));
    FCrba=(1/(4*pi*lambda1*kB*T)^0.5).*exp(-(-(E+0.5*VV(i))+(eg1-(alphar)*VV(i))-lambda1).^2./(4*lambda1*kB*T));
    FClab=(1/(4*pi*lambda1*kB*T)^0.5).*exp(-((E-0.5*VV(i))-(eg1+(alphal)*VV(i))-lambda1).^2./(4*lambda1*kB*T));
    FClba=(1/(4*pi*lambda1*kB*T)^0.5).*exp(-(-(E-0.5*VV(i))+(eg1+(alphal)*VV(i))-lambda1).^2./(4*lambda1*kB*T));
    %kf=k*exp((alpham*Vh(iV))/(2*(N-1)*kB*T));
    %kb=k*exp(-(alpham*Vh(iV))/(2*(N-1)*kB*T));
    Vd=VV(i);
    UL=(0.5*Vd);
    flab=1./(1+exp((E-UL)./(kB*T)));
    flba=1./(1+exp(-(E-UL)./(kB*T)));
    frba=1./(1+exp(-(E+UL)./(kB*T)));
    frab=1./(1+exp((E+UL)./(kB*T)));
    %Rate integration
    kabl1=dE*sum(flab.*FClab)*gl;
    kbal1=dE*(sum(flba.*FClba))*gl;
    kabr1=dE*sum((frab.*FCrab))*gr;
    kbar1=dE*(sum(frba.*FCrba))*gr;
    %%%������������ʷ��̾���
    %%%��������ڵ��ת������
 
    K(1,1)=-(kbal1+kbar1);
    K(1,2)=kabl1+kabr1;
    K(2,1)=1;
    K(2,2)=1;
        
    B(N+1,1)=1;
    %%%%��ⲻͬλ�����ռ�ݸ���
    P=K\B;
     Current=-q*((P(2)*kabl1)-(kbal1*P(1)));
%           Current=-10^8*q*((P(2)*kabl1)-(kbal1*P(1)));
IV=[VV(i),Current];
result1(i,:)=[IV,P'];
end
end

