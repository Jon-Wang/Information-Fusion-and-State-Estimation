%*************************************%
%sensor-3
%*************************************%
path(path,'C:\Documents and Settings\Administrator\桌面\现代时间序列分析方法')
clc;close all;clear;Bushu=100;
%*************************************%
T=1.5;
Qw=2;
Qv3=[3 0;0 1];
fai=[1 T;0 1];gama=[0.5*T^2;T];H03=[1 0;0 1];H13=[1 0;0 1];
%--------------------%
randn('seed',4);w=sqrt(Qw)*randn(1,Bushu+10);
randn('seed',5);v3=sqrt(Qv3)*randn(2,Bushu+10);
x(:,2)=[0;0];x(:,1)=[0;0];
y3(:,1)=v3(:,1);y3(:,2)=H03*x(:,2)+H13*x(:,1)+v3(:,2);
for i=3:Bushu+10
    x(:,i)=fai*x(:,i-1)+gama*w(i-1);
    y3(:,i)=H03*x(:,i)+H13*x(:,i-1)+v3(:,i);
end
%--------------------%
A0=eye(2);
A1=[-1 -T;0 -1];
B0=[0;0];
B1=[0.5*T*T;T];
B2=[0.5*T*T;T];
%--------------------%
R(1:2,1:2)=B0*Qw*B0'+B1*Qw*B1'+B2*Qw*B2'+A0*Qv3*A0'+A1*Qv3*A1';
R(:,3:4)=B1*Qw*B0'+B2*Qw*B1'+A1*Qv3*A0';
% R(:,5:6)=B2*Qw*B0';    %R(:,5:6)=0
%*************************************%
[Qe,d]=TS_GW(Bushu,R);
%*************************************%
D1=d(1:2,199:200);Qe=Qe(1:2,199:200);
%----------白噪声估计----------%
e3(:,1)=[0;0];e3(:,2)=[0;0];
for t=3:Bushu+2
    e3(:,t)=A0*y3(:,t)+A1*y3(:,t-1)-D1*e3(:,t-1);
end
for t=3:Bushu+1
    v3jian(:,t)=Qv3*A0'*inv(Qe)*e3(:,t)+Qv3*(-D1*A0+A1)'*inv(Qe)*e3(:,t+1);
    w3jian(t)=Qw*B0'*inv(Qe)*e3(:,t);
end
%----------误差方差阵----------%
P3(:,:,1)=0.1*eye(2);
for i=1:Bushu
    PP3(:,:,i+1)=fai* P3(:,:,i)*fai'+gama*Qw*gama';%预报误差方差阵
    Qeps3(:,:,i+1)=(H03+H13*inv(fai))*PP3(:,:,i+1)*(H03+H13*inv(fai))'+(H13*inv(fai)*gama)*Qw*(H13*inv(fai)*gama)'+Qv3;
    k3(:,:,i+1)=(PP3(:,:,i+1)*(H03+H13*inv(fai))')*inv(Qeps3(:,:,i+1));
    P3(:,:,i+1)=(fai-k3(:,:,i+1)*H03*fai-k3(:,:,i+1)*H13)*P3(:,:,i)*(fai-k3(:,:,i+1)*H03*fai-k3(:,:,i+1)*H13)'+...
                k3(:,:,i+1)*Qv3*k3(:,:,i+1)'+((eye(2)-k3(:,:,i+1)*H03)*gama)*Qw*((eye(2)-k3(:,:,i+1)*H03)*gama)';%滤波误差方差阵
end 
%------------------------------%
for t=3:Bushu+1
    x3jian(:,t)=inv(H03*fai+H13)*y3(:,t)-inv(H03*fai+H13)*v3jian(:,t);
end
%------------------------------%

t=1:Bushu;
figure;
subplot(2,2,1);plot(t,v3jian(1,t),'b:');
subplot(2,2,2);plot(t,v3jian(2,t),'b:');
subplot(2,2,3);plot(t,w3jian(t),'b:');
figure;
subplot(2,2,1);plot(t,x(1,t),'b',t,x3jian(1,t),'r:');
subplot(2,2,2);plot(t,x(2,t),'b',t,x3jian(2,t),'r:');
