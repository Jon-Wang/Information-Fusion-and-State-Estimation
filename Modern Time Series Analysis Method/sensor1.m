%*************************************%
%sensor-1
%*************************************%
path(path,'C:\Documents and Settings\Administrator\桌面\现代时间序列分析方法')
clc;close all;clear;Bushu=100;
%*************************************%
T=1.5;
Qw=2;
Qv1=[1 0;0 1];
fai=[1 T;0 1];gama=[0.5*T^2;T];H01=[1 2;0 1];
%--------------------%
randn('seed',4);w=sqrt(Qw)*randn(1,Bushu+10);
randn('seed',1);v1=sqrt(Qv1)*randn(2,Bushu+10);
x(:,2)=[0;0];x(:,1)=[0;0];
y1(:,1)=v1(:,1);y1(:,2)=H01*x(:,2)+v1(:,2);
for i=3:Bushu+10
    x(:,i)=fai*x(:,i-1)+gama*w(i-1);
    y1(:,i)=H01*x(:,i)+v1(:,i);
end
%--------------------%
A0=eye(2);
A1=[-1 -T;0 -1];
B0=[0;0];
B1=[0.5*T*T+2*T;T];
%--------------------%
R(1:2,1:2)=B0*Qw*B0'+B1*Qw*B1'+A0*Qv1*A0'+A1*Qv1*A1';
R(:,3:4)=B1*Qw*B0'+A1*Qv1*A0';
%*************************************%
[Qe,d]=TS_GW(Bushu,R);
%*************************************%
D1=d(1:2,199:200);Qe=Qe(1:2,199:200);
%----------白噪声估计----------%
e(:,1)=[0;0];e(:,2)=[0;0];
for t=3:Bushu+1
    e(:,t)=A0*y1(:,t)+A1*y1(:,t-1)-D1*e(:,t-1);
    v1jian(:,t)=Qv1*A0'*inv(Qe)*e(:,t);
end
%----------误差方差阵----------%
P1(:,:,1)=0.1*eye(2);
for i=1:Bushu
    PP1(:,:,i+1)=fai* P1(:,:,i)*fai'+gama*Qw*gama';%预报误差方差阵
    Qeps1(:,:,i+1)=H01*PP1(:,:,i+1)*H01'+Qv1;
    k1(:,:,i+1)=PP1(:,:,i+1)*H01'*inv(Qeps1(:,:,i+1));
    P1(:,:,i+1)=(eye(2)-k1(:,:,i+1)*H01)*PP1(:,:,i+1);%滤波误差方差阵
end 
%------------------------------%
for t=3:Bushu+1
    x1jian(:,t)=inv(H01)*y1(:,t)-inv(H01)*v1jian(:,t);
end
t=1:Bushu;
figure;
subplot(2,2,1);plot(t,v1jian(1,t),'b:');
subplot(2,2,2);plot(t,v1jian(2,t),'b:');
subplot(2,2,3);plot(t,x(1,t),'b',t,x1jian(1,t),'r:');
subplot(2,2,4);plot(t,x(2,t),'b',t,x1jian(2,t),'r:');

