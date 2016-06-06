%*************************************%
%sensor-2
%*************************************%
path(path,'C:\Documents and Settings\Administrator\桌面\现代时间序列分析方法')
clc;close all;clear;Bushu=100;
%*************************************%
T=1.5;
Qw=2;
Qv2=[2 0;0 1];
fai=[1 T;0 1];gama=[0.5*T^2;T];H12=[1 1;0 1];
%--------------------%
randn('seed',4);w=sqrt(Qw)*randn(1,Bushu+10);
randn('seed',1);v2=sqrt(Qv2)*randn(2,Bushu+10);
x(:,2)=[0;0];x(:,1)=[0;0];
y2(:,1)=v2(:,1);y2(:,2)=H12*x(:,1)+v2(:,2);
for i=3:Bushu+10
    x(:,i)=fai*x(:,i-1)+gama*w(i-1);
    y2(:,i)=H12*x(:,i-1)+v2(:,i);
end
%--------------------%
A0=eye(2);
A1=[-1 -T;0 -1];
B0=[0;0];B1=[0;0];
B2=[0.5*T*T+T;T];
%--------------------%
R(1:2,1:2)=B0*Qw*B0'+B1*Qw*B1'+B2*Qw*B2'+A0*Qv2*A0'+A1*Qv2*A1';
R(:,3:4)=B1*Qw*B0'+B2*Qw*B1'+A1*Qv2*A0';
%*************************************%
[Qe,d]=TS_GW(Bushu,R);
%*************************************%
D1=d(1:2,199:200);Qe=Qe(1:2,199:200);
%----------白噪声估计----------%
e2(:,1)=[0;0];e2(:,2)=[0;0];
for t=3:Bushu+2
    e2(:,t)=A0*y2(:,t)+A1*y2(:,t-1)-D1*e2(:,t-1);
end
for t=3:Bushu+1
    v2jian(:,t)=Qv2*A0'*inv(Qe)*e2(:,t)+Qv2*(-D1*A0+A1)'*inv(Qe)*e2(:,t+1);
end
%----------误差方差阵----------%
P2(:,:,1)=0.1*eye(2);
for i=1:Bushu
    PP2(:,:,i+1)=fai* P2(:,:,i)*fai'+gama*Qw*gama';%预报误差方差阵
    Qeps2(:,:,i+1)=(H12*inv(fai))*PP2(:,:,i+1)*(H12*inv(fai))'+(H12*inv(fai)*gama)*Qw*(H12*inv(fai)*gama)'+Qv2;
    k2(:,:,i+1)=(PP2(:,:,i+1)*(inv(fai))'*H12')*inv(Qeps2(:,:,i+1));
    P2(:,:,i+1)=(fai-k2(:,:,i+1)*H12)*P2(:,:,i)*(fai-k2(:,:,i+1)*H12)'+k2(:,:,i+1)*Qv2*k2(:,:,i+1)'+gama*Qw*gama';%滤波误差方差阵
end 
%------------------------------%
for t=3:Bushu+1
    x2jian(:,t)=inv(H12)*y2(:,t)-inv(H12)*v2jian(:,t);
end
t=1:Bushu;
figure;
subplot(2,2,1);plot(t,v2jian(1,t),'b:');
subplot(2,2,2);plot(t,v2jian(2,t),'b:');
subplot(2,2,3);plot(t,x(1,t),'b',t,x2jian(1,t),'r:');
subplot(2,2,4);plot(t,x(2,t),'b',t,x2jian(2,t),'r:');


