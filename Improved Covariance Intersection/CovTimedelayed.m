%-----------------------------------------%
%改进的CI融合估值器
%20160418
%-----------------------------------------%
clc;clear all;
Bushu=300;
T=1.5;
Qw=2;
Rv1=1;Rv2=0.54;Rv3=0.8;

randn('seed',1);w=sqrt(Qw)*randn(1,Bushu+10);
randn('seed',1);v1=sqrt(Rv1)*randn(1,Bushu+10);

randn('seed',5);v3=sqrt(Rv3)*randn(1,Bushu+10);
%--------------给出系统模型---------------%
fai=[1 T;0 1];gama=[0.5*T^2;T];H01=[1 0.5];H03=[1.5 0];H13=[0.6 1];
x(:,1)=[0;0];x(:,2)=[0;0];
y1(2)=H01*x(:,2)+v1(2);y3(2)=H03*x(:,2)+H13*x(:,1)+v3(2);
for i=3:Bushu+10
    x(:,i)=fai*x(:,i-1)+gama*w(i-1);
    y1(i)=H01*x(:,i)+v1(i);y3(:,i)=H03*x(:,i)+H13*x(:,i-1)+v3(i);
end
%------------------局部Kalman滤波器---------------%
%---------------------传感器1----------------%
PP1(:,:,1)=0.1*eye(2);P1(:,:,1)=0.1*eye(2);x1jian(:,1)=zeros(2,1);
for i=1:Bushu+5
    PP1(:,:,i+1)=fai* P1(:,:,i)*fai'+gama*Qw*gama';%预报误差方差阵
    Qeps1(i+1)=H01*PP1(:,:,i+1)*H01'+Rv1;
    k1(:,i+1)=(PP1(:,:,i+1)*H01')*inv(Qeps1(i+1));
    P1(:,:,i+1)=PP1(:,:,i+1)-k1(:,i+1)*Qeps1(i+1)*k1(:,i+1)';%滤波误差方差阵
    x1jianp(:,i+1)=fai*x1jian(:,i);%预报
    eps1(i+1)=y1(i+1)-H01*x1jianp(:,i+1);
    x1jian(:,i+1)=x1jianp(:,i+1)+k1(:,i+1)*eps1(i+1);%滤波
end 
% t=1:Bushu;
% figure
% subplot(2,2,1);plot(t,x(1,t),'b',t,x1jian(1,t),'r:');
% subplot(2,2,2);plot(t,x(2,t),'b',t,x1jian(2,t),'r:');

%-----------------传感器3----------------%
PP3(:,:,1)=0.1*eye(2);P3(:,:,1)=0.1*eye(2);x3jian(:,1)=zeros(2,1);
for i=1:Bushu+5
    PP3(:,:,i+1)=fai* P3(:,:,i)*fai'+gama*Qw*gama';%预报误差方差阵
    Qeps3(i+1)=H03*PP3(:,:,i+1)*H03'+H13*P3(:,:,i)*H13'+H03*fai*P3(:,:,i)*H13'+H13*P3(:,:,i)*fai'*H03'+Rv3;
    k3(:,i+1)=(PP3(:,:,i+1)*H03'+fai*P3(:,:,i)*H13')*inv(Qeps3(i+1));
    P3(:,:,i+1)=PP3(:,:,i+1)-k3(:,i+1)*Qeps3(i+1)*k3(:,i+1)';%滤波误差方差阵
    x3jianp(:,i+1)=fai*x3jian(:,i);%预报
    eps3(i+1)=y3(i+1)-H03*x3jianp(:,i+1)-H13*x3jian(:,i);
    x3jian(:,i+1)=x3jianp(:,i+1)+k3(:,i+1)*eps3(i+1);%滤波
end
% t=1:Bushu;
% figure
% subplot(2,2,1);plot(t,x(1,t),'b',t,x3jian(1,t),'r:');
% subplot(2,2,2);plot(t,x(2,t),'b',t,x3jian(2,t),'r:');

%-----------------迹----------------%
for i=1:Bushu+5
    a(i)=trace(P1(:,:,i));
    c(i)=trace(P3(:,:,i));
end

%-----------------滤波误差互协方差阵----------------%
P13(:,:,1)=eye(2);
for i=1:Bushu  
    PP13(:,:,i+1)=fai* P13(:,:,i)*fai'+gama*Qw*gama';%预报误差互协方差阵
    %---滤波误差互协方差阵----%
    P13(:,:,i+1)=[eye(2)-k1(:,i+1)*H01]*PP13(:,:,i+1)*[eye(2)-k3(:,i+1)*H03]'-[eye(2)-k1(:,i+1)*H01]*fai*P13(:,:,i)*H13'*k3(:,i+1)';P31(:,:,i+1)=P13(:,:,i+1)';
end
%-----------------按矩阵加权---------------%
for i=1:Bushu
    Psigma(:,:,i)=[P1(:,:,i),P13(:,:,i);
                  P13(:,:,i)',P3(:,:,i)];
end
e=[eye(2),eye(2)]';
for i=1:Bushu
    A(:,:,i)=inv(Psigma(:,:,i))*e*inv(e'*inv(Psigma(:,:,i))*e);
    Pm(:,:,i)=inv(e'*inv(Psigma(:,:,i))*e);%误差方差阵
    xmjian(:,i)=A(1:2,:,i)'*x1jian(:,i)+A(3:4,:,i)'*x3jian(:,i);
end
% t=1:Bushu;
% figure
% subplot(2,2,1);plot(t,x(1,t),'b',t,xmjian(1,t),'r:');
% subplot(2,2,2);plot(t,x(2,t),'b',t,xmjian(2,t),'r:');
%-----------------SCI----------------%
deta=0.0001;pp1=P1(:,:,Bushu);pp3=P3(:,:,Bushu);
[w13,Pci13]=y13_0618(deta,pp1,pp3)%pp1和pp3形成pci1

Pci13_=Pci13*(w13*w13*inv(pp1)*pp1*inv(pp1)+w13*(1-w13)*inv(pp1)*P13(:,:,Bushu)*inv(pp3)+...
       w13*(1-w13)*inv(pp3)*P31(:,:,Bushu)*inv(pp1)+(1-w13)*(1-w13)*inv(pp3)*pp3*inv(pp3))*Pci13;
 
% P_13=[0.15 0.15;0.3 1];P_31=P_13';
rho=0.7;
P_13=rho*chol(pp1)'*chol(pp3);P_31=P_13';
Pcic13=Pci13*(w13*w13*inv(pp1)*pp1*inv(pp1)+w13*(1-w13)*inv(pp1)*P_13*inv(pp3)+...
       w13*(1-w13)*inv(pp3)*P_31*inv(pp1)+(1-w13)*(1-w13)*inv(pp3)*pp3*inv(pp3))*Pci13;
   
for i=1:Bushu
     xci13(:,i)=Pci13*(w13*inv(pp1)*x1jian(:,i)+(1-w13)*inv(pp3)*x3jian(:,i));
end
t=1:Bushu;
figure
subplot(2,2,1);plot(t,x(1,t),'b',t,xci13(1,t),'r:');
subplot(2,2,2);plot(t,x(2,t),'b',t,xci13(2,t),'r:');axis([0,Bushu,-45,20]);
%-----------------椭圆半径----------------%
P1_ni=inv(P1(:,:,Bushu));P3_ni=inv(P3(:,:,Bushu));
Pm_ni=inv(Pm(:,:,Bushu));
Pci13_ni=inv(Pci13);Pci13__ni=inv(Pci13_);
Pcic13_ni=inv(Pcic13);
theta=0:pi/100:2*pi;
r1=1./sqrt(P1_ni(1,1)*cos(theta).^2+(P1_ni(1,2)+P1_ni(2,1))*cos(theta).*sin(theta)+P1_ni(2,2)*sin(theta).^2);
r3=1./sqrt(P3_ni(1,1)*cos(theta).^2+(P3_ni(1,2)+P3_ni(2,1))*cos(theta).*sin(theta)+P3_ni(2,2)*sin(theta).^2);

rm=1./sqrt(Pm_ni(1,1)*cos(theta).^2+(Pm_ni(1,2)+Pm_ni(2,1))*cos(theta).*sin(theta)+Pm_ni(2,2)*sin(theta).^2);
rci13=1./sqrt(Pci13_ni(1,1)*cos(theta).^2+(Pci13_ni(1,2)+Pci13_ni(2,1))*cos(theta).*sin(theta)+Pci13_ni(2,2)*sin(theta).^2);
rci13_=1./sqrt(Pci13__ni(1,1)*cos(theta).^2+(Pci13__ni(1,2)+Pci13__ni(2,1))*cos(theta).*sin(theta)+Pci13__ni(2,2)*sin(theta).^2);
rcic13=1./sqrt(Pcic13_ni(1,1)*cos(theta).^2+(Pcic13_ni(1,2)+Pcic13_ni(2,1))*cos(theta).*sin(theta)+Pcic13_ni(2,2)*sin(theta).^2);
%--------------作图----------------%
t=1:Bushu;
figure 
hold on;
polar(theta,r1,'b');
polar(theta,r3,'b-.');
polar(theta,rm,'r');
polar(theta,rci13,'k');
polar(theta,rci13_,'k-.');
polar(theta,rcic13,'m');
