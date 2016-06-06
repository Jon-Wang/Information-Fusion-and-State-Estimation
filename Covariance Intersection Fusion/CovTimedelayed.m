%-----------------------------------------%
%多传感器多重观测滞后CI融合
%-----------------------------------------%
clc;clear all;
Bushu=300;
T=1.5;
Qw=2;
Rv1=1;Rv2=0.54;Rv3=2;
% N=200;
% for j=1:N
randn('seed',1);w=sqrt(Qw)*randn(1,Bushu+10);
randn('seed',1);v1=sqrt(Rv1)*randn(1,Bushu+10);
randn('seed',1);v2=sqrt(Rv2)*randn(1,Bushu+10);
randn('seed',5);v3=sqrt(Rv3)*randn(1,Bushu+10);
%--------------给出系统模型---------------%
fai=[1 T;0 1];gama=[0.5*T^2;T];H01=[1 0];H12=[0.4 0.8];H03=[1.5 0];H13=[0.6 1];
% fai=[1 1.2;0.3 1];gama=[1;3];H01=[1 0];H12=[0.4 1];H03=[1.5 0];H13=[0.6 1];
x(:,2)=[0;0];x(:,1)=[0;0];
y1(2)=H01*x(:,2)+v1(2);y2(2)=H12*x(:,1)+v2(2);y3(2)=H03*x(:,2)+H13*x(:,1)+v3(2);
for i=3:Bushu+10
    x(:,i)=fai*x(:,i-1)+gama*w(i-1);
    y1(i)=H01*x(:,i)+v1(i);y2(i)=H12*x(:,i-1)+v2(i);y3(:,i)=H03*x(:,i)+H13*x(:,i-1)+v3(i);
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

%-----------------传感器2----------------%
PP2(:,:,1)=0.1*eye(2);P2(:,:,1)=0.1*eye(2);x2jian(:,1)=zeros(2,1);
for i=1:Bushu+5
    PP2(:,:,i+1)=fai* P2(:,:,i)*fai'+gama*Qw*gama';%预报误差方差阵
    Qeps2(i+1)=H12*P2(:,:,i)*H12'+Rv2;
    k2(:,i+1)=(fai*P2(:,:,i)*H12')*inv(Qeps2(i+1));
    P2(:,:,i+1)=PP2(:,:,i+1)-k2(:,i+1)*Qeps2(i+1)*k2(:,i+1)';%滤波误差方差阵
    x2jianp(:,i+1)=fai*x2jian(:,i);%预报
    eps2(i+1)=y2(i+1)-H12*x2jian(:,i);
    x2jian(:,i+1)=x2jianp(:,i+1)+k2(:,i+1)*eps2(i+1);%滤波
end 
% t=1:Bushu;
% figure
% subplot(2,2,1);plot(t,x(1,t),'b',t,x2jian(1,t),'r:');
% subplot(2,2,2);plot(t,x(2,t),'b',t,x2jian(2,t),'r:');

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
    b(i)=trace(P2(:,:,i));
    c(i)=trace(P3(:,:,i));
end

%-----------------滤波误差互协方差阵----------------%
P12(:,:,1)=eye(2);P13(:,:,1)=eye(2);P23(:,:,1)=eye(2);
for i=1:Bushu  
    PP12(:,:,i+1)=fai* P12(:,:,i)*fai'+gama*Qw*gama';%预报误差互协方差阵
    PP13(:,:,i+1)=fai* P13(:,:,i)*fai'+gama*Qw*gama';
    PP23(:,:,i+1)=fai* P23(:,:,i)*fai'+gama*Qw*gama';
    %---滤波误差互协方差阵----%
    P12(:,:,i+1)=[eye(2)-k1(:,i+1)*H01]*PP12(:,:,i+1)-[eye(2)-k1(:,i+1)*H01]*fai*P12(:,:,i)*H12'*k2(:,i+1)';P21(:,:,i+1)=P12(:,:,i+1)';
    P13(:,:,i+1)=[eye(2)-k1(:,i+1)*H01]*PP13(:,:,i+1)*[eye(2)-k3(:,i+1)*H03]'-[eye(2)-k1(:,i+1)*H01]*fai*P13(:,:,i)*H13'*k3(:,i+1)';P31(:,:,i+1)=P13(:,:,i+1)';
    P23(:,:,i+1)=PP23(:,:,i+1)*[eye(2)-k3(:,i+1)*H03]'+k1(:,i+1)*H12*P23(:,:,i)*H13'*k3(:,i+1)'-fai*P23(:,:,i)*H13'*k3(:,i+1)'-k2(:,i+1)*H12*P23(:,:,i)*fai'*[eye(2)-k3(:,i+1)*H03];P32(:,:,i+1)=P23(:,:,i+1)';
end

%-----------------按矩阵加权---------------%
for i=1:Bushu
    Psigma(:,:,i)=[P1(:,:,i),P12(:,:,i),P13(:,:,i);
                  P12(:,:,i)',P2(:,:,i),P23(:,:,i);
                  P13(:,:,i)',P23(:,:,i)',P3(:,:,i)];
end
e=[eye(2),eye(2),eye(2)]';
for i=1:Bushu
    A(:,:,i)=inv(Psigma(:,:,i))*e*inv(e'*inv(Psigma(:,:,i))*e);
    Pf(:,:,i)=inv(e'*inv(Psigma(:,:,i))*e);%误差方差阵
    xfjian(:,i)=A(1:2,:,i)'*x1jian(:,i)+A(3:4,:,i)'*x2jian(:,i)+A(5:6,:,i)'*x3jian(:,i);
end
% t=1:Bushu;
% figure
% subplot(2,2,1);plot(t,x(1,t),'b',t,xfjian(1,t),'r:');
% subplot(2,2,2);plot(t,x(2,t),'b',t,xfjian(2,t),'r:');
%-----------------exf---------------%
exf(300)=[0];
exf(1)=(xfjian(1,1)-x(1,1))*(xfjian(1,1)-x(1,1));
for i=2:Bushu
   exf(i)=(xfjian(1,i)-x(1,i))*(xfjian(1,i)-x(1,i))+exf(i-1);
end   
% t=1:Bushu;
% figure
% subplot(2,2,1);plot(t,exf(1,t),'b');
% subplot(2,2,2);plot(t,exf(2,t),'b');

%-----------------BCI----------------%
f=f_bcifun(w);
Aeq=[1 1 1];Beq=[1]; LB=[0 0 0]';UB=[1 1 1]';A=[];B=[];
w0=[0.4;0.4;0.2];
w=fmincon(@f_bcifun,w0,A,B,Aeq,Beq,LB,UB)
PBci=inv(w(1)*inv(P1(:,:,Bushu))+w(2)*inv(P2(:,:,Bushu))+w(3)*inv(P3(:,:,Bushu)));%CI融合估值误差方差阵
JBci=trace(PBci);%CI融合估值误差方差阵的迹
for i=1:Bushu
    xBci(:,i)=PBci*(w(1)*inv(P1(:,:,Bushu))*x1jian(:,i)+w(2)*inv(P2(:,:,Bushu))*x2jian(:,i)+w(3)*inv(P3(:,:,Bushu))*x3jian(:,i));
end
% %-----------------BCI----------------%
% % PBci_=PBci*(w(1)*w(1)*inv(p11(:,:,Bushu))*p11(:,:,Bushu)*inv(p11(:,:,Bushu))+w(1)*w(2)*inv(p11(:,:,Bushu))*p12(:,:,Bushu)*inv(p22(:,:,Bushu))+...
% %     w(1)*w(3)*inv(p11(:,:,Bushu))*p13(:,:,Bushu)*inv(p33(:,:,Bushu))+...
% %     w(2)*w(1)*inv(p22(:,:,Bushu))*p12(:,:,Bushu)'*inv(p11(:,:,Bushu))+w(2)*w(2)*inv(p22(:,:,Bushu))*p22(:,:,Bushu)*inv(p22(:,:,Bushu))+...
% %     w(2)*w(3)*inv(p22(:,:,Bushu))*p23(:,:,Bushu)*inv(p33(:,:,Bushu))+...        
% %     w(3)*w(1)*inv(p33(:,:,Bushu))*p13(:,:,Bushu)'*inv(p11(:,:,Bushu))+w(3)*w(2)*inv(p33(:,:,Bushu))*p23(:,:,Bushu)'*inv(p22(:,:,Bushu))+...
% %     w(3)*w(3)*inv(p33(:,:,Bushu))*p33(:,:,Bushu)*inv(p33(:,:,Bushu)))*PBci;
% % JBci_=trace(PBci_);
% P1_ni=inv(P1(:,:,Bushu));P2_ni=inv(P2(:,:,Bushu));P3_ni=inv(P3(:,:,Bushu));
% PBci_ni=inv(PBci);%PBci_ni_=inv(PBci_);
% theta=0:pi/100:2*pi;
%-----------------SCI----------------%
deta=0.0001;pp1=P1(:,:,Bushu);pp2=P2(:,:,Bushu);pp3=P3(:,:,Bushu);
% [w12,Pci12]=y12_0618(pp1,pp2,deta)%P1和P2形成pci1
% [w123,Pci123,trPci123]=y123_0618(Pci12,pp3,deta)%P3和pci1融合成psci123
[w13,Pci13]=y13_0618(deta,pp1,pp3)%pp1和pp3形成pci1
[w132,Pci132,trPci132]=y132_0618(deta,Pci13,pp2)%pp2和pci1融合成psci132
for i=1:Bushu
%      xci12(:,i)=Pci12*(w12*inv(P1(:,:,Bushu))*x1jian(:,i)+(1-w12)*inv(P2(:,:,Bushu))*x2jian(:,i));
%      xci123(:,i)=Pci123*(w123*inv(Pci12)*xci12(:,i)+(1-w123)*inv(P3(:,:,Bushu))*x3jian(:,i));%xci123
     xci13(:,i)=Pci13*(w13*inv(pp1)*x1jian(:,i)+(1-w13)*inv(pp3)*x3jian(:,i));
     xci132(:,i)=Pci132*(w132*inv(Pci13)*xci13(:,i)+(1-w132)*inv(pp2)*x2jian(:,i));%xci132
end
t=1:Bushu;
figure
plot(t,x(1,t),'b',t,xci132(1,t),'r:');
figure
plot(t,x(2,t),'b',t,xci132(2,t),'r:');axis([0,Bushu,-45,20]);
%-----------------exci---------------%
exci(300)=[0];
exci(1)=(xci132(1,1)-x(1,1))*(xci132(1,1)-x(1,1));
for i=2:Bushu
   exci(i)=(xci132(1,i)-x(1,i))*(xci132(1,i)-x(1,i))+exci(i-1);
end   
t=1:Bushu;
figure
plot(t,exf(t),'b',t,exci(t),'r');axis([0,300,0,180]);

Pci132_=Pci132*(w132*w13*w132*w13*inv(pp1)*pp1*inv(pp1)+w132*w13*w132*(1-w13)*inv(pp1)*P13(:,:,Bushu)*inv(pp3)+...
        w132*w13*(1-w132)*inv(pp1)*P12(:,:,Bushu)*inv(pp2)+w132*(1-w13)*w132*w13*inv(pp3)*P31(:,:,Bushu)*inv(pp1)+w132*(1-w13)*w132*(1-w13)*inv(pp3)+...
        w132*(1-w13)*(1-w132)*inv(pp3)*P32(:,:,Bushu)*inv(pp2)+(1-w132)*w132*w13*inv(pp2)*P21(:,:,Bushu)*inv(pp1)+(1-w132)*w132*(1-w13)*inv(pp2)*P23(:,:,Bushu)*inv(pp3)+...
        (1-w132)*(1-w132)*inv(pp2)*pp2*inv(pp2))*Pci132;
%-----------------椭圆半径----------------%
P1_ni=inv(P1(:,:,Bushu));P2_ni=inv(P2(:,:,Bushu));P3_ni=inv(P3(:,:,Bushu));
% Pci123_ni=inv(Pci123);
Pci132_ni=inv(Pci132);Pci132__ni=inv(Pci132_);
% PBci_ni=inv(PBci);%PBci_ni_=inv(PBci_);
theta=0:pi/100:2*pi;
r1=1./sqrt(P1_ni(1,1)*cos(theta).^2+(P1_ni(1,2)+P1_ni(2,1))*cos(theta).*sin(theta)+P1_ni(2,2)*sin(theta).^2);
r2=1./sqrt(P2_ni(1,1)*cos(theta).^2+(P2_ni(1,2)+P2_ni(2,1))*cos(theta).*sin(theta)+P2_ni(2,2)*sin(theta).^2);
r3=1./sqrt(P3_ni(1,1)*cos(theta).^2+(P3_ni(1,2)+P3_ni(2,1))*cos(theta).*sin(theta)+P3_ni(2,2)*sin(theta).^2);
% rBci=1./sqrt(PBci_ni(1,1)*cos(theta).^2+(PBci_ni(1,2)+PBci_ni(2,1))*cos(theta).*sin(theta)+PBci_ni(2,2)*sin(theta).^2);
% rBci_=1./sqrt(PBci_ni_(1,1)*cos(theta).^2+(PBci_ni_(1,2)+PBci_ni_(2,1))*cos(theta).*sin(theta)+PBci_ni_(2,2)*sin(theta).^2);
% rci123=1./sqrt(Pci123_ni(1,1)*cos(theta).^2+(Pci123_ni(1,2)+Pci123_ni(2,1))*cos(theta).*sin(theta)+Pci123_ni(2,2)*sin(theta).^2);
rci132=1./sqrt(Pci132_ni(1,1)*cos(theta).^2+(Pci132_ni(1,2)+Pci132_ni(2,1))*cos(theta).*sin(theta)+Pci132_ni(2,2)*sin(theta).^2);
rci132_=1./sqrt(Pci132__ni(1,1)*cos(theta).^2+(Pci132__ni(1,2)+Pci132__ni(2,1))*cos(theta).*sin(theta)+Pci132__ni(2,2)*sin(theta).^2);
%-----------------MSE曲线----------------%
%     for i=1:Bushu
%         ErroP1(j,i)=(x1jian(:,i)-x(:,i))'*(x1jian(:,i)-x(:,i));
%         ErroP2(j,i)=(x2jian(:,i)-x(:,i))'*(x2jian(:,i)-x(:,i));
%         ErroP3(j,i)=(x3jian(:,i)-x(:,i))'*(x3jian(:,i)-x(:,i));
%         ErroPc(j,i)=(xcjian(:,i)-x(:,i))'*(xcjian(:,i)-x(:,i));
%         ErroPm(j,i)=(xmjian(:,i)-x(:,i))'*(xmjian(:,i)-x(:,i));
%         ErroPd(j,i)=(xdjian(:,i)-x(:,i))'*(xdjian(:,i)-x(:,i));
%         ErroPs(j,i)=(xsjian(:,i)-x(:,i))'*(xsjian(:,i)-x(:,i));
%         ErroPBci(j,i)=(xBci(:,i)-x(:,i))'*(xBci(:,i)-x(:,i));
%     end
%     %  end
%      MSE1=sum(ErroP1)/N;
%      MSE2=sum(ErroP2)/N;
%      MSE3=sum(ErroP3)/N;
%      MSEc=sum(ErroPc)/N;
%      MSEm=sum(ErroPm)/N;
%      MSEd=sum(ErroPd)/N;
%      MSEs=sum(ErroPs)/N;
%      MSEci=sum(ErroPBci)/N;
%--------------作图----------------%
t=1:Bushu;
figure 
hold on;
polar(theta,r1,'b');
polar(theta,r2,'k');
polar(theta,r3,'r');
polar(theta,rci132,'k');
polar(theta,rci132_,'k-.');axis([-2.3,2.3,-2.3,2.3]);
% polar(theta,rci123,'b-.');
% polar(theta,rBci,'b-.');
%polar(theta,rBci_,'b-.');%  figure
%  t=50:50:Bushu;
% plot(t,MSE1(t),'r-h',t,MSE2(t),'r-s',t,MSE3(t),'r-x',t,MSEc(t),'k-^',t,MSEm(t),'g-o',t,MSEd(t),'g-d',t,MSEs(t),'g-*',t,MSEci(t),'r-v'); 
% legend('MSE1','MSE2','MSE3','MSEc','MSEm','MSEd','MSEs','MSEci');
% hold on
% t=1:Bushu;
% line([50,Bushu],[a(Bushu),a(Bushu)]);line([50,Bushu],[b(Bushu),b(Bushu)]);line([50,Bushu],[c(Bushu),c(Bushu)]);
% line([50,Bushu],[jiPm(Bushu),jiPm(Bushu)]);line([50,Bushu],[jiPd(Bushu),jiPd(Bushu)]);line([50,Bushu],[jiPs(Bushu),jiPs(Bushu)]);
% line([50,Bushu],[JBci,JBci]);line([50,Bushu],[JBci_,JBci_]);
