%-----------------------------------------%
%带相关噪声多传感器多重观测滞后CI融合
%-----------------------------------------%
clc;clear all;
Bushu=100;
% alpha1=2;alpha2=1.5;alpha3=0.5;
% Qxi1=2;Qxi2=0.5;Qxi3=3;Qw=0.03;
alpha1=2;alpha2=1.5;alpha3=0.5;
Qxi1=2;Qxi2=3.5;Qxi3=3;Qw=0.03;
randn('seed',6);w=sqrt(Qw)*randn(1,Bushu+10);
randn('seed',2);xi1=sqrt(Qxi1)*randn(1,Bushu+10);
randn('seed',5);xi2=sqrt(Qxi2)*randn(1,Bushu+10);
randn('seed',4);xi3=sqrt(Qxi3)*randn(1,Bushu+10);
v1=alpha1*w+xi1;v2=alpha2*w+xi2;v3=alpha3*w+xi3;
%-------------------系统模型-------------------%
fai0=[0.8 0;0.1 -0.5];fai1=[-0.4 0;0.6 0.5];fai2=[0.4 -0.6;0 0.5];
gama=[1;-0.5];
% H01=[1 0];H11=[0 2];H21=[1 0.5];
% H02=[0 0];H12=[1 1];H22=[0.4 0.8];
% H03=[1.5 0];H13=[0.6 1];H23=[0 0];
H01=[4 1];H11=[0 2];H21=[1 1];
H02=[5 0.9];H12=[1 1];H22=[0.4 0.8];
H03=[5 0.8];H13=[0.6 1];H23=[0 0];
x(:,3)=[0;0];x(:,2)=[0;0];x(:,1)=[0;0];
y1(1)=H01*x(:,1)+v1(1);y1(2)=H01*x(:,2)+H11*x(:,1)+v1(2);y1(3)=H01*x(:,3)+H11*x(:,2)+H21*x(:,1)+v1(3);
y2(1)=H02*x(:,1)+v2(1);y2(2)=H02*x(:,2)+H12*x(:,1)+v2(2);y2(3)=H02*x(:,3)+H12*x(:,2)+H22*x(:,1)+v2(3);
y3(1)=H03*x(:,1)+v3(1);y3(2)=H03*x(:,2)+H13*x(:,1)+v3(2);y3(3)=H03*x(:,3)+H13*x(:,2)+H23*x(:,1)+v3(3);
for i=4:Bushu+10
    x(:,i)=fai0*x(:,i-1)+fai1*x(:,i-2)+fai2*x(:,i-3)+gama*w(i-1);
    y1(i)=H01*x(:,i)+H11*x(:,i-1)+H21*x(:,i-2)+v1(i);
    y2(i)=H02*x(:,i)+H12*x(:,i-1)+H22*x(:,i-2)+v2(i);
    y3(i)=H03*x(:,i)+H13*x(:,i-1)+H23*x(:,i-2)+v3(i);
end
%--------------------------------------------%
S1=alpha1*Qw;S2=alpha2*Qw;S3=alpha3*Qw;
Qv1=alpha1*alpha1*Qw+Qxi1;Qv2=alpha2*alpha2*Qw+Qxi2;Qv3=alpha3*alpha3*Qw+Qxi3;
J1=gama*S1*inv(Qv1);J2=gama*S2*inv(Qv2);J3=gama*S3*inv(Qv3);
Qw_1=gama*(Qw-S1*inv(Qv1)*S1')*gama';Qw_2=gama*(Qw-S2*inv(Qv2)*S2')*gama';Qw_3=gama*(Qw-S3*inv(Qv3)*S3')*gama';
%------------------局部Kalman滤波器---------------%
%---------------------传感器1----------------%
x1jian(:,1)=zeros(2,1);x1jianp(:,1)=zeros(2,1);x1jianpp(:,1)=zeros(2,1);
p11(:,:,1)=0.1*eye(2);p12(:,:,1)=0.1*eye(2);p15(:,:,1)=0.1*eye(2);p13(:,:,1)=0.1*eye(2);p16(:,:,1)=0.1*eye(2);p19(:,:,1)=0.1*eye(2);
for i=1:Bushu+5
    Qeps1(i+1)=H01*p11(:,:,i)*H01'+H01*p12(:,:,i)*H11'+H01*p13(:,:,i)*H21'+...
               H11*p12(:,:,i)'*H01'+H11*p15(:,:,i)*H11'+H11*p16(:,:,i)*H21'+...
               H21*p13(:,:,i)'*H01'+H21*p16(:,:,i)'*H11'+H21*p19(:,:,i)*H21'+Qv1;
    k1(:,i+1)=p11(:,:,i)*H01'*inv(Qeps1(i+1))+p12(:,:,i)*H11'*inv(Qeps1(i+1))+p13(:,:,i)*H21'*inv(Qeps1(i+1));%滤波增益         
    k1pp(:,i+1)=p12(:,:,i)'*H01'*inv(Qeps1(i+1))+p15(:,:,i)*H11'*inv(Qeps1(i+1))+p16(:,:,i)*H21'*inv(Qeps1(i+1));%平滑增益
    k1pp2(:,i+1)=p13(:,:,i)'*H01'*inv(Qeps1(i+1))+p16(:,:,i)'*H11'*inv(Qeps1(i+1))+p19(:,:,i)*H21'*inv(Qeps1(i+1));%2步平滑增益
    
    p16(:,:,i+1)=p12(:,:,i)-k1(:,i+1)*Qeps1(i+1)*k1pp(:,i+1)';
    p15(:,:,i+1)=p11(:,:,i)-k1(:,i+1)*Qeps1(i+1)*k1(:,i+1)';
    p19(:,:,i+1)=p15(:,:,i)-k1pp(:,i+1)*Qeps1(i+1)*k1pp(:,i+1)';
    p111(:,:,i+1)=p13(:,:,i)'- k1pp2(:,i+1)*Qeps1(i+1)*k1(:,i+1)';
    p112(:,:,i+1)=p16(:,:,i)'- k1pp2(:,i+1)*Qeps1(i+1)*k1pp(:,i+1)';
    p113(:,:,i+1)=p19(:,:,i)- k1pp2(:,i+1)*Qeps1(i+1)*k1pp2(:,i+1)';
    p110(:,:,i+1)=p111(:,:,i+1)*(fai0-J1*H01)+p112(:,:,i+1)*(fai1-J1*H11)+p113(:,:,i+1)*(fai2-J1*H21);
    p114(:,:,i+1)=p13(:,:,i)-k1(:,i+1)*Qeps1(i+1)*k1pp2(:,i+1)';
    p12(:,:,i+1)=(p15(:,:,i+1)*(fai0-J1*H01)'+p16(:,:,i+1)*(fai1-J1*H11)'+p114(:,:,i+1)*(fai2-J1*H21))';
    p115(:,:,i+1)=p16(:,:,i)-k1pp(:,i+1)*Qeps1(i+1)*k1pp2(:,i+1)';
    p13(:,:,i+1)=(p12(:,:,i+1)*(fai0-J1*H01)'+p19(:,:,i+1)*(fai1-J1*H11)'+p115(:,:,i+1)*(fai2-J1*H21))';
    p11(:,:,i+1)=(fai0-J1*H01)*p12(:,:,i+1)'+(fai1-J1*H11)*p13(:,:,i+1)'+(fai2-J1*H21)*p110(:,:,i+1)+Qw_1;
    
    eps1(i+1)=y1(i+1)-H01*x1jianp(:,i)-H11*x1jian(:,i)-H21*x1jianpp(:,i);
    x1jianpp2(:,i+1)=x1jianpp(:,i)+k1pp2(:,i+1)*eps1(i+1);%2步平滑
    x1jianpp(:,i+1)=x1jian(:,i)+k1pp(:,i+1)*eps1(i+1);%平滑
    x1jian(:,i+1)=x1jianp(:,i)+k1(:,i+1)*eps1(i+1);%滤波
    x1jianp(:,i+1)=(fai0-J1*H01)*x1jian(:,i+1)+(fai1-J1*H11)*x1jianpp(:,i+1)+(fai2-J1*H21)*x1jianpp2(:,i+1)+J1*y1(i+1);%预报
end 
t=1:Bushu;
figure
subplot(2,2,1);plot(t,x(1,t),'b',t,x1jian(1,t),'r:');
subplot(2,2,2);plot(t,x(2,t),'b',t,x1jian(2,t),'r:');
% -----------------传感器2----------------%
x2jian(:,1)=zeros(2,1);x2jianp(:,1)=zeros(2,1);x2jianpp(:,1)=zeros(2,1);
p21(:,:,1)=0.1*eye(2);p22(:,:,1)=0.1*eye(2);p25(:,:,1)=0.1*eye(2);p23(:,:,1)=0.1*eye(2);p26(:,:,1)=0.1*eye(2);p29(:,:,1)=0.1*eye(2);
for i=1:Bushu+5
    Qeps2(i+1)=H02*p21(:,:,i)*H02'+H02*p22(:,:,i)*H12'+H02*p23(:,:,i)*H22'+...
               H12*p22(:,:,i)'*H02'+H12*p25(:,:,i)*H12'+H12*p26(:,:,i)*H22'+...
               H22*p23(:,:,i)'*H02'+H22*p26(:,:,i)'*H12'+H22*p29(:,:,i)*H22'+Qv2;
    k2(:,i+1)=p21(:,:,i)*H02'*inv(Qeps2(i+1))+p22(:,:,i)*H12'*inv(Qeps2(i+1))+p23(:,:,i)*H22'*inv(Qeps2(i+1));%滤波增益         
    k2pp(:,i+1)=p22(:,:,i)'*H02'*inv(Qeps2(i+1))+p25(:,:,i)*H12'*inv(Qeps2(i+1))+p26(:,:,i)*H22'*inv(Qeps2(i+1));%平滑增益
    k2pp2(:,i+1)=p23(:,:,i)'*H02'*inv(Qeps2(i+1))+p26(:,:,i)'*H12'*inv(Qeps2(i+1))+p29(:,:,i)*H22'*inv(Qeps2(i+1));%2步平滑增益
    
    p26(:,:,i+1)=p22(:,:,i)-k2(:,i+1)*Qeps2(i+1)*k2pp(:,i+1)';
    p25(:,:,i+1)=p21(:,:,i)-k2(:,i+1)*Qeps2(i+1)*k2(:,i+1)';
    p29(:,:,i+1)=p25(:,:,i)-k2pp(:,i+1)*Qeps2(i+1)*k2pp(:,i+1)';
    p211(:,:,i+1)=p23(:,:,i)'- k2pp2(:,i+1)*Qeps2(i+1)*k2(:,i+1)';
    p212(:,:,i+1)=p26(:,:,i)'- k2pp2(:,i+1)*Qeps2(i+1)*k2pp(:,i+1)';
    p213(:,:,i+1)=p29(:,:,i)- k2pp2(:,i+1)*Qeps2(i+1)*k2pp2(:,i+1)';
    p210(:,:,i+1)=p211(:,:,i+1)*(fai0-J2*H02)+p212(:,:,i+1)*(fai1-J2*H12)+p213(:,:,i+1)*(fai2-J2*H22);
    p214(:,:,i+1)=p23(:,:,i)-k2(:,i+1)*Qeps2(i+1)*k2pp2(:,i+1)';
    p22(:,:,i+1)=(p25(:,:,i+1)*(fai0-J2*H02)'+p26(:,:,i+1)*(fai1-J2*H12)'+p214(:,:,i+1)*(fai2-J2*H22))';
    p215(:,:,i+1)=p26(:,:,i)-k2pp(:,i+1)*Qeps2(i+1)*k2pp2(:,i+1)';
    p23(:,:,i+1)=(p22(:,:,i+1)*(fai0-J2*H02)'+p29(:,:,i+1)*(fai1-J2*H12)'+p215(:,:,i+1)*(fai2-J2*H22))';
    p21(:,:,i+1)=(fai0-J2*H02)*p22(:,:,i+1)'+(fai1-J2*H12)*p23(:,:,i+1)'+(fai2-J2*H22)*p210(:,:,i+1)+Qw_2;
    
    eps2(i+1)=y2(i+1)-H02*x2jianp(:,i)-H12*x2jian(:,i)-H22*x2jianpp(:,i);
    x2jianpp2(:,i+1)=x2jianpp(:,i)+k2pp2(:,i+1)*eps2(i+1);%2步平滑
    x2jianpp(:,i+1)=x2jian(:,i)+k2pp(:,i+1)*eps2(i+1);%平滑
    x2jian(:,i+1)=x2jianp(:,i)+k2(:,i+1)*eps2(i+1);%滤波
    x2jianp(:,i+1)=(fai0-J2*H02)*x2jian(:,i+1)+(fai1-J2*H12)*x2jianpp(:,i+1)+(fai2-J2*H22)*x2jianpp2(:,i+1)+J2*y2(i+1);%预报
end 
t=1:Bushu;
figure
subplot(2,2,1);plot(t,x(1,t),'b',t,x2jian(1,t),'r:');
subplot(2,2,2);plot(t,x(2,t),'b',t,x2jian(2,t),'r:');
%-----------------传感器3----------------%
x3jian(:,1)=zeros(2,1);x3jianp(:,1)=zeros(2,1);x3jianpp(:,1)=zeros(2,1);
p31(:,:,1)=0.1*eye(2);p32(:,:,1)=0.1*eye(2);p35(:,:,1)=0.1*eye(2);p33(:,:,1)=0.1*eye(2);p36(:,:,1)=0.1*eye(2);p39(:,:,1)=0.1*eye(2);
for i=1:Bushu+5
    Qeps3(i+1)=H03*p31(:,:,i)*H03'+H03*p32(:,:,i)*H13'+H03*p33(:,:,i)*H23'+...
               H13*p32(:,:,i)'*H03'+H13*p35(:,:,i)*H13'+H13*p36(:,:,i)*H23'+...
               H23*p33(:,:,i)'*H03'+H23*p36(:,:,i)'*H13'+H23*p39(:,:,i)*H23'+Qv3;
    k3(:,i+1)=p31(:,:,i)*H03'*inv(Qeps3(i+1))+p32(:,:,i)*H13'*inv(Qeps3(i+1))+p33(:,:,i)*H23'*inv(Qeps3(i+1));%滤波增益         
    k3pp(:,i+1)=p32(:,:,i)'*H03'*inv(Qeps3(i+1))+p35(:,:,i)*H13'*inv(Qeps3(i+1))+p36(:,:,i)*H23'*inv(Qeps3(i+1));%平滑增益
    k3pp2(:,i+1)=p33(:,:,i)'*H03'*inv(Qeps3(i+1))+p36(:,:,i)'*H13'*inv(Qeps3(i+1))+p39(:,:,i)*H23'*inv(Qeps3(i+1));%2步平滑增益
    
    p36(:,:,i+1)=p32(:,:,i)-k3(:,i+1)*Qeps3(i+1)*k3pp(:,i+1)';
    p35(:,:,i+1)=p31(:,:,i)-k3(:,i+1)*Qeps3(i+1)*k3(:,i+1)';
    p39(:,:,i+1)=p35(:,:,i)-k3pp(:,i+1)*Qeps3(i+1)*k3pp(:,i+1)';
    p311(:,:,i+1)=p33(:,:,i)'- k3pp2(:,i+1)*Qeps3(i+1)*k3(:,i+1)';
    p312(:,:,i+1)=p36(:,:,i)'- k3pp2(:,i+1)*Qeps3(i+1)*k3pp(:,i+1)';
    p313(:,:,i+1)=p39(:,:,i)- k3pp2(:,i+1)*Qeps3(i+1)*k3pp2(:,i+1)';
    p310(:,:,i+1)=p311(:,:,i+1)*(fai0-J3*H03)+p312(:,:,i+1)*(fai1-J3*H13)+p313(:,:,i+1)*(fai2-J3*H23);
    p314(:,:,i+1)=p33(:,:,i)-k3(:,i+1)*Qeps3(i+1)*k3pp2(:,i+1)';
    p32(:,:,i+1)=(p35(:,:,i+1)*(fai0-J3*H03)'+p36(:,:,i+1)*(fai1-J3*H13)'+p314(:,:,i+1)*(fai2-J3*H23))';
    p315(:,:,i+1)=p36(:,:,i)-k3pp(:,i+1)*Qeps3(i+1)*k3pp2(:,i+1)';
    p33(:,:,i+1)=(p32(:,:,i+1)*(fai0-J3*H03)'+p39(:,:,i+1)*(fai1-J3*H13)'+p315(:,:,i+1)*(fai2-J3*H23))';
    p31(:,:,i+1)=(fai0-J3*H03)*p32(:,:,i+1)'+(fai1-J3*H13)*p33(:,:,i+1)'+(fai2-J3*H23)*p310(:,:,i+1)+Qw_3;
    
    eps3(i+1)=y3(i+1)-H03*x3jianp(:,i)-H13*x3jian(:,i)-H23*x3jianpp(:,i);
    x3jianpp2(:,i+1)=x3jianpp(:,i)+k3pp2(:,i+1)*eps3(i+1);%2步平滑
    x3jianpp(:,i+1)=x3jian(:,i)+k3pp(:,i+1)*eps3(i+1);%平滑
    x3jian(:,i+1)=x3jianp(:,i)+k3(:,i+1)*eps3(i+1);%滤波
    x3jianp(:,i+1)=(fai0-J3*H03)*x3jian(:,i+1)+(fai1-J3*H13)*x3jianpp(:,i+1)+(fai2-J3*H23)*x3jianpp2(:,i+1)+J3*y3(i+1);%预报
end 
t=1:Bushu;
figure
subplot(2,2,1);plot(t,x(1,t),'b',t,x3jian(1,t),'r:');
subplot(2,2,2);plot(t,x(2,t),'b',t,x3jian(2,t),'r:');
%-----------------迹----------------%
trace_p1=trace(p15(:,:,Bushu))
trace_p2=trace(p25(:,:,Bushu))
trace_p3=trace(p35(:,:,Bushu))
% %-----------------滤波误差互协方差阵----------------%
% P12(:,:,1)=eye(2);P13(:,:,1)=eye(2);P23(:,:,1)=eye(2);
% for i=1:Bushu  
%     PP12(:,:,i+1)=fai* P12(:,:,i)*fai'+gama*Qw*gama';%预报误差互协方差阵
%     PP13(:,:,i+1)=fai* P13(:,:,i)*fai'+gama*Qw*gama';
%     PP23(:,:,i+1)=fai* P23(:,:,i)*fai'+gama*Qw*gama';
%     %---滤波误差互协方差阵----%
%     P12(:,:,i+1)=[eye(2)-k1(:,i+1)*H01]*PP12(:,:,i+1)-[eye(2)-k1(:,i+1)*H01]*fai*P12(:,:,i)*H12'*k2(:,i+1)';P21(:,:,i+1)=P12(:,:,i+1)';
%     P13(:,:,i+1)=[eye(2)-k1(:,i+1)*H01]*PP13(:,:,i+1)*[eye(2)-k3(:,i+1)*H03]'-[eye(2)-k1(:,i+1)*H01]*fai*P13(:,:,i)*H13'*k3(:,i+1)';P31(:,:,i+1)=P13(:,:,i+1)';
%     P23(:,:,i+1)=PP23(:,:,i+1)*[eye(2)-k3(:,i+1)*H03]'+k1(:,i+1)*H12*P23(:,:,i)*H13'*k3(:,i+1)'-fai*P23(:,:,i)*H13'*k3(:,i+1)'-k2(:,i+1)*H12*P23(:,:,i)*fai'*[eye(2)-k3(:,i+1)*H03];P32(:,:,i+1)=P23(:,:,i+1)';
% end
% 
% %-----------------按矩阵加权---------------%
% for i=1:Bushu
%     Psigma(:,:,i)=[P1(:,:,i),P12(:,:,i),P13(:,:,i);
%                   P12(:,:,i)',P2(:,:,i),P23(:,:,i);
%                   P13(:,:,i)',P23(:,:,i)',P3(:,:,i)];
% end
% e=[eye(2),eye(2),eye(2)]';
% for i=1:Bushu
%     A(:,:,i)=inv(Psigma(:,:,i))*e*inv(e'*inv(Psigma(:,:,i))*e);
%     Pf(:,:,i)=inv(e'*inv(Psigma(:,:,i))*e);%误差方差阵
%     xfjian(:,i)=A(1:2,:,i)'*x1jian(:,i)+A(3:4,:,i)'*x2jian(:,i)+A(5:6,:,i)'*x3jian(:,i);
% end
% % t=1:Bushu;
% % figure
% % subplot(2,2,1);plot(t,x(1,t),'b',t,xfjian(1,t),'r:');
% % subplot(2,2,2);plot(t,x(2,t),'b',t,xfjian(2,t),'r:');
% %-----------------exf---------------%
% exf(300)=[0];
% exf(1)=(xfjian(1,1)-x(1,1))*(xfjian(1,1)-x(1,1));
% for i=2:Bushu
%    exf(i)=(xfjian(1,i)-x(1,i))*(xfjian(1,i)-x(1,i))+exf(i-1);
% end   
% % t=1:Bushu;
% % figure
% % subplot(2,2,1);plot(t,exf(1,t),'b');
% % subplot(2,2,2);plot(t,exf(2,t),'b');
% 
% %-----------------SCI----------------%
% deta=0.0001;pp1=P1(:,:,Bushu);pp2=P2(:,:,Bushu);pp3=P3(:,:,Bushu);
% [w13,Pci13]=y13_0618(deta,pp1,pp3)%pp1和pp3形成pci1
% [w132,Pci132,trPci132]=y132_0618(deta,Pci13,pp2)%pp2和pci1融合成psci132
% for i=1:Bushu
%      xci13(:,i)=Pci13*(w13*inv(pp1)*x1jian(:,i)+(1-w13)*inv(pp3)*x3jian(:,i));
%      xci132(:,i)=Pci132*(w132*inv(Pci13)*xci13(:,i)+(1-w132)*inv(pp2)*x2jian(:,i));%xci132
% end
% t=1:Bushu;
% figure
% plot(t,x(1,t),'b',t,xci132(1,t),'r:');
% figure
% plot(t,x(2,t),'b',t,xci132(2,t),'r:');axis([0,Bushu,-45,20]);
% %-----------------exci---------------%
% exci(300)=[0];
% exci(1)=(xci132(1,1)-x(1,1))*(xci132(1,1)-x(1,1));
% for i=2:Bushu
%    exci(i)=(xci132(1,i)-x(1,i))*(xci132(1,i)-x(1,i))+exci(i-1);
% end   
% t=1:Bushu;
% figure
% plot(t,exf(t),'b',t,exci(t),'r');axis([0,300,0,180]);
% 
% Pci132_=Pci132*(w132*w13*w132*w13*inv(pp1)*pp1*inv(pp1)+w132*w13*w132*(1-w13)*inv(pp1)*P13(:,:,Bushu)*inv(pp3)+...
%         w132*w13*(1-w132)*inv(pp1)*P12(:,:,Bushu)*inv(pp2)+w132*(1-w13)*w132*w13*inv(pp3)*P31(:,:,Bushu)*inv(pp1)+w132*(1-w13)*w132*(1-w13)*inv(pp3)+...
%         w132*(1-w13)*(1-w132)*inv(pp3)*P32(:,:,Bushu)*inv(pp2)+(1-w132)*w132*w13*inv(pp2)*P21(:,:,Bushu)*inv(pp1)+(1-w132)*w132*(1-w13)*inv(pp2)*P23(:,:,Bushu)*inv(pp3)+...
%         (1-w132)*(1-w132)*inv(pp2)*pp2*inv(pp2))*Pci132;
% %-----------------椭圆半径----------------%
% P1_ni=inv(P1(:,:,Bushu));P2_ni=inv(P2(:,:,Bushu));P3_ni=inv(P3(:,:,Bushu));
% Pci132_ni=inv(Pci132);Pci132__ni=inv(Pci132_);
% 
% theta=0:pi/100:2*pi;
% r1=1./sqrt(P1_ni(1,1)*cos(theta).^2+(P1_ni(1,2)+P1_ni(2,1))*cos(theta).*sin(theta)+P1_ni(2,2)*sin(theta).^2);
% r2=1./sqrt(P2_ni(1,1)*cos(theta).^2+(P2_ni(1,2)+P2_ni(2,1))*cos(theta).*sin(theta)+P2_ni(2,2)*sin(theta).^2);
% r3=1./sqrt(P3_ni(1,1)*cos(theta).^2+(P3_ni(1,2)+P3_ni(2,1))*cos(theta).*sin(theta)+P3_ni(2,2)*sin(theta).^2);
% 
% rci132=1./sqrt(Pci132_ni(1,1)*cos(theta).^2+(Pci132_ni(1,2)+Pci132_ni(2,1))*cos(theta).*sin(theta)+Pci132_ni(2,2)*sin(theta).^2);
% rci132_=1./sqrt(Pci132__ni(1,1)*cos(theta).^2+(Pci132__ni(1,2)+Pci132__ni(2,1))*cos(theta).*sin(theta)+Pci132__ni(2,2)*sin(theta).^2);
% %--------------作图----------------%
% t=1:Bushu;
% figure 
% hold on;
% polar(theta,r1,'b');
% polar(theta,r2,'k');
% polar(theta,r3,'r');
% polar(theta,rci132,'k');
% polar(theta,rci132_,'k-.');%axis([-2.3,2.3,-2.3,2.3]);
