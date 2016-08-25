%-----------------------------------------%
%名称：带有色噪声多传感器时滞系统CI融合估值器
%作者：王军
%时间：2016-07-28
%-----------------------------------------%
clc;clear all;
Bushu=200;
B1=0.6;B2=0.3;B3=0.3;
Qxi1=1;Qxi2=2;Qxi3=3;Qw=0.07;
randn('seed',2);w=sqrt(Qw)*randn(1,Bushu+10);
randn('seed',2);xi1=sqrt(Qxi1)*randn(1,Bushu+10);
randn('seed',5);xi2=sqrt(Qxi2)*randn(1,Bushu+10);
randn('seed',4);xi3=sqrt(Qxi3)*randn(1,Bushu+10);
%-------------------系统模型-------------------%
fai0=[0.65 0;-0.25 0.5];fai1=[-0.1 0.2;0.1 -0.1];fai2=[0.45 0.5;0 0.5];
gama=[1.4;-1.2];
h01=[0.2 1];h11=[0.1 2];h21=[0.1 0.8];h31=[0 0];
h02=[1 0.9];h12=[0.1 0];h22=[0 0.8];h32=[0 0];
h03=[-3 0.9];h13=[0.6 0];h23=[0.1 0];h33=[0 0];

H01=h01*fai0+h11-B1*h01;H11=h01*fai1+h21-B1*h11;H21=h01*fai2+h31-B1*h21;
H02=h02*fai0+h12-B2*h02;H12=h02*fai1+h22-B2*h12;H22=h02*fai2+h32-B2*h22;
H03=h03*fai0+h13-B3*h03;H13=h03*fai1+h23-B3*h13;H23=h03*fai2+h33-B3*h23;

v1=h01*gama*w+xi1;v2=h02*gama*w+xi2;v3=h03*gama*w+xi3;
%--------------------------------------------%
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
S1=Qw*gama'*h01';S2=Qw*gama'*h02';S3=Qw*gama'*h03';
Qv1=h01*gama*Qw*gama'*h01'+Qxi1;Qv2=h02*gama*Qw*gama'*h02'+Qxi2;Qv3=h03*gama*Qw*gama'*h03'+Qxi3;
%------------------局部Kalman滤波器---------------%
%----------------------传感器1--------------------%
x1jian(:,1)=zeros(2,1);x1jianp(:,1)=zeros(2,1);x1jianpp(:,1)=zeros(2,1);
p11(:,:,1)=0.1*eye(2);p12(:,:,1)=0.1*eye(2);p15(:,:,1)=0.1*eye(2);p13(:,:,1)=0.1*eye(2);p16(:,:,1)=0.1*eye(2);p19(:,:,1)=0.1*eye(2);
for i=1:Bushu+5
    Qeps1(i+1)=H01*p11(:,:,i)*H01'+H01*p12(:,:,i)*H11'+H01*p13(:,:,i)*H21'+...
               H11*p12(:,:,i)'*H01'+H11*p15(:,:,i)*H11'+H11*p16(:,:,i)*H21'+...
               H21*p13(:,:,i)'*H01'+H21*p16(:,:,i)'*H11'+H21*p19(:,:,i)*H21'+Qv1;
    k1(:,i+1)=p11(:,:,i)*H01'*inv(Qeps1(i+1))+p12(:,:,i)*H11'*inv(Qeps1(i+1))+p13(:,:,i)*H21'*inv(Qeps1(i+1));      %滤波增益         
    k1pp(:,i+1)=p12(:,:,i)'*H01'*inv(Qeps1(i+1))+p15(:,:,i)*H11'*inv(Qeps1(i+1))+p16(:,:,i)*H21'*inv(Qeps1(i+1));   %平滑增益
    k1pp2(:,i+1)=p13(:,:,i)'*H01'*inv(Qeps1(i+1))+p16(:,:,i)'*H11'*inv(Qeps1(i+1))+p19(:,:,i)*H21'*inv(Qeps1(i+1)); %2步平滑增益
    
    pw1(i+1)=Qw-S1*inv(Qeps1(i+1))*S1';
    
    p16(:,:,i+1)=p12(:,:,i)-k1(:,i+1)*Qeps1(i+1)*k1pp(:,i+1)';
    p15(:,:,i+1)=p11(:,:,i)-k1(:,i+1)*Qeps1(i+1)*k1(:,i+1)';
    p19(:,:,i+1)=p15(:,:,i)-k1pp(:,i+1)*Qeps1(i+1)*k1pp(:,i+1)';
    p111(:,:,i+1)=p13(:,:,i)'- k1pp2(:,i+1)*Qeps1(i+1)*k1(:,i+1)';
    p112(:,:,i+1)=p16(:,:,i)'- k1pp2(:,i+1)*Qeps1(i+1)*k1pp(:,i+1)';
    p113(:,:,i+1)=p19(:,:,i)- k1pp2(:,i+1)*Qeps1(i+1)*k1pp2(:,i+1)';
    p110(:,:,i+1)=p111(:,:,i+1)*fai0'+p112(:,:,i+1)*fai1'+p113(:,:,i+1)*fai2'-k1pp2(:,i+1)*S1'*gama';
    p114(:,:,i+1)=p13(:,:,i)-k1(:,i+1)*Qeps1(i+1)*k1pp2(:,i+1)';
    p12(:,:,i+1)=(p15(:,:,i+1)*fai0'+p16(:,:,i+1)*fai1'+p114(:,:,i+1)*fai2'-k1(:,i+1)*S1'*gama')';
    p115(:,:,i+1)=p16(:,:,i)-k1pp(:,i+1)*Qeps1(i+1)*k1pp2(:,i+1)';
    p13(:,:,i+1)=(p12(:,:,i+1)*fai0'+p19(:,:,i+1)*fai1'+p115(:,:,i+1)*fai2'-k1pp(:,i+1)*S1'*gama')';
    p11(:,:,i+1)=fai0*p12(:,:,i+1)'+fai1*p13(:,:,i+1)'+fai2*p110(:,:,i+1)-gama*S1*k1(:,i+1)'*fai0'-...
                 gama*S1*k1pp(:,i+1)'*fai1'-gama*S1*k1pp2(:,i+1)'*fai2'+gama*pw1(i+1)*gama';
    
    eps1(i+1)=y1(i+1)-H01*x1jianp(:,i)-H11*x1jian(:,i)-H21*x1jianpp(:,i);
    x1jianpp2(:,i+1)=x1jianpp(:,i)+k1pp2(:,i+1)*eps1(i+1);  %2步平滑
    x1jianpp(:,i+1)=x1jian(:,i)+k1pp(:,i+1)*eps1(i+1);      %平滑
    x1jian(:,i+1)=x1jianp(:,i)+k1(:,i+1)*eps1(i+1);         %滤波
    x1jianp(:,i+1)=fai0*x1jian(:,i+1)+fai1*x1jianpp(:,i+1)+fai2*x1jianpp2(:,i+1)+gama*S1*inv(Qeps1(i+1))*eps1(i+1); %预报
end 
t=1:Bushu;
figure
subplot(2,2,1);plot(t,x(1,t),'b',t,x1jianp(1,t),'r:');
subplot(2,2,2);plot(t,x(2,t),'b',t,x1jianp(2,t),'r:');
% --------------------传感器2---------------------%
x2jian(:,1)=zeros(2,1);x2jianp(:,1)=zeros(2,1);x2jianpp(:,1)=zeros(2,1);
p21(:,:,1)=0.1*eye(2);p22(:,:,1)=0.1*eye(2);p25(:,:,1)=0.1*eye(2);p23(:,:,1)=0.1*eye(2);p26(:,:,1)=0.1*eye(2);p29(:,:,1)=0.1*eye(2);
for i=1:Bushu+5
    Qeps2(i+1)=H02*p21(:,:,i)*H02'+H02*p22(:,:,i)*H12'+H02*p23(:,:,i)*H22'+...
               H12*p22(:,:,i)'*H02'+H12*p25(:,:,i)*H12'+H12*p26(:,:,i)*H22'+...
               H22*p23(:,:,i)'*H02'+H22*p26(:,:,i)'*H12'+H22*p29(:,:,i)*H22'+Qv2;
    k2(:,i+1)=p21(:,:,i)*H02'*inv(Qeps2(i+1))+p22(:,:,i)*H12'*inv(Qeps2(i+1))+p23(:,:,i)*H22'*inv(Qeps2(i+1));      %滤波增益         
    k2pp(:,i+1)=p22(:,:,i)'*H02'*inv(Qeps2(i+1))+p25(:,:,i)*H12'*inv(Qeps2(i+1))+p26(:,:,i)*H22'*inv(Qeps2(i+1));   %平滑增益
    k2pp2(:,i+1)=p23(:,:,i)'*H02'*inv(Qeps2(i+1))+p26(:,:,i)'*H12'*inv(Qeps2(i+1))+p29(:,:,i)*H22'*inv(Qeps2(i+1)); %2步平滑增益
    
    pw2(i+1)=Qw-S2*inv(Qeps2(i+1))*S2';
    
    p26(:,:,i+1)=p22(:,:,i)-k2(:,i+1)*Qeps2(i+1)*k2pp(:,i+1)';
    p25(:,:,i+1)=p21(:,:,i)-k2(:,i+1)*Qeps2(i+1)*k2(:,i+1)';
    p29(:,:,i+1)=p25(:,:,i)-k2pp(:,i+1)*Qeps2(i+1)*k2pp(:,i+1)';
    p211(:,:,i+1)=p23(:,:,i)'- k2pp2(:,i+1)*Qeps2(i+1)*k2(:,i+1)';
    p212(:,:,i+1)=p26(:,:,i)'- k2pp2(:,i+1)*Qeps2(i+1)*k2pp(:,i+1)';
    p213(:,:,i+1)=p29(:,:,i)- k2pp2(:,i+1)*Qeps2(i+1)*k2pp2(:,i+1)';
    p210(:,:,i+1)=p211(:,:,i+1)*fai0'+p212(:,:,i+1)*fai1'+p213(:,:,i+1)*fai2'-k2pp2(:,i+1)*S2'*gama';
    p214(:,:,i+1)=p23(:,:,i)-k2(:,i+1)*Qeps2(i+1)*k2pp2(:,i+1)';
    p22(:,:,i+1)=(p25(:,:,i+1)*fai0'+p26(:,:,i+1)*fai1'+p214(:,:,i+1)*fai2'-k2(:,i+1)*S2'*gama')';
    p215(:,:,i+1)=p26(:,:,i)-k2pp(:,i+1)*Qeps2(i+1)*k2pp2(:,i+1)';
    p23(:,:,i+1)=(p22(:,:,i+1)*fai0'+p29(:,:,i+1)*fai1'+p215(:,:,i+1)*fai2'-k2pp(:,i+1)*S2'*gama')';
    p21(:,:,i+1)=fai0*p22(:,:,i+1)'+fai1*p23(:,:,i+1)'+fai2*p210(:,:,i+1)-gama*S2*k2(:,i+1)'*fai0'-...
                 gama*S2*k2pp(:,i+1)'*fai1'-gama*S2*k2pp2(:,i+1)'*fai2'+gama*pw2(i+1)*gama';
    
    eps2(i+1)=y2(i+1)-H02*x2jianp(:,i)-H12*x2jian(:,i)-H22*x2jianpp(:,i);
    x2jianpp2(:,i+1)=x2jianpp(:,i)+k2pp2(:,i+1)*eps2(i+1);  %2步平滑
    x2jianpp(:,i+1)=x2jian(:,i)+k2pp(:,i+1)*eps2(i+1);      %平滑
    x2jian(:,i+1)=x2jianp(:,i)+k2(:,i+1)*eps2(i+1);         %滤波
    x2jianp(:,i+1)=fai0*x2jian(:,i+1)+fai1*x2jianpp(:,i+1)+fai2*x2jianpp2(:,i+1)+gama*S2*inv(Qeps2(i+1))*eps2(i+1); %预报
end 
t=1:Bushu;
figure
subplot(2,2,1);plot(t,x(1,t),'b',t,x2jianp(1,t),'r:');
subplot(2,2,2);plot(t,x(2,t),'b',t,x2jianp(2,t),'r:');
%------------------传感器3--------------------%
x3jian(:,1)=zeros(2,1);x3jianp(:,1)=zeros(2,1);x3jianpp(:,1)=zeros(2,1);
p31(:,:,1)=0.1*eye(2);p32(:,:,1)=0.1*eye(2);p35(:,:,1)=0.1*eye(2);p33(:,:,1)=0.1*eye(2);p36(:,:,1)=0.1*eye(2);p39(:,:,1)=0.1*eye(2);
for i=1:Bushu+5
    Qeps3(i+1)=H03*p31(:,:,i)*H03'+H03*p32(:,:,i)*H13'+H03*p33(:,:,i)*H23'+...
               H13*p32(:,:,i)'*H03'+H13*p35(:,:,i)*H13'+H13*p36(:,:,i)*H23'+...
               H23*p33(:,:,i)'*H03'+H23*p36(:,:,i)'*H13'+H23*p39(:,:,i)*H23'+Qv3;
    k3(:,i+1)=p31(:,:,i)*H03'*inv(Qeps3(i+1))+p32(:,:,i)*H13'*inv(Qeps3(i+1))+p33(:,:,i)*H23'*inv(Qeps3(i+1));      %滤波增益         
    k3pp(:,i+1)=p32(:,:,i)'*H03'*inv(Qeps3(i+1))+p35(:,:,i)*H13'*inv(Qeps3(i+1))+p36(:,:,i)*H23'*inv(Qeps3(i+1));   %平滑增益
    k3pp2(:,i+1)=p33(:,:,i)'*H03'*inv(Qeps3(i+1))+p36(:,:,i)'*H13'*inv(Qeps3(i+1))+p39(:,:,i)*H23'*inv(Qeps3(i+1)); %2步平滑增益
    
    pw3(i+1)=Qw-S3*inv(Qeps3(i+1))*S3';

    p36(:,:,i+1)=p32(:,:,i)-k3(:,i+1)*Qeps3(i+1)*k3pp(:,i+1)';
    p35(:,:,i+1)=p31(:,:,i)-k3(:,i+1)*Qeps3(i+1)*k3(:,i+1)';
    p39(:,:,i+1)=p35(:,:,i)-k3pp(:,i+1)*Qeps3(i+1)*k3pp(:,i+1)';
    p311(:,:,i+1)=p33(:,:,i)'- k3pp2(:,i+1)*Qeps3(i+1)*k3(:,i+1)';
    p312(:,:,i+1)=p36(:,:,i)'- k3pp2(:,i+1)*Qeps3(i+1)*k3pp(:,i+1)';
    p313(:,:,i+1)=p39(:,:,i)- k3pp2(:,i+1)*Qeps3(i+1)*k3pp2(:,i+1)';
    p310(:,:,i+1)=p311(:,:,i+1)*fai0'+p312(:,:,i+1)*fai1'+p313(:,:,i+1)*fai2'-k3pp2(:,i+1)*S3'*gama';
    p314(:,:,i+1)=p33(:,:,i)-k3(:,i+1)*Qeps3(i+1)*k3pp2(:,i+1)';
    p32(:,:,i+1)=(p35(:,:,i+1)*fai0'+p36(:,:,i+1)*fai1'+p314(:,:,i+1)*fai2'-k3(:,i+1)*S3'*gama')';
    p315(:,:,i+1)=p36(:,:,i)-k3pp(:,i+1)*Qeps3(i+1)*k3pp2(:,i+1)';
    p33(:,:,i+1)=(p32(:,:,i+1)*fai0'+p39(:,:,i+1)*fai1'+p315(:,:,i+1)*fai2'-k3pp(:,i+1)*S3'*gama')';
    p31(:,:,i+1)=fai0*p32(:,:,i+1)'+fai1*p33(:,:,i+1)'+fai2*p310(:,:,i+1)-gama*S3*k3(:,i+1)'*fai0'-...
                 gama*S3*k3pp(:,i+1)'*fai1'-gama*S3*k3pp2(:,i+1)'*fai2'+gama*pw3(i+1)*gama';
    
    eps3(i+1)=y3(i+1)-H03*x3jianp(:,i)-H13*x3jian(:,i)-H23*x3jianpp(:,i);
    x3jianpp2(:,i+1)=x3jianpp(:,i)+k3pp2(:,i+1)*eps3(i+1);  %2步平滑
    x3jianpp(:,i+1)=x3jian(:,i)+k3pp(:,i+1)*eps3(i+1);      %平滑
    x3jian(:,i+1)=x3jianp(:,i)+k3(:,i+1)*eps3(i+1);         %滤波
    x3jianp(:,i+1)=fai0*x3jian(:,i+1)+fai1*x3jianpp(:,i+1)+fai2*x3jianpp2(:,i+1)+gama*S3*inv(Qeps3(i+1))*eps3(i+1); %预报
end 
t=1:Bushu;
figure
subplot(2,2,1);plot(t,x(1,t),'b',t,x3jianp(1,t),'r:');
subplot(2,2,2);plot(t,x(2,t),'b',t,x3jianp(2,t),'r:');
%-------------------迹--------------------%
trace_p1=trace(p11(:,:,Bushu))
trace_p2=trace(p21(:,:,Bushu))
trace_p3=trace(p31(:,:,Bushu))
%-----------------滤波误差互协方差阵----------------%
Qv12=h01*gama*Qw*gama'*h02';Qv21=h02*gama*Qw*gama'*h01';
Qv13=h01*gama*Qw*gama'*h03';Qv31=h03*gama*Qw*gama'*h01';
Qv23=h02*gama*Qw*gama'*h03';Qv32=h03*gama*Qw*gama'*h02';

P121(:,:,1)=0.1*eye(2);P122(:,:,1)=0.1*eye(2);P123(:,:,1)=0.1*eye(2);P125(:,:,1)=0.1*eye(2);P126(:,:,1)=0.1*eye(2);P129(:,:,1)=0.1*eye(2);
P131(:,:,1)=0.1*eye(2);P132(:,:,1)=0.1*eye(2);P133(:,:,1)=0.1*eye(2);P135(:,:,1)=0.1*eye(2);P136(:,:,1)=0.1*eye(2);P139(:,:,1)=0.1*eye(2);
P231(:,:,1)=0.1*eye(2);P232(:,:,1)=0.1*eye(2);P233(:,:,1)=0.1*eye(2);P235(:,:,1)=0.1*eye(2);P236(:,:,1)=0.1*eye(2);P239(:,:,1)=0.1*eye(2);
for i=1:Bushu  
    Qeps12(i+1)=H01*P121(:,:,i)*H02'+H01*P122(:,:,i)*H12'+H01*P123(:,:,i)*H22'+...
               H11*P122(:,:,i)'*H02'+H11*P125(:,:,i)*H12'+H11*P126(:,:,i)*H22'+...
               H21*P123(:,:,i)'*H02'+H21*P126(:,:,i)'*H12'+H21*P129(:,:,i)*H22'+Qv12;
    P126(:,:,i+1)=P122(:,:,i)+k1(:,i+1)*Qeps12(i+1)*k2pp(:,i+1)'-k1(:,i+1)*(H01*P122(:,:,i)+H11*P125(:,:,i)+H21*P126(:,:,i)')-(P121(:,:,i)*H02'*k2pp(:,i+1)'+P122(:,:,i)*H12'*k2pp(:,i+1)'+P123(:,:,i)*H22'*k2pp(:,i+1)');
    P125(:,:,i+1)=P121(:,:,i)+k1(:,i+1)*Qeps12(i+1)*k2(:,i+1)'-k1(:,i+1)*(H01*P121(:,:,i)+H11*P122(:,:,i)'+H21*P123(:,:,i)')-(P121(:,:,i)*H02'*k2(:,i+1)'+P122(:,:,i)*H12'*k2(:,i+1)'+P123(:,:,i)*H22'*k2(:,i+1)');
    P129(:,:,i+1)=P125(:,:,i)+k1pp(:,i+1)*Qeps12(i+1)*k2pp(:,i+1)'-k1pp(:,i+1)*(H01*P122(:,:,i)+H11*P125(:,:,i)+H21*P126(:,:,i)')-(P122(:,:,i)'*H02'*k2pp(:,i+1)'+P125(:,:,i)*H12'*k2pp(:,i+1)'+P126(:,:,i)*H22'*k2pp(:,i+1)');
    P1211(:,:,i+1)=P123(:,:,i)'+ k1pp2(:,i+1)*Qeps12(i+1)*k2(:,i+1)'-k1pp2(:,i+1)*(H01*P121(:,:,i)+H11*P122(:,:,i)'+H21*P123(:,:,i)')-(P123(:,:,i)'*H02'*k2(:,i+1)'+P126(:,:,i)'*H12'*k2(:,i+1)'+P129(:,:,i)*H22'*k2(:,i+1)');
    P1212(:,:,i+1)=P126(:,:,i)'+ k1pp2(:,i+1)*Qeps12(i+1)*k2pp(:,i+1)'-k1pp2(:,i+1)*(H01*P122(:,:,i)+H11*P125(:,:,i)'+H21*P126(:,:,i)')-(P123(:,:,i)'*H02'*k2pp(:,i+1)'+P126(:,:,i)'*H12'*k2pp(:,i+1)'+P129(:,:,i)*H22'*k2pp(:,i+1)');
    P1213(:,:,i+1)=P129(:,:,i)+ k1pp2(:,i+1)*Qeps12(i+1)*k2pp2(:,i+1)'-k1pp2(:,i+1)*(H01*P123(:,:,i)+H11*P126(:,:,i)+H21*P129(:,:,i))-(P123(:,:,i)'*H02'*k2pp2(:,i+1)'+P126(:,:,i)'*H12'*k2pp2(:,i+1)'+P129(:,:,i)*H22'*k2pp2(:,i+1)');
    P1210(:,:,i+1)=P1211(:,:,i+1)*fai0'+P1212(:,:,i+1)*fai1'+P1213(:,:,i+1)*fai2'+(-k1pp2(:,i+1)*S1'-(P123(:,:,i)'*H02'+P126(:,:,i)'*H12'+P129(:,:,i)*H22')*(inv(Qeps2(i+1)))'*S2'+k1pp2(:,i+1)*Qeps12(i+1)*(inv(Qeps2(i+1)))'*S2')*gama';
    P1214(:,:,i+1)=P123(:,:,i)+k1(:,i+1)*Qeps12(i+1)*k2pp2(:,i+1)'-k1(:,i+1)*(H01*P123(:,:,i)+H11*P126(:,:,i)+H21*P129(:,:,i))-(P121(:,:,i)*H02'*k2pp2(:,i+1)'+P122(:,:,i)*H12'*k2pp2(:,i+1)'+P123(:,:,i)*H22'*k2pp2(:,i+1)');
    P122(:,:,i+1)=(P125(:,:,i+1)*fai0'+P126(:,:,i+1)*fai1'+P1214(:,:,i+1)*fai2'+(-k1(:,i+1)*S1'-(P121(:,:,i)*H02'+P122(:,:,i)*H12'+P123(:,:,i)*H22')*(inv(Qeps2(i+1)))'*S2'+k1(:,i+1)*Qeps12(i+1)*(inv(Qeps2(i+1)))'*S2')*gama')';
    P1215(:,:,i+1)=P126(:,:,i)+k1pp(:,i+1)*Qeps12(i+1)*k2pp2(:,i+1)'-k1pp(:,i+1)*(H01*P123(:,:,i)+H11*P126(:,:,i)+H21*P129(:,:,i))-(P122(:,:,i)'*H02'*k2pp2(:,i+1)'+P125(:,:,i)*H12'*k2pp2(:,i+1)'+P126(:,:,i)*H22'*k2pp2(:,i+1)');
    P123(:,:,i+1)=(P122(:,:,i+1)*fai0'+P129(:,:,i+1)*fai1'+P1215(:,:,i+1)*fai2'+(-k1pp(:,i+1)*S1'-(P122(:,:,i)'*H02'+P125(:,:,i)*H12'+P126(:,:,i)*H22')*(inv(Qeps2(i+1)))'*S2'+k1pp(:,i+1)*Qeps12(i+1)*(inv(Qeps2(i+1)))'*S2')*gama')';
    P121(:,:,i+1)=fai0*P122(:,:,i+1)'+fai1*P123(:,:,i+1)'+fai2*P1210(:,:,i+1)+...
        gama*(-S2*k2(:,i+1)'-S1*inv(Qeps1(i+1))*(H01*P121(:,:,i)+H11*P122(:,:,i)'+H21*P123(:,:,i)')+S1*inv(Qeps1(i+1))*Qeps12(i+1)*k2(:,i+1)')*fai0'+...;
        gama*(-S2*k2pp(:,i+1)'-S1*inv(Qeps1(i+1))*(H01*P122(:,:,i)+H11*P125(:,:,i)+H21*P126(:,:,i)')+S1*inv(Qeps1(i+1))*Qeps12(i+1)*k2pp(:,i+1)')*fai1'+...
        gama*(-S2*k2pp2(:,i+1)'-S1*inv(Qeps1(i+1))*(H01*P123(:,:,i)+H11*P126(:,:,i)+H21*P129(:,:,i))+S1*inv(Qeps1(i+1))*Qeps12(i+1)*k2pp2(:,i+1)')*fai2'+...
        gama*(Qw+S1*inv(Qeps1(i+1))*Qeps12(i+1)*(inv(Qeps2(i+1)))'*S2'-S2*inv(Qeps2(i+1))*S2'-S1*inv(Qeps1(i+1))*S1')*gama';
end

for i=1:Bushu  
    Qeps13(i+1)=H01*P131(:,:,i)*H03'+H01*P132(:,:,i)*H13'+H01*P133(:,:,i)*H23'+...
               H11*P132(:,:,i)'*H03'+H11*P135(:,:,i)*H13'+H11*P136(:,:,i)*H23'+...
               H21*P133(:,:,i)'*H03'+H21*P136(:,:,i)'*H13'+H21*P139(:,:,i)*H23'+Qv13;
    P136(:,:,i+1)=P132(:,:,i)+k1(:,i+1)*Qeps13(i+1)*k3pp(:,i+1)'-k1(:,i+1)*(H01*P132(:,:,i)+H11*P135(:,:,i)+H21*P136(:,:,i)')-(P131(:,:,i)*H03'*k3pp(:,i+1)'+P132(:,:,i)*H13'*k3pp(:,i+1)'+P133(:,:,i)*H23'*k3pp(:,i+1)');
    P135(:,:,i+1)=P131(:,:,i)+k1(:,i+1)*Qeps13(i+1)*k3(:,i+1)'-k1(:,i+1)*(H01*P131(:,:,i)+H11*P132(:,:,i)'+H21*P133(:,:,i)')-(P131(:,:,i)*H03'*k3(:,i+1)'+P132(:,:,i)*H13'*k3(:,i+1)'+P133(:,:,i)*H23'*k3(:,i+1)');
    P139(:,:,i+1)=P135(:,:,i)+k1pp(:,i+1)*Qeps13(i+1)*k3pp(:,i+1)'-k1pp(:,i+1)*(H01*P132(:,:,i)+H11*P135(:,:,i)+H21*P136(:,:,i)')-(P132(:,:,i)'*H03'*k3pp(:,i+1)'+P135(:,:,i)*H13'*k3pp(:,i+1)'+P136(:,:,i)*H23'*k3pp(:,i+1)');
    P1311(:,:,i+1)=P133(:,:,i)'+ k1pp2(:,i+1)*Qeps13(i+1)*k3(:,i+1)'-k1pp2(:,i+1)*(H01*P131(:,:,i)+H11*P132(:,:,i)'+H21*P133(:,:,i)')-(P133(:,:,i)'*H03'*k3(:,i+1)'+P136(:,:,i)'*H13'*k3(:,i+1)'+P139(:,:,i)*H23'*k3(:,i+1)');
    P1312(:,:,i+1)=P136(:,:,i)'+ k1pp2(:,i+1)*Qeps13(i+1)*k3pp(:,i+1)'-k1pp2(:,i+1)*(H01*P132(:,:,i)+H11*P135(:,:,i)'+H21*P136(:,:,i)')-(P133(:,:,i)'*H03'*k3pp(:,i+1)'+P136(:,:,i)'*H13'*k3pp(:,i+1)'+P139(:,:,i)*H23'*k3pp(:,i+1)');
    P1313(:,:,i+1)=P139(:,:,i)+ k1pp2(:,i+1)*Qeps13(i+1)*k3pp2(:,i+1)'-k1pp2(:,i+1)*(H01*P133(:,:,i)+H11*P136(:,:,i)+H21*P139(:,:,i))-(P133(:,:,i)'*H03'*k3pp2(:,i+1)'+P136(:,:,i)'*H13'*k3pp2(:,i+1)'+P139(:,:,i)*H23'*k3pp2(:,i+1)');
    P1310(:,:,i+1)=P1311(:,:,i+1)*fai0'+P1312(:,:,i+1)*fai1'+P1313(:,:,i+1)*fai2'+(-k1pp2(:,i+1)*S1'-(P133(:,:,i)'*H03'+P136(:,:,i')*H13'+P139(:,:,i)*H23')*(inv(Qeps3(i+1)))'*S3'+k1pp2(:,i+1)*Qeps13(i+1)*(inv(Qeps3(i+1)))'*S3')*gama';
    P1314(:,:,i+1)=P133(:,:,i)+k1(:,i+1)*Qeps13(i+1)*k3pp2(:,i+1)'-k1(:,i+1)*(H01*P133(:,:,i)+H11*P136(:,:,i)+H21*P139(:,:,i))-(P131(:,:,i)*H03'*k3pp2(:,i+1)'+P132(:,:,i)*H13'*k3pp2(:,i+1)'+P133(:,:,i)*H23'*k3pp2(:,i+1)');
    P132(:,:,i+1)=(P135(:,:,i+1)*fai0'+P136(:,:,i+1)*fai1'+P1314(:,:,i+1)*fai2'+(-k1(:,i+1)*S1'-(P131(:,:,i)*H03'+P132(:,:,i)*H13'+P133(:,:,i)*H23')*(inv(Qeps3(i+1)))'*S3'+k1(:,i+1)*Qeps13(i+1)*(inv(Qeps3(i+1)))'*S3')*gama')';
    P1315(:,:,i+1)=P136(:,:,i)+k1pp(:,i+1)*Qeps13(i+1)*k3pp2(:,i+1)'-k1pp(:,i+1)*(H01*P133(:,:,i)+H11*P136(:,:,i)+H21*P139(:,:,i))-(P132(:,:,i)'*H03'*k3pp2(:,i+1)'+P135(:,:,i)*H13'*k3pp2(:,i+1)'+P136(:,:,i)*H23'*k3pp2(:,i+1)');
    P133(:,:,i+1)=(P132(:,:,i+1)*fai0'+P139(:,:,i+1)*fai1'+P1315(:,:,i+1)*fai2'+(-k1pp(:,i+1)*S1'-(P132(:,:,i)'*H03'+P135(:,:,i)*H13'+P136(:,:,i)*H23')*(inv(Qeps3(i+1)))'*S3'+k1pp(:,i+1)*Qeps13(i+1)*(inv(Qeps3(i+1)))'*S3')*gama')';
    P131(:,:,i+1)=fai0*P132(:,:,i+1)'+fai1*P133(:,:,i+1)'+fai2*P1310(:,:,i+1)+...
        gama*(-S3*k3(:,i+1)'-S1*inv(Qeps1(i+1))*(H01*P131(:,:,i)+H11*P132(:,:,i)'+H21*P133(:,:,i)')+S1*inv(Qeps1(i+1))*Qeps13(i+1)*k3(:,i+1)')*fai0'+...;
        gama*(-S3*k3pp(:,i+1)'-S1*inv(Qeps1(i+1))*(H01*P132(:,:,i)+H11*P135(:,:,i)+H21*P136(:,:,i)')+S1*inv(Qeps1(i+1))*Qeps13(i+1)*k3pp(:,i+1)')*fai1'+...
        gama*(-S3*k3pp2(:,i+1)'-S1*inv(Qeps1(i+1))*(H01*P133(:,:,i)+H11*P136(:,:,i)+H21*P139(:,:,i))+S1*inv(Qeps1(i+1))*Qeps13(i+1)*k3pp2(:,i+1)')*fai2'+...
        gama*(Qw+S1*inv(Qeps1(i+1))*Qeps13(i+1)*(inv(Qeps3(i+1)))'*S3'-S3*inv(Qeps3(i+1))*S3'-S1*inv(Qeps1(i+1))*S1')*gama';
end

for i=1:Bushu  
    Qeps23(i+1)=H02*P231(:,:,i)*H03'+H02*P232(:,:,i)*H13'+H02*P233(:,:,i)*H23'+...
               H12*P232(:,:,i)'*H03'+H12*P235(:,:,i)*H13'+H12*P236(:,:,i)*H23'+...
               H22*P233(:,:,i)'*H03'+H22*P236(:,:,i)'*H13'+H22*P239(:,:,i)*H23'+Qv23;
    P236(:,:,i+1)=P232(:,:,i)+k2(:,i+1)*Qeps23(i+1)*k3pp(:,i+1)'-k2(:,i+1)*(H02*P232(:,:,i)+H12*P235(:,:,i)+H22*P236(:,:,i)')-(P231(:,:,i)*H03'*k3pp(:,i+1)'+P232(:,:,i)*H13'*k3pp(:,i+1)'+P233(:,:,i)*H23'*k3pp(:,i+1)');
    P235(:,:,i+1)=P231(:,:,i)+k2(:,i+1)*Qeps23(i+1)*k3(:,i+1)'-k2(:,i+1)*(H02*P231(:,:,i)+H12*P232(:,:,i)'+H22*P233(:,:,i)')-(P231(:,:,i)*H03'*k3(:,i+1)'+P232(:,:,i)*H13'*k3(:,i+1)'+P233(:,:,i)*H23'*k3(:,i+1)');
    P239(:,:,i+1)=P235(:,:,i)+k2pp(:,i+1)*Qeps23(i+1)*k3pp(:,i+1)'-k2pp(:,i+1)*(H02*P232(:,:,i)+H12*P235(:,:,i)+H22*P236(:,:,i)')-(P232(:,:,i)'*H03'*k3pp(:,i+1)'+P235(:,:,i)*H13'*k3pp(:,i+1)'+P236(:,:,i)*H23'*k3pp(:,i+1)');
    P2311(:,:,i+1)=P233(:,:,i)'+ k2pp2(:,i+1)*Qeps23(i+1)*k3(:,i+1)'-k2pp2(:,i+1)*(H02*P231(:,:,i)+H12*P232(:,:,i)'+H22*P233(:,:,i)')-(P233(:,:,i)'*H03'*k3(:,i+1)'+P236(:,:,i)'*H13'*k3(:,i+1)'+P239(:,:,i)*H23'*k3(:,i+1)');
    P2312(:,:,i+1)=P236(:,:,i)'+ k2pp2(:,i+1)*Qeps23(i+1)*k3pp(:,i+1)'-k2pp2(:,i+1)*(H02*P232(:,:,i)+H12*P235(:,:,i)'+H22*P236(:,:,i)')-(P233(:,:,i)'*H03'*k3pp(:,i+1)'+P236(:,:,i)'*H13'*k3pp(:,i+1)'+P239(:,:,i)*H23'*k3pp(:,i+1)');
    P2313(:,:,i+1)=P239(:,:,i)+ k2pp2(:,i+1)*Qeps23(i+1)*k3pp2(:,i+1)'-k2pp2(:,i+1)*(H02*P233(:,:,i)+H12*P236(:,:,i)+H22*P239(:,:,i))-(P233(:,:,i)'*H03'*k3pp2(:,i+1)'+P236(:,:,i)'*H13'*k3pp2(:,i+1)'+P239(:,:,i)*H23'*k3pp2(:,i+1)');
    P2310(:,:,i+1)=P2311(:,:,i+1)*fai0'+P2312(:,:,i+1)*fai1'+P2313(:,:,i+1)*fai2'+(-k2pp2(:,i+1)*S2'-(P233(:,:,i)'*H03'+P236(:,:,i')*H13'+P239(:,:,i)*H23')*(inv(Qeps3(i+1)))'*S3'+k2pp2(:,i+1)*Qeps23(i+1)*(inv(Qeps3(i+1)))'*S3')*gama';
    P2314(:,:,i+1)=P233(:,:,i)+k2(:,i+1)*Qeps23(i+1)*k3pp2(:,i+1)'-k2(:,i+1)*(H02*P233(:,:,i)+H12*P236(:,:,i)+H22*P239(:,:,i))-(P231(:,:,i)*H03'*k3pp2(:,i+1)'+P232(:,:,i)*H13'*k3pp2(:,i+1)'+P233(:,:,i)*H23'*k3pp2(:,i+1)');
    P232(:,:,i+1)=(P235(:,:,i+1)*fai0'+P236(:,:,i+1)*fai1'+P2314(:,:,i+1)*fai2'+(-k2(:,i+1)*S2'-(P231(:,:,i)*H03'+P232(:,:,i)*H13'+P233(:,:,i)*H23')*(inv(Qeps3(i+1)))'*S3'+k2(:,i+1)*Qeps23(i+1)*(inv(Qeps3(i+1)))'*S3')*gama')';
    P2315(:,:,i+1)=P236(:,:,i)+k2pp(:,i+1)*Qeps23(i+1)*k3pp2(:,i+1)'-k2pp(:,i+1)*(H02*P233(:,:,i)+H12*P236(:,:,i)+H22*P239(:,:,i))-(P232(:,:,i)'*H03'*k3pp2(:,i+1)'+P235(:,:,i)*H13'*k3pp2(:,i+1)'+P236(:,:,i)*H23'*k3pp2(:,i+1)');
    P233(:,:,i+1)=(P232(:,:,i+1)*fai0'+P239(:,:,i+1)*fai1'+P2315(:,:,i+1)*fai2+(-k2pp(:,i+1)*S2'-(P232(:,:,i)'*H03'+P235(:,:,i)*H13'+P236(:,:,i)*H23')*(inv(Qeps3(i+1)))'*S3'+k2pp(:,i+1)*Qeps23(i+1)*(inv(Qeps3(i+1)))'*S3')*gama')';
    P231(:,:,i+1)=fai0*P232(:,:,i+1)'+fai1*P233(:,:,i+1)'+fai2*P2310(:,:,i+1)+...
        gama*(-S3*k3(:,i+1)'-S2*inv(Qeps2(i+1))*(H02*P231(:,:,i)+H12*P232(:,:,i)'+H22*P233(:,:,i)')+S2*inv(Qeps2(i+1))*Qeps23(i+1)*k3(:,i+1)')*fai0'+...;
        gama*(-S3*k3pp(:,i+1)'-S2*inv(Qeps2(i+1))*(H02*P232(:,:,i)+H12*P235(:,:,i)+H22*P236(:,:,i)')+S2*inv(Qeps2(i+1))*Qeps23(i+1)*k3pp(:,i+1)')*fai1'+...
        gama*(-S3*k3pp2(:,i+1)'-S2*inv(Qeps2(i+1))*(H02*P233(:,:,i)+H12*P236(:,:,i)+H22*P239(:,:,i))+S2*inv(Qeps2(i+1))*Qeps23(i+1)*k3pp2(:,i+1)')*fai2'+...
        gama*(Qw+S2*inv(Qeps2(i+1))*Qeps23(i+1)*(inv(Qeps3(i+1)))'*S3'-S3*inv(Qeps3(i+1))*S3'-S2*inv(Qeps2(i+1))*S2')*gama';
end
%-----------------按矩阵加权---------------%
for i=2:Bushu+1
    Psigma(:,:,i)=[p11(:,:,i-1),P121(:,:,i-1),P131(:,:,i-1);
                  P121(:,:,i-1)',p21(:,:,i-1),P231(:,:,i-1);
                  P131(:,:,i-1)',P231(:,:,i-1)',p31(:,:,i-1)];
end
e=[eye(2),eye(2),eye(2)]';
for i=3:Bushu+1
    A(:,:,i)=inv(Psigma(:,:,i))*e*inv(e'*inv(Psigma(:,:,i))*e);
    Pm(:,:,i-1)=inv(e'*inv(Psigma(:,:,i))*e);%误差方差阵
    xfjianp(:,i-1)=A(1:2,:,i)'*x1jianp(:,i-1)+A(3:4,:,i)'*x2jianp(:,i-1)+A(5:6,:,i)'*x3jianp(:,i-1);
end
trace_Pm=trace(Pm(:,:,Bushu))
t=1:Bushu;
figure
subplot(2,2,1);plot(t,x(1,t),'b',t,xfjianp(1,t),'r:');
subplot(2,2,2);plot(t,x(2,t),'b',t,xfjianp(2,t),'r:');
%-----------------SCI----------------%
deta=0.0001;pp1=p11(:,:,Bushu);pp2=p21(:,:,Bushu);pp3=p31(:,:,Bushu);
[w13,Pci13]=y13_0618(deta,pp1,pp3)                  %pp1和pp3形成pci1
[w132,Pci132,trPci132]=y132_0618(deta,Pci13,pp2)    %pp2和pci1融合成psci132
for i=1:Bushu
     xci13(:,i)=Pci13*(w13*inv(pp1)*x1jianp(:,i)+(1-w13)*inv(pp3)*x3jianp(:,i));
     xci132(:,i)=Pci132*(w132*inv(Pci13)*xci13(:,i)+(1-w132)*inv(pp2)*x2jianp(:,i));%xci132
end
t=1:Bushu;
figure
subplot(2,2,1);plot(t,x(1,t),'b',t,xci132(1,t),'r:');
subplot(2,2,2);plot(t,x(2,t),'b',t,xci132(2,t),'r:');
%-----------------ex---------------%
% ex1(:,:,1)=(x1jianp(:,1)-x(:,1))*(x1jianp(:,1)-x(:,1))';
% ex2(:,:,1)=(x2jianp(:,1)-x(:,1))*(x2jianp(:,1)-x(:,1))';
% ex3(:,:,1)=(x3jianp(:,1)-x(:,1))*(x3jianp(:,1)-x(:,1))';
% exf(:,:,1)=(xfjianp(:,1)-x(:,1))*(xfjianp(:,1)-x(:,1))';
% exc(:,:,1)=(xci132(:,1)-x(:,1))*(xci132(:,1)-x(:,1))';
% for i=2:Bushu
%    ex1(:,:,i)=(x1jianp(:,i)-x(:,i))*(x1jianp(:,i)-x(:,i))'+ex1(:,:,i-1); 
%    ex2(:,:,i)=(x2jianp(:,i)-x(:,i))*(x2jianp(:,i)-x(:,i))'+ex2(:,:,i-1);  
%    ex3(:,:,i)=(x3jianp(:,i)-x(:,i))*(x3jianp(:,i)-x(:,i))'+ex3(:,:,i-1);
%    exf(:,:,i)=(xfjianp(:,i)-x(:,i))*(xfjianp(:,i)-x(:,i))'+exf(:,:,i-1);
%    exc(:,:,i)=(xci132(:,i)-x(:,i))*(xci132(:,i)-x(:,i))'+exc(:,:,i-1);
% end   
% for i=2:Bushu
%     e11(i)=ex1(1,1,i);e12(i)=ex2(1,1,i);e13(i)=ex1(1,1,i);e1f(i)=exf(1,1,i);e1c(i)=exc(1,1,i);
%     e21(i)=ex1(2,2,i);e22(i)=ex2(2,2,i);e23(i)=ex1(2,2,i);e2f(i)=exf(2,2,i);e2c(i)=exc(2,2,i);
% end
% t=2:Bushu;
% figure
% subplot(2,2,1);plot(t,e11(t),'m',t,e12(t),'b',t,e13(t),'r',t,e1f(t),'k',t,e1c(t),'k:');
% subplot(2,2,2);plot(t,e21(t),'m',t,e22(t),'b',t,e23(t),'r',t,e2f(t),'k',t,e2c(t),'k:');

ex11(1)=(x1jianp(1,1)-x(1,1))*(x1jianp(1,1)-x(1,1));
for i=2:Bushu
   ex11(i)=(x1jianp(1,i)-x(1,i))*(x1jianp(1,i)-x(1,i))+ex11(i-1);
end   
ex12(1)=(x2jianp(1,1)-x(1,1))*(x2jianp(1,1)-x(1,1));
for i=2:Bushu
   ex12(i)=(x2jianp(1,i)-x(1,i))*(x2jianp(1,i)-x(1,i))+ex12(i-1);
end   
ex13(1)=(x3jianp(1,1)-x(1,1))*(x3jianp(1,1)-x(1,1));
for i=2:Bushu
   ex13(i)=(x3jianp(1,i)-x(1,i))*(x3jianp(1,i)-x(1,i))+ex13(i-1);
end   
ex1f(1)=(xfjianp(1,1)-x(1,1))*(xfjianp(1,1)-x(1,1));
for i=2:Bushu
   ex1f(i)=(xfjianp(1,i)-x(1,i))*(xfjianp(1,i)-x(1,i))+ex1f(i-1);
end   
ex1c(1)=(xci132(1,1)-x(1,1))*(xci132(1,1)-x(1,1));
for i=2:Bushu
   ex1c(i)=(xci132(1,i)-x(1,i))*(xci132(1,i)-x(1,i))+ex1c(i-1);
end   
t=1:Bushu;
figure
plot(t,ex11(t),'m',t,ex12(t),'b',t,ex13(t),'r',t,ex1f(t),'k',t,ex1c(t),'k:');


ex21(1)=(x1jianp(2,1)-x(2,1))*(x1jianp(2,1)-x(2,1));
for i=2:Bushu
   ex21(i)=(x1jianp(2,i)-x(2,i))*(x1jianp(2,i)-x(2,i))+ex21(i-1);
end   
ex22(1)=(x2jianp(2,1)-x(2,1))*(x2jianp(2,1)-x(2,1));
for i=2:Bushu
   ex22(i)=(x2jianp(2,i)-x(2,i))*(x2jianp(2,i)-x(2,i))+ex22(i-1);
end   
ex23(1)=(x3jianp(2,1)-x(2,1))*(x3jianp(2,1)-x(2,1));
for i=2:Bushu
   ex23(i)=(x3jianp(2,i)-x(2,i))*(x3jianp(2,i)-x(2,i))+ex23(i-1);
end   
ex2f(1)=(xfjianp(2,1)-x(2,1))*(xfjianp(2,1)-x(2,1));
for i=2:Bushu
   ex2f(i)=(xfjianp(2,i)-x(2,i))*(xfjianp(2,i)-x(2,i))+ex2f(i-1);
end   
ex2c(1)=(xci132(2,1)-x(2,1))*(xci132(2,1)-x(2,1));
for i=2:Bushu
   ex2c(i)=(xci132(2,i)-x(2,i))*(xci132(2,i)-x(2,i))+ex2c(i-1);
end   
t=1:Bushu;
figure
plot(t,ex21(t),'m',t,ex22(t),'b',t,ex23(t),'r',t,ex2f(t),'k',t,ex2c(t),'k:');
%--------------实际误差方差阵-------------% 
Pci132_=Pci132*(w132*w13*w132*w13*inv(pp1)*pp1*inv(pp1)+w132*w13*w132*(1-w13)*inv(pp1)*P131(:,:,Bushu)*inv(pp3)+...
        w132*w13*(1-w132)*inv(pp1)*P121(:,:,Bushu)*inv(pp2)+w132*(1-w13)*w132*w13*inv(pp3)*P131(:,:,Bushu)'*inv(pp1)+w132*(1-w13)*w132*(1-w13)*inv(pp3)+...
        w132*(1-w13)*(1-w132)*inv(pp3)*P231(:,:,Bushu)'*inv(pp2)+(1-w132)*w132*w13*inv(pp2)*P121(:,:,Bushu)'*inv(pp1)+(1-w132)*w132*(1-w13)*inv(pp2)*P231(:,:,Bushu)*inv(pp3)+...
        (1-w132)*(1-w132)*inv(pp2)*pp2*inv(pp2))*Pci132;
%-----------------椭圆半径-----------------%
P1_ni=inv(pp1);P2_ni=inv(pp2);P3_ni=inv(pp3);
Pci132_ni=inv(Pci132);
Pci132__ni=inv(Pci132_);
Pm_ni=inv(Pm(:,:,Bushu));
theta=0:pi/100:2*pi;
r1=1./sqrt(P1_ni(1,1)*cos(theta).^2+(P1_ni(1,2)+P1_ni(2,1))*cos(theta).*sin(theta)+P1_ni(2,2)*sin(theta).^2);
r2=1./sqrt(P2_ni(1,1)*cos(theta).^2+(P2_ni(1,2)+P2_ni(2,1))*cos(theta).*sin(theta)+P2_ni(2,2)*sin(theta).^2);
r3=1./sqrt(P3_ni(1,1)*cos(theta).^2+(P3_ni(1,2)+P3_ni(2,1))*cos(theta).*sin(theta)+P3_ni(2,2)*sin(theta).^2);
rm=1./sqrt(Pm_ni(1,1)*cos(theta).^2+(Pm_ni(1,2)+Pm_ni(2,1))*cos(theta).*sin(theta)+Pm_ni(2,2)*sin(theta).^2);
rci132=1./sqrt(Pci132_ni(1,1)*cos(theta).^2+(Pci132_ni(1,2)+Pci132_ni(2,1))*cos(theta).*sin(theta)+Pci132_ni(2,2)*sin(theta).^2);
rci132_=1./sqrt(Pci132__ni(1,1)*cos(theta).^2+(Pci132__ni(1,2)+Pci132__ni(2,1))*cos(theta).*sin(theta)+Pci132__ni(2,2)*sin(theta).^2);
%--------------作图----------------%
t=1:Bushu;
figure 
hold on;
polar(theta,r1,'m');
polar(theta,r2,'b');
polar(theta,r3,'r');
polar(theta,rm,'k');
polar(theta,rci132,'k:');
polar(theta,rci132_,'k-.');