%%model-1 拟静力学部分
clear 
close all
clc
global alfa A Zn F K_oj K_ij f D omega_i m I dm Q_1j Q_2j alfa_i alfa_o Ic rate_ad mu_o mu_i
%参数输入
f=[0.523 0.523];%沟道曲率
dm=53.9;        %节圆直径 单位：mm
D=7.3;         %钢球直径 单位：mm
alfa=15;        %初始接触角 /°度
Zn=18;
den=7800;%钢球密度 kg/mm^3
m=pi*den*(D*1e-3)^3/6;%钢球质量 kg
I=m*(D*1e-3)^2/10;%钢球的转动惯量 kg*m^2
Ic=I+0.25*m*(dm*1e-3)^2;%钢球绕轴承轴线的转动惯量
K_ij =7.5806e+05;%单位 N*mm^1.5
K_oj =7.7985e+05;%单位 N*mm^1.5
A=(f(1)+f(2)-1)*D;%曲率中心距离 mm
%系统输入
dF=1;W_x=[];W_y=[];W_z=[];W_c=[];rate=[];alfa_oi=[];rate_ad=1;
mu=[];
% for k=1:3
 for F=500;
omega_i=10000*pi/30;%内圈转速 rad/s
T=2*pi/omega_i;%内圈转动周期 s
deta_a=0.0101;
X=[0.0484 0.1622 0.0014 0.0014 deta_a];%初值向量：5*1
options =optimoptions('fsolve','Algorithm','Levenberg-Marquardt');
[Y,fval,exitflag]=fsolve(@fun,X,options);
if exitflag==1
    disp('结果收敛');
else
    disp('所设初值无效');
end
Q_1j=K_oj*Y(3)^1.5;%单位/N
Q_2j=K_ij*Y(4)^1.5;%单位/N
A_2j=A*cosd(alfa);
cos_1j=Y(2)/((f(1)-0.5)*D+Y(3));
cos_2j=(A_2j-Y(2))/((f(2)-0.5)*D+Y(4));
alfa_i=acos(cos_2j);
alfa_o=acos(cos_1j);
alfa_oi(1,dF)=alfa_i;
alfa_oi(2,dF)=alfa_o;
W=2*[450 1600 -250 227.9*0.9 2*pi/Zn 82*0.9 0];%角速度初值向量：7*1
tspan = 0:1e-5:0.3;
% options = odeset('RelTol',1e-3,'AbsTol',1e-6);%'Maxstep',1e-20,
[t,y] = ode15s(@fun1,tspan,W,options);
% num=size(y,1);
% rate_ad=y(num,4)/227
% end
%   num=size(y,1);
%   uox=0.5*y(num,4)*1e-3*(dm)+0.5*D*1e-3*(y(num,3)*cos(alfa_o)-y(num,2)*sin(alfa_o)+y(num,4)*cos(alfa_o))
%   uix=0.5*1e-3*(-y(num,4)+omega_i)*(dm)+0.5*D*1e-3*(y(num,3)*cos(alfa_i)-y(num,2)*sin(alfa_i)-(-y(num,4)+omega_i)*cos(alfa_i))
%   wso=y(num,3)*sin(alfa_o)+y(num,2)*cos(alfa_o)%-y(num,4)*sin(alfa_o)
%   wsi=y(num,3)*sin(alfa_i)+y(num,2)*cos(alfa_i)%+(y(num,4)-omega_i)*sin(alfa_i)
%   W_R=sqrt((y(num,1))^2+(y(num,2))^2+(y(num,3))^2)
%   beta_xz=atan(y(num,2)/y(num,3))*180/pi
%   beta_yz=atan(y(num,1)/y(num,3))*180/pi
mu(1,dF)=mu_o;
mu(2,dF)=mu_i;
num=length(t);
rate(1,dF)=y(num,4)/omega_i;
W_x(1:num,dF)=y(:,1);
W_y(1:num,dF)=y(:,2);
W_z(1:num,dF)=y(:,3);
W_c(1:num,dF)=y(:,4);
dF=dF+1;
end

%%附码
% alfa_o=alfa_oi(1,k);
% alfa_i=alfa_oi(2,k);
%保持架速度
% F_y=20:20:600;
% F_the=0:20:600;
% num2=size(F_the);
% num1=size(W_c,1);
% theory=227.5*ones(num2);
% plot(F_the,theory,'-r')
% hold on
% grid on
% plot(F_y,W_c(num1,:),'*')
% %滑动速度提取
num1=size(W_x,1);
num2=size(W_x,2);
uox=[];uix=[];
for k=1:num2
    alfa_o=alfa_oi(1,k);
alfa_i=alfa_oi(2,k);
  W=[W_x(num1,k) W_y(num1,k) W_z(num1,k) W_c(num1,k)];
  uox(k)=0.5*W(4)*1e-3*(dm)+0.5*D*1e-3*(W(3)*cos(alfa_o)-W(2)*sin(alfa_o)+W(4)*cos(alfa_o));
  uix(k)=0.5*1e-3*(-W(4)+omega_i)*(dm)+0.5*D*1e-3*(W(3)*cos(alfa_i)-W(2)*sin(alfa_i)-(-W(4)+omega_i)*cos(alfa_i));
end
F_y=25:25:600;
plot(F_y,uox,'-*');
hold on
plot(F_y,uix,'-*');
grid on
% % %自旋速度提取
% num1=size(W_x,1);
% num2=size(W_x,2);
% wso=[];wsi=[];
% for k=1:num2
%   alfa_o=alfa_oi(1,k);
%   alfa_i=alfa_oi(2,k);
%   W=[W_x(num1,k) W_y(num1,k) W_z(num1,k) W_c(num1,k)];
%   wso(k)=W(3)*sin(alfa_o)+W(2)*cos(alfa_o)+W(4)*sin(alfa_o);
%   wsi(k)=W(3)*sin(alfa_i)+W(2)*cos(alfa_i)+(W(4)-omega_i)*sin(alfa_i);
% end
% F_y=20:20:600;
% plot(F_y,-wso,'b-*');
% hold on
% plot(F_y,-wsi,'r-*')
% grid on
% % % % %姿态角变化
% num1=size(W_x,1);
% num2=size(W_x,2);
% beta_xz=[];beta_yz=[];  W_R=[];
% for k=1:num2
%   W=[W_x(num1,k) W_y(num1,k) W_z(num1,k) W_c(num1,k)];
%   W_R(k)=sqrt(W(1)^2+W(2)^2+W(3)^2);
%   beta_xz(k)=atan(W(2)/W(3))*180/pi;
%   beta_yz(k)=atan(W(1)/W(3))*180/pi;
% end
% F_y=20:20:600;
% plot(F_y,-beta_xz,'b-*');
% hold on
% plot(F_y,-beta_yz,'r-*');
% grid on
% figure(2)
% plot(F_y,W_R,'r-*');
%%滑滚比的计算
num1=size(W_x,1);
num2=size(W_x,2);
uox=[];uix=[];u1=[];u2=[];S1=[];S2=[];
for k=1:num2
    alfa_o=alfa_oi(1,k);
alfa_i=alfa_oi(2,k);
  W=[W_x(num1,k) W_y(num1,k) W_z(num1,k) W_c(num1,k)];
  uox(k)=0.5*W(4)*1e-3*(dm)+0.5*D*1e-3*(W(3)*cos(alfa_o)...
      -W(2)*sin(alfa_o)+W(4)*cos(alfa_o));
  uix(k)=0.5*1e-3*(-W(4)+omega_i)*(dm)+0.5*D*1e-3*(W(3)*cos(alfa_i)...
      -W(2)*sin(alfa_i)-(-W(4)+omega_i)*cos(alfa_i));
  u1(k)=0.5*W(4)*1e-3*(dm+D*cos(alfa_o))+0.5*D*1e-3*(W(2)*sin(alfa_o)...
      -W(3)*cos(alfa_o));
  u2(k)=0.5*1e-3*(-W(4)+omega_i)*(dm-D*cos(alfa_i))+0.5*D...
      *1e-3*(W(2)*sin(alfa_i)-W(3)*cos(alfa_i));
  S1(k)=2*uox(k)/u1(k);
  S2(k)=2*uix(k)/u2(k);
end
plot(-S1,-mu(1,:),'-o');
hold on
plot(S2,mu(2,:),'-s');