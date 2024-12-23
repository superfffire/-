%%
%% 刚度需要自行导入，用势能法计算的刚度即可，转角也需要自行导入，导入转角时注释掉代码中的y_c即可。
%%
clc
clear all
%% 齿轮参数
zp=80;
zg=80;
%% 质量
mp = 1;
mg = 1;
m=3;                                  %模数
Ip=0.002230;                             %主动轮转动惯量
Ig=0.002230;                            %从动轮转动惯量
rp=0.5*m*zp*cos(20*pi/180)/1000;       %主动轮基圆半径
rg=0.5*m*zg*cos(20*pi/180)/1000;       %从动轮基圆半径
%% 等效质量
mh=Ip*Ig/(Ig*rp^2+Ip*rg^2);
%% 正常刚度导入
% kw2=xlsread('K_betai_K_M.xlsx',1, 'B2:B302');%此处导入自己计算的时变啮合刚度刚度
kw2 = readmatrix('K_betai_K_M.csv', 'Range', 'B2:B1002');
kx2=y_z();%压力角转转角
%% 轴承参数
nball=8;
Dbearing=0.03765;          %%轴承滚道直径
dbearing=0.0087;           %%滚子直径
%% 外力
Tg=40;
Tp=40;
%% 转速
w1=8*2*pi;%主动轮
w2=w1*zp/zg;%从动轮
%% 调用ode45
T=2*pi/zp;
options = odeset('RelTol',1e-4);
x0=[0;0 ;0 ;0 ;0 ;w1; 0 ;w2];
[t,x] = ode45('motion_model',[0 36*T],x0,flag,zp,kw2,mp,mg,Ip,Ig,rg,rp,mh,Tg,Tp,kx2); 
%% 刚度导入判断
n=6000;
k=[1 n]';
%% 动态啮合力
F=[1 n]';
for i=1:n
    [k(i),F(i)]=shuntaigangdu(kx2,zp,kw2,rp,rg,x,i,mh);
end
%% 
yp=x(end-n+1:end,1);
thetal=x(end-5000+1:end,5);
 figure(1);
 plot(thetal,x(end-5000+1:end,1));
 %% 刚度图
 figure(2);
 plot(thetal,k(end-5000+1:end));
 %% 啮合力
 figure(3);
 plot(thetal,F(end-5000+1:end));
%% 相图
 figure(4);
 plot(x(9216:end,1),x(9216:end,2)); 

