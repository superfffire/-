function [km,F]  = shuntaigangdu(kx2,zp,kw2,rp,rg,x,i,mh)
%SHUNTAIGANGDU 此处显示有关此函数的摘要
%   此处显示详细说明
%% 齿轮副瞬态刚度
kwx1=mod(x(i,5),2*pi/zp);
%实时刚度计算
km=interp1(kx2,kw2,kwx1,'makima');
km
%% 轴承1刚度
% wtheta2=0*pi/180;
% nbearing=w1/2/pi;
% fbearing=nbearing*0.5*nball*(1-dbearing/Dbearing);
% Tbearing=1/fbearing;
%kpy=0.075*10^8*sin(2*pi/Tbearing*t);
kpy=0.075*10^8;
%% 轴承2刚度
% wtheta2=0*pi/180;
% nbearing=w1/2/pi*zp/zg;
% fbearing=nbearing*0.5*nball*(1-dbearing/Dbearing);
% Tbearing=1/fbearing;
%kgy=0.075*10^8*sin(2*pi/Tbearing*t);
kgy=0.075*10^8;
%% 阻尼
alphac=0.0000625;
bate=0.0000008;
cm=alphac*mh+bate*km;
%% 误差
e=0;
edot=0;
%% 动态啮合力
Fk=km.*(rp.*x(i,5)-x(i,1)+x(i,3)-rg.*x(i,7)-e);
Fc=cm.*(rp.*x(i,6)-x(i,2)+x(i,4)-rg.*x(i,8)-edot);
F=Fk+Fc;
end

