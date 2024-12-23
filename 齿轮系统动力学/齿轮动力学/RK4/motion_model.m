function dx=motion_model(t,x,~,zp,kw2,mp,mg,Ip,Ig,rg,rp,mh,Tg,Tp,kx2)
%% 齿轮副瞬态刚度
kwx1=mod(x(5),2*pi/zp);

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
alpha=0.02;%阻尼比
cpy=alpha*sqrt(kpy/mp);
cgy=alpha*sqrt(kgy/mg);
%%
e=0;
edot=0;
%% F
Fk=km.*(rp.*x(5)-x(1)+x(3)-rg.*x(7)-e);
Fc=cm.*(rp.*x(6)-x(2)+x(4)-rg.*x(8)-edot);
F=Fk+Fc;
%% 动力学方程
dx = [x(2);
    (F-cpy.*x(2)-kpy.*x(1))./mp;
     x(4);
    (-F-cgy.*x(4)-kgy.*x(3))./mg;
     x(6);
    (Tp-rp.*F)./Ip;
     x(8);
    (rg.*F-Tg)./Ig
    ];