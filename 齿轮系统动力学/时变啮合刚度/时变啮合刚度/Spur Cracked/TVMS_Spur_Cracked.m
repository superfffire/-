tic
%% Clear Workspace and Command Window
clc
clear
%% Enter Gear Values
m=2e-3;
N1=25;
N2=30;
alfa=20*pi/180;
h_a_star=1;
c_star=0.25;
E=206e9;
Nu=0.3;
L=20e-3;
r_int1=6.5e-3;%hubradius(m,N1,h_a_star,c_star);
r_int2=7.6e-3;%hubradius(m,N2,h_a_star,c_star);
%% Enter Crack Values
q_0=0.5e-3;
alfa_c=70*pi/180;
%% Calculate Initial Values
G=E/(2*(1+Nu));
r1=m*N1/2;
r_a1=r1+h_a_star*m;
r_f1=r1-h_a_star*m;
r_d1=r1-(c_star+h_a_star)*m;
r_b1=r1*cos(alfa);
r2=m*N2/2;
r_a2=r2+h_a_star*m;
r_f2=r2-h_a_star*m;
r_d2=r2-(c_star+h_a_star)*m;
r_b2=r2*cos(alfa);
eps_alfa=(sqrt(r_a2^2-r_b2^2)+sqrt(r_a1^2-r_b1^2)-(r1+r2)*sin(alfa))/(pi*m*cos(alfa));
inv_alfa=tan(alfa)-alfa;
teta_b1=pi/(2*N1)+inv_alfa;
teta_b2=pi/(2*N2)+inv_alfa;
h_f1=r_d1/r_int1;
h_f2=r_d2/r_int2;
%% Calculate Meshing Degrees
syms x
if r_b1<r_f1
    alfa_01=bsm2(r_b1*((x+teta_b1)*sin(x)+cos(x))-r_f1,0,pi/2,1e-4);
    alfa_12=bsm2(r_b2*((x+teta_b2)*sin(x)+cos(x))-r_a2,0,pi/2,1e-4);
else
    alfa_01=0;
    alfa_12=bsm2(r_b2*((x+teta_b2)*sin(x)+cos(x))-(r_a2-(r_b1-r_f1)),0,pi/2,1e-4);
end
if r_b2<r_f2
    alfa_02=bsm2(r_b2*((x+teta_b2)*sin(x)+cos(x))-r_f2,0,pi/2,1e-4);
    alfa_11=bsm2(r_b1*((x+teta_b1)*sin(x)+cos(x))-r_a1,0,pi/2,1e-4);
else
    alfa_02=0;
    alfa_11=bsm2(r_b1*((x+teta_b1)*sin(x)+cos(x))-(r_a1-(r_b2-r_f2)),0,pi/2,1e-4);
end
%% Calculate Root arcs and Degrees
syms x
if r_b1<r_d1
    beta_01=bsm2(r_b1*((x+teta_b1)*sin(x)+cos(x))-r_f1,0,pi/2,1e-4);
    L_d1=r_b1*((beta_01+teta_b1)*cos(beta_01)-sin(beta_01));
else
    beta_01=0;
    L_d1=r_b1*teta_b1;
end
if r_b2<r_d2
    beta_02=bsm2(r_b2*((x+teta_b2)*sin(x)+cos(x))-r_f2,0,pi/2,1e-4);
    L_d2=r_b2*((beta_02+teta_b2)*cos(beta_02)-sin(beta_02));
else
    beta_02=0;
    L_d2=r_b2*teta_b2;
end
teta_f1=atan(L_d1/r_d1);
teta_f2=atan(L_d2/r_d2);
S_f1=2*teta_f1*r_d1;
S_f2=2*teta_f2*r_d2;
%% Calculate Gear 1 Stiffnesses (Cracked)
if r_d1<=r_b1
    [K_a1c K_b1c K_s1c K_f1c]=toothmesh1(E,L,G,r1,r_b1,r_d1,teta_f1,S_f1,h_f1,teta_b1,alfa_01,alfa_11,beta_01,q_0,alfa_c);
else
    [K_a1c K_b1c K_s1c K_f1c]=toothmesh2(E,L,G,r1,r_b1,r_d1,teta_f1,S_f1,h_f1,teta_b1,alfa_01,alfa_11,beta_01,q_0,alfa_c);
end
%% Calculate Gear 1 Stiffnesses (Healthy)
if r_d1<=r_b1
    [K_a1 K_b1 K_s1 K_f1]=toothmesh1(E,L,G,r1,r_b1,r_d1,teta_f1,S_f1,h_f1,teta_b1,alfa_01,alfa_11,beta_01,0,0);
else
    [K_a1 K_b1 K_s1 K_f1]=toothmesh2(E,L,G,r1,r_b1,r_d1,teta_f1,S_f1,h_f1,teta_b1,alfa_01,alfa_11,beta_01,0,0);
end
%% Calculate Gear 2 Stiffnesses
if r_d2<=r_b2
    [K_a2 K_b2 K_s2 K_f2]=toothmesh1(E,L,G,r2,r_b2,r_d2,teta_f2,S_f2,h_f2,teta_b2,alfa_02,alfa_12,beta_02,0,0);
else
    [K_a2 K_b2 K_s2 K_f2]=toothmesh2(E,L,G,r2,r_b2,r_d2,teta_f2,S_f2,h_f2,teta_b2,alfa_02,alfa_12,beta_02,0,0);
end
K_a2=fliplr(K_a2);
K_b2=fliplr(K_b2);
K_s2=fliplr(K_s2);
K_f2=fliplr(K_f2);
%% Calculate Meshing Stiffness
K_A=1./(1./K_a1+1./K_a2);
K_B=1./(1./K_b1+1./K_b2);
K_F=1./(1./K_f1+1./K_f2);
K_S=1./(1./K_s1+1./K_s2);
K_h=(pi*E*L)/(4*(1-Nu^2));
K_h=ones(1,length(K_a1))*K_h;
K=1./(1./K_h + 1./K_A + 1./K_B + 1./K_S + 1./K_F);

K_Ac=1./(1./K_a1c+1./K_a2);
K_Bc=1./(1./K_b1c+1./K_b2);
K_Fc=1./(1./K_f1c+1./K_f2);
K_Sc=1./(1./K_s1c+1./K_s2);
Kc=1./(1./K_h + 1./K_Ac + 1./K_Bc + 1./K_Sc + 1./K_Fc);

PTH=length(K);
K2=[zeros(1,PTH) K zeros(1,PTH)];
K2c=[zeros(1,PTH) Kc zeros(1,PTH)];
Pb=floor(PTH/eps_alfa);
dbeta=360*eps_alfa/(N1*(PTH-1));
for i=1:PTH
    K_M(i)=K2c(i+PTH)+K2(i+PTH+Pb)+K2(i+PTH-Pb);
    betai(i)=dbeta*(i-1);
end
%% Plot Figures
figure(1)
plot(betai,K)
figure(2)
plot(betai,K_M)
toc