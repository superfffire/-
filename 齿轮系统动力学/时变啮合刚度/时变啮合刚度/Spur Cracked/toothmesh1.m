function [K_a K_b K_s K_f]=toothmesh1(E,L,G,r,r_b,r_d,teta_f,S_f,h_f,teta_b,alfa_0,alfa_1,beta_0,q_0,alfa_c)
C_pcf = [ -5.574e-5 -1.9986e-3 -2.3015e-4  4.7702e-3  0.0271  6.8045;
          60.111e-5  28.100e-3 -83.431e-4 -9.9256e-3  0.1624  0.9086;
         -50.952e-5  185.50e-3  0.0538e-4  53.300e-3  0.2895  0.9236;
         -6.2042e-5  9.0889e-3 -4.0964e-4  7.8297e-3 -0.1472  0.6904];
Ai=C_pcf(:,1);    Bi=C_pcf(:,2);    Ci=C_pcf(:,3);
Di=C_pcf(:,4);    Ei=C_pcf(:,5);    Fi=C_pcf(:,6);
for i=1:4
    Xi(i)=Ai(i)/teta_f^2 + Bi(i)*h_f^2 + Ci(i)*h_f/teta_f + Di(i)/teta_f + Ei(i)*h_f + Fi(i);
end
L_star=Xi(1);    M_star=Xi(2);    P_star=Xi(3);    Q_star=Xi(4);
h0=r_b*teta_b;
h_q=h0-q_0*sin(alfa_c);
dbeta=(alfa_1-alfa_0)/1000;
for i=1:1001
    beta=(i-1)*dbeta+alfa_0;
    h=r_b*((beta+teta_b)*cos(beta)-sin(beta));
    d=r_b*((beta+teta_b)*sin(beta)+cos(beta))-r;
    u_f=d+r-r_d-h*tan(beta);
    invK_f=(cos(beta)^2/(E*L))*(L_star*(u_f/S_f)^2 + M_star*(u_f/S_f) + P_star*(1+Q_star*(tan(beta))^2));
    invK_a=0;
    invK_b=0;
    invK_s=0;
    dbetaj=(beta-beta_0)/1000;
    for j=1:1000
        beta1=(j-1)*dbetaj+beta_0;
        h1=r_b*((beta1+teta_b)*cos(beta1)-sin(beta1));
        x1=r_b*((beta1+teta_b)*sin(beta1)+cos(beta1))-r;
        beta2=(j)*dbetaj+beta_0;
        h2=r_b*((beta2+teta_b)*cos(beta2)-sin(beta2));
        x2=r_b*((beta2+teta_b)*sin(beta2)+cos(beta2))-r;
        dx=x2-x1;
        if h2>h_q
            h2=h_q;
        end
        A_x=(h1+h2)*L;
        I_x=1/12*(h1+h2)^3*L;
        invK_a=invK_a + (sin(beta)^2/(E*A_x))*dx;
        invK_b=invK_b + (((d-x1)*cos(beta)-h*sin(beta))^2/(E*I_x))*dx;
        invK_s=invK_s + (1.2*cos(beta)^2/(G*A_x))*dx;
    end
    dx=(r_b-r_d)/1000;
    for j=1:1000
        x1=r_d+(j-1)*dx-r;
        h1=h0;
        x2=r_d+(j)*dx-r;
        h2=h0;
        dx=x2-x1;
        if h2>h_q
            h2=h_q;
        end
        A_x=(h1+h2)*L;
        I_x=1/12*(h1+h2)^3*L;
        invK_a=invK_a + (sin(beta)^2/(E*A_x))*dx;
        invK_b=invK_b + (((d-x1)*cos(beta)-h*sin(beta))^2/(E*I_x))*dx;
        invK_s=invK_s + (1.2*cos(beta)^2/(G*A_x))*dx;
    end
    K_a(i)=1/invK_a;
    K_b(i)=1/invK_b;
    K_s(i)=1/invK_s;
    K_f(i)=1/invK_f;
end
end