function zdash= modelopt32bb_2(t,z,OptPvalue2)


  S1=z(1);
  I1=z(2);
  R1=z(3);
  S2=z(4);
  I2=z(5);
  X=z(6);

Lambda_h = OptPvalue2(1) ;
 Lambda_v = OptPvalue2(2) ;
 sigma_h = OptPvalue2(3) ;
 delta_h = OptPvalue2(4) ;
 mu_h = OptPvalue2(5) ;
mu_v0 = OptPvalue2(6) ;
 mu_v1 = OptPvalue2(7) ;
 rho_h = OptPvalue2(8) ;
 p_h = OptPvalue2(9) ;
 p_v = OptPvalue2(10) ;
b_0 = OptPvalue2(11) ;
 beta_0 = OptPvalue2(12);
w1 = OptPvalue2(13);
 w2 = OptPvalue2(14);
 c1 = OptPvalue2(15);
 c2 = OptPvalue2(16);
 Malcas = OptPvalue2(17);
cuppa = OptPvalue2(18);
r = OptPvalue2(19);

loadedData1=load("optpara1new.mat");
c=loadedData1.AAA(1,1);
d=loadedData1.AAA(1,2);

loadedData=load('Togo_2new.mat');

A = loadedData.pars1(1,3);
B = loadedData.pars1(1,2);


if t<=51
    if t<=45
cuppa = A;
r = B;
b_0=0.5;
    else
cuppa = c;
r = d;
b_0=0.5;
    end 
end


bbeta=b_0;
 bmuv=b_0;

N1=S1+I1+R1;

beta=beta_0*(1-X*bbeta);
lambda_h=(p_h*beta*I2)/N1;
lambda_v=(p_v*beta*I1)/N1;
Xi_h=delta_h+sigma_h+mu_h;

mu_v=mu_v0+mu_v1*X*bmuv;
N2=S2+I2;

% [t,(-r+(b_0)*(w1*I1+w2*N2))];
% Game theoretic model of ITN
S1dash = Lambda_h-lambda_h*S1+rho_h*R1-mu_h*S1;
I1dash = lambda_h*S1-Xi_h*I1;
R1dash = sigma_h*I1-(rho_h+mu_h)*R1;
S2dash = Lambda_v-lambda_v*S2-mu_v*S2;
I2dash = lambda_v*S2-mu_v*I2; 
Xdash = cuppa*X*(1-X)*(-r+(1-bbeta)*(w1*I1+w2*N2));


 zdash=[S1dash;I1dash;R1dash;S2dash;I2dash;Xdash];
 end