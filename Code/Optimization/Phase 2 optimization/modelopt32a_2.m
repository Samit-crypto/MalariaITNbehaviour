function udash= modelopt32a_2(t,u,OptPvalue1)

  S1=u(1);
  I1=u(2);
  R1=u(3);
  S2=u(4);
  I2=u(5);
  X=u(6);

Lambda_h = OptPvalue1(1) ;
 Lambda_v = OptPvalue1(2) ;
 sigma_h = OptPvalue1(3) ;
 delta_h = OptPvalue1(4) ;
 mu_h = OptPvalue1(5) ;
mu_v0 = OptPvalue1(6) ;
 mu_v1 = OptPvalue1(7) ;
 rho_h = OptPvalue1(8) ;
 p_h = OptPvalue1(9) ;
 p_v = OptPvalue1(10) ;
b_0 = OptPvalue1(11) ;
 beta_0 = OptPvalue1(12);
w1 = OptPvalue1(13);
 w2 = OptPvalue1(14);
 c1 = OptPvalue1(15);
 c2 = OptPvalue1(16);
 Malcas = OptPvalue1(17);
cuppa = OptPvalue1(18);
r = OptPvalue1(19);
  
   
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


 udash=[S1dash;I1dash;R1dash;S2dash;I2dash;Xdash];
 end