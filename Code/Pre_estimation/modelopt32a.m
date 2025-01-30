function udash= modelopt32a(t,u,para)

  S1=u(1);
  I1=u(2);
  R1=u(3);
  S2=u(4);
  I2=u(5);
  X=u(6);

Lambda_h=para(1);
 Lambda_v= para(2);
  sigma_h=para(3);
   delta_h=para(4);
  mu_h=para(5); 
 mu_v0= para(6); 
  mu_v1=para(7);  
  rho_h=para(8);
  p_h=para(9);  
  p_v=para(10);
  cuppa=para(11);

  r=para(13);
%    r_i=para(15);
  b_0=para(14);
  beta_0=para(15);
  w1=para(16);
  w2=para(17);
 
   n= para(18);
  T=para(19);

  
  


%     bbeta = (((2^n+1)/(2^(n+1)))*(((2^n-1)/(2^n+1))+(1/(1+(mod(t,T)/(T/2))^n)))*b_0);
               bbeta=b_0;
%       bmuv = (((2^n)+1/2^(n))*(((-1)/((2^n)+1))+(1/(1+(mod(t,T)/(T/2))^n)))*b_0);
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
% S2dash = Lambda_v-lambda_v*S2-mu_v*S2-a*b_1*S2;
% I2dash = lambda_v*S2-mu_v*I2-a*b_1*I2; 

S2dash = Lambda_v-lambda_v*S2-mu_v*S2;
I2dash = lambda_v*S2-mu_v*I2; 
Xdash = cuppa*X*(1-X)*(-r+(1-bbeta)*(w1*I1+w2*N2));


 udash=[S1dash;I1dash;R1dash;S2dash;I2dash;Xdash];
 end