function [u11,V1]=Masteropt32a()

T1 = 160;
% tau=[0:T1];


tau=1:T1;


para=paraopt32();
u0= [1000 100 0 2500 1000 0];

delta_h=para(4);
p_h=para(9);
 b_0= para(14);
beta_0=para(15); 

 options = odeset('RelTol',1e-4,'AbsTol',1e-4);
[t,u1] = ode45(@modelopt32a,tau, u0,options,para);
  
S1sol=u1(:,1);
I1sol=u1(:,2);
R1sol=u1(:,3);
S2sol=u1(:,4);
I2sol=u1(:,5);
Xsol=u1(:,6);
% 
 N1sol=S1sol+I1sol+R1sol;

for i=1:T1
 bbeta(i)=b_0;
 beta(i)=beta_0*(1-Xsol(i)*bbeta(i));
 
 D(i)=((p_h*(beta(i)*I2sol(i))*S1sol(i))/N1sol(i));

end
B=(delta_h*u1(:,2)*1000)./N1sol;
 u11=[u1 B];
V1=D';

 end

