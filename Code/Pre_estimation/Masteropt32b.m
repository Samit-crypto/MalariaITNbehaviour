function [u22,V2]=Masteropt32b(m)

T1 = 40;
% tau=[0:T1];


tau=1:T1;


para=paraopt32();

delta_h=para(4);
p_h=para(9);
 b_0= para(14);
beta_0=para(15); 

 options = odeset('RelTol',1e-4,'AbsTol',1e-4);
[t,u2] = ode45(@modelopt32a,tau, m,options,para);
  
S1sol=u2(:,1);
I1sol=u2(:,2);
R1sol=u2(:,3);
S2sol=u2(:,4);
I2sol=u2(:,5);
Xsol=u2(:,6);
% 
 N1sol=S1sol+I1sol+R1sol;

for i=1:T1
 bbeta(i)=b_0;
 beta(i)=beta_0*(1-Xsol(i)*bbeta(i));
 
 D(i)=((p_h*(beta(i)*I2sol(i))*S1sol(i))/N1sol(i));

end
B=(delta_h*u2(:,2)*1000)./N1sol;
 u22=[u2 B];
V2=D';

 end

