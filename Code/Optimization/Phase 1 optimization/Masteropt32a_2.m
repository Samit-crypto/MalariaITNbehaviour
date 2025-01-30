function [A1,B1]=Masteropt32a_2()

T1 = 160;
% tau=[0:T1];


tau=1:T1;

OptPvalue1 = paraValues32a();
% para=paraopt38();
u0= [1000 100 0 2500 1000 0];


 p_h=OptPvalue1(9) ;
 b_0= OptPvalue1(11);
beta_0=OptPvalue1(12); 

 options = odeset('RelTol',1e-4,'AbsTol',1e-4);
[t,A1] = ode45(@modelopt32a_2,tau,u0,options,OptPvalue1);
  
S1sol=A1(:,1);
I1sol=A1(:,2);
R1sol=A1(:,3);
S2sol=A1(:,4);
I2sol=A1(:,5);
Xsol=A1(:,6);
% 
 N1sol=S1sol+I1sol+R1sol;

for i=1:T1
 bbeta(i)=b_0;
 beta(i)=beta_0*(1-Xsol(i)*bbeta(i));
 
 D(i)=((p_h*(beta(i)*I2sol(i))*S1sol(i))/N1sol(i));

end

B1=D';

 end

