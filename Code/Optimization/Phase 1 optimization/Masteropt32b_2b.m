
function [B,C]=Masteropt32b_2b(m, p, q )

T1 = 52;

tau1=1:T1;

 OptPvalue2 = paraValues32b(p,q);

p_h=OptPvalue2(9) ;
 b_0= OptPvalue2(11);
beta_0=OptPvalue2(12); 

 options = odeset('RelTol',1e-4,'AbsTol',1e-4);
[t, B] = ode45(@modelopt32b_2b,tau1, m,options,OptPvalue2);
  
S1sol=B(:,1);
I1sol=B(:,2);
R1sol=B(:,3);
S2sol=B(:,4);
I2sol=B(:,5);
Xsol=B(:,6);
% 
 N1sol=S1sol+I1sol+R1sol;

k=41;
for i=1:T1
if i<=k
bbeta(i)=0.5;
else
    bbeta(i)=b_0;
end
 beta(i)=beta_0*(1-Xsol(i)*bbeta(i));
 
 D(i)=((p_h*(beta(i)*I2sol(i))*S1sol(i))/N1sol(i));
end

C=D';



 end

