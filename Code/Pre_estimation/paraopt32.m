%  reparametrization in 6 month framework
function para = paraopt32()

% Lambda_h=10^3/(55*365);    % 10^3/(55*365) per day human birth rate
Lambda_h=(9.8*10^3*180)/(55*365);            % 10^3/(55) per year

% Lambda_v=10^4/14;    % 10^4/14 mosquito birth rate
Lambda_v=(10^4*180)/14;    % (10^4*365)/14 per year mosquito birth rate

% sigma_h=1/285;     % 1/285 (0.0035) per day rate at which infectious human acquire
                               % immunity [14/10^4 (0.0014) 17/10^3 (0.0170)] per day
                   
sigma_h=(180)*0.0049;     % immunity [14/10^4 (0.0014) 17/10^3 (0.0170)]*365 per year
% delta_h=32.9/(365*10^3);     % 32.9/(365*10^3)  [0 41/10^5(0.00014)] per day; disease induced mortality rate
delta_h=.000071*180;         %[0 41/10^5(0.00014)]*365 per year

% mu_h=1/(55*365);        % 1/(55*365) [1/72 1/35]*1/365 per day human natural death rate  
mu_h=180/(39*365);        % 1/(55*365) human natural death rate  

mu_v0=180*0.0476;       %1/14 per day mosquito natural death rate, [1/21 (0.0476), 1/14 (0.0714)]
mu_v1=180*0.0476;       % 1/14 per day dependent death rate of mosquito, [1/21 (0.0476), 5/10(0.2)]

% rho_h=1/(5*365);  % 1/(5*365) rate human loose malaria immunity
rho_h=.00040*180;  % [0.000055, .0110]*365 per year rate human loose malaria immunity
p_h=0.013;    % 22/10^3; [.01 .27] disease transmission probability from infectious human
             
p_v= .072;  % 48/10^2   [.072 0.64] disease transmission probability from infectious mosquitoes
             
cuppa=0.5;  % daily sampling rate of individuals per day
b_0=0.5;      % efficacy bed net

beta_0=200.6;   % average rate of mosquito biting rate per day, 0.30 for SO of nonauto model
w1 =1/1000;       % sensitivity parameter w.r.t infected humans
w2 =1/60000;       % 1/6000 sensitivity parameter w.r.t number of mosquitoes

r=0.08;
n=6;
 T=3*365;   %ITN Replacement period
% T=a;



 para(1)=Lambda_h;
  para(2)=Lambda_v;
  para(3)=sigma_h;
  para(4)= delta_h;
  para(5)=mu_h; 
  para(6)=mu_v0; 
  para(7)=mu_v1;  
  para(8)=rho_h;
  para(9)=p_h;  
  para(10)=p_v;
  para(11)=cuppa;
%    para(12)=Tau;
%    para(12)=L;
  para(13)=r;
%    para(15)=r_i;
  para(14)=b_0;
  para(15)=beta_0;
  para(16)=w1;
  para(17)=w2;
  %para(20)=b_1;
  para(18)=n;
  para(19)=T;
  
end
