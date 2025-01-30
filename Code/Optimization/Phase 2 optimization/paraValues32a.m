%Rep Congo

function OptPvalue1 = paraValues32a()

Lambda_h=(10.2*10^3*180)/(55*365);      % 10^3/(55*365) human birth rate
Lambda_v=(10^4*180)/14;    % 10^4/14 mosquito birth rate
sigma_h=(180)*0.0049;     % 1/285  [.0014  .0170] rate at which infectious human acquire immunity
delta_h=.000071*180;     % 32.9/(365*10^3)  [0 41/10^5(0.00014)]; disease induced mortality rate
mu_h=180/(39*365);       % 1/(55*365) [1/72*365 1/35*365]human natural death rate  
mu_v0=180*0.0476;      %1/14 [1/21 1/14] mosquito natural death rate, 1/64 for SO of nonauto model
mu_v1=180*0.0476;      % 1/14  dependent death rate of mosquito, 0.21 for SO of nonauto model
rho_h=.00040*180;  % 1/(5*365) [0.000055, .0110]rate human loose malaria immunity
p_h=0.013;    % 22/10^3  [.01 .27]; disease transmission probability from infectious human
              % to susceptible mosquitoes, .01 for SO of nonauto model
p_v=0.072 ;  % 48/10^2 [.072 0.64]disease transmission probability from infectious mosquitoes
              % to susceptible humans, .10 for SO of nonauto model
 cuppa=0.55;  % daily sampling rate of individuals per day

b_0=0.5;      % efficacy bed net
 beta_0=215.6;   % average rate of mosquito biting rate per day, 0.30 for SO of nonauto model
w1 =1/1000;       % sensitivity parameter w.r.t infected humans
w2 =1/60000;       % 1/6000 sensitivity parameter w.r.t number of mosquitoes
r=0.1;
% c1=0.2;
% c2=0.1;
Malcas=352.5;


OptPvalue1(1) = Lambda_h;
OptPvalue1(2) = Lambda_v;
OptPvalue1(3) = sigma_h;
OptPvalue1(4) = delta_h;
OptPvalue1(5) = mu_h;
OptPvalue1(6) = mu_v0;
OptPvalue1(7) = mu_v1;
OptPvalue1(8) = rho_h;
OptPvalue1(9) = p_h;
OptPvalue1(10) = p_v;
OptPvalue1(11) = b_0;
OptPvalue1(12) = beta_0;
OptPvalue1(13) = w1;
OptPvalue1(14) = w2;
% OptPvalue1(15) = c1;
% OptPvalue1(16) = c2;
OptPvalue1(17) = Malcas;
OptPvalue1(18) = cuppa;
OptPvalue1(19) = r;


