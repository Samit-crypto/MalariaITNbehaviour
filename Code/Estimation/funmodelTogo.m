function [zdot]=funmodelTogo(t,z,p1,p2,p3)
                 

 
beta0=p1;   % average rate of mosquito biting rate per day, 0.30 for SO of nonauto model
r=p2;
cuppa=p3;
% w1=p4;
% w2=p5;

Lambda_h=(10.2*10^3*180)/(55*365);     % 10^3/(55*365) human birth rate
Lambda_v=(10^4*180)/14;    % 10^4/14 mosquito birth rate
sigma_h=(180)*0.0049;    % 1/285  [.0014  .0170] rate at which infectious human acquire immunity
delta_h=.000071*180;     % 32.9/(365*10^3)  [0 41/10^5(0.00014)]; disease induced mortality rate
mu_h=180/(39*365);        % 1/(55*365) [1/72*365 1/35*365]human natural death rate  
mu_v0=180*0.0476;       %1/14 [1/21 1/14] mosquito natural death rate, 1/64 for SO of nonauto model
mu_v1=180*0.0476;       % 1/14  dependent death rate of mosquito, 0.21 for SO of nonauto model
rho_h=.00040*180;  % 1/(5*365) [0.000055, .0110]rate human loose malaria immunity
p_h=0.013;    % 22/10^3  [.01 .27]; disease transmission probability from infectious human
              % to susceptible mosquitoes, .01 for SO of nonauto model
p_v=0.072 ;  % 48/10^2 [.072 0.64]disease transmission probability from infectious mosquitoes
              % to susceptible humans, .10 for SO of nonauto model
% cuppa=0.01;  % daily sampling rate of individuals per day

b_0=0.5;      % efficacy bed net
% beta_0=0.99;   % average rate of mosquito biting rate per day, 0.30 for SO of nonauto model
w1 =1/1000;       % sensitivity parameter w.r.t infected humans
w2 =1/60000;       % 1/6000 sensitivity parameter w.r.t number of mosquitoes
% r=2.5;


bbeta=b_0;
bmuv=b_0;

N1=z(1)+z(2)+z(3);

beta_0=beta0*(1-z(6)*bbeta);
lambda_h=(p_h*beta_0*z(5))/N1;
lambda_v=(p_v*beta_0*z(2))/N1;
Xi_h=delta_h+sigma_h+mu_h;

mu_v=mu_v0+mu_v1*z(6)*bmuv;
N2=z(4)+z(5);

    zdot = [Lambda_h-lambda_h*z(1)+rho_h*z(3)-mu_h*z(1);
             lambda_h*z(1)-Xi_h*z(2);
            sigma_h*z(2)-(rho_h+mu_h)*z(3);
             Lambda_v-lambda_v*z(4)-mu_v*z(4);
             lambda_v*z(4)-mu_v*z(5);
             cuppa*z(6)*(1-z(6))*(-r+(1-bbeta)*(w1*z(2)+w2*N2))];
      