% This *.m file will be called by master file 'optimizamaster.m' for optimization. 
% This function will call *.m files of ODE model. It will also call 'gfit.m' to 
%determine the goodness-of-fit.

function [mle, XX] =optimizeprogITNTogo2(parameters,Inciddat,ITNdat,runningPlot,numinit)

global znew  TP  ITNP  XX  AA  pars


%sample size
nout = length(Inciddat);

% No. of objects to be estimated (beta, Initial Condition: S,L,I,A)
q = length(parameters);

%Initial Conditions
init=parameters(1:numinit);

% Vector of objects to be estimated
pars = parameters((numinit+1):q);
%display(pars, 'pars')

beta0 = pars(1)
r = pars(2)
cuppa = pars(3)
repos = pars(4)



%parameter values
p_h = 0.013;
b0 = 0.5;

 
% Setup the ODE solver
%time span
ts= 1:160;

options = odeset('RelTol',1e-4,'AbsTol',1e-4);


% Call the solver
z0=[init(1), init(2), init(3),init(4),init(5),init(6)];
[t, zout]=ode45(@(t,z)funmodelTogo(t,z,beta0,r,cuppa),ts,z0,options);
A=[t,zout];


zout(end,6)=0.057;
zout(end,2)=440;

znew = zout(end,:);


% ts1= 1:40;
ts1= 1:51;

[t, zout1]=ode45(@(t,z)funmodelTogo(t,z,beta0,r,cuppa),ts1,znew,options);
B=[t,zout1];

AA = [A;B];
Nh=AA(:,2)+AA(:,3)+AA(:,4);


    for j = 1:length(AA)
        TPP(j,1) =( ((p_h*beta0*(1-AA(j,7).*b0)).*AA(j,6)).*AA(j,2))./Nh(j,1);
        ITNPP(j,1) = 100*AA(j,7);
        
    end

 

  
%total population
%Totalpop = 106825181/2000;
TP = poissrnd(repos*(TPP)); %%BLACK%%
ITNP = poissrnd(0.72*(ITNPP)); %%BLACK%%




% mle = -sum(log(0.000001+poisspdf(Inciddat,TP(end-39:end,1))))-sum(log(0.000001+poisspdf(ITNdat,ITNP(end-39:end,1))));
% AIC = 2*q+2*mle;
% MLhood = mle ;
% 
% % measurement of Goodness of fit (by MSE)
% if runningPlot==1
% goodofit1 = gfit(Inciddat,TP(end-39:end,1),'9');
% goodofit2 = gfit(ITNdat,ITNP(end-39:end,1),'9');
% end
% 
% XX = [goodofit1,goodofit2];
% disp(XX)
% 
% if runningPlot==1
%     figure(1)
%     subplot(121),plot(1:nout,Inciddat,'k*-',1:nout,TP(end-39:end,1),'ro-');
%     subplot(122), plot(1:nout,ITNdat,'k*-',1:nout,ITNP(end-39:end,1),'g--');
% end 
%      
% end


% mle = -sum(log(0.0000001+poisspdf(Inciddat,TP(end-50:end-10,1))))-sum(log(0.0000001+poisspdf(ITNdat,ITNP(end-50:end-10,1))));
mle = -sum(log(0.0000001+poisspdf(Inciddat,TP(end-50:end-6,1))))-sum(log(0.0000001+poisspdf(ITNdat,ITNP(end-50:end-6,1))));

AIC = 2*q+2*mle;
MLhood = mle ;

% measurement of Goodness of fit (by MSE)
if runningPlot==1
% goodofit1 = gfit(Inciddat,TP(end-50:end-10,1),'9');
% goodofit2 = gfit(ITNdat,ITNP(end-50:end-10,1),'9');

goodofit1 = gfit(Inciddat,TP(end-50:end-6,1),'9');
goodofit2 = gfit(ITNdat,ITNP(end-50:end-6,1),'9');

end

XX = [goodofit1,goodofit2];
disp(XX)

if stoppingCondition(XX)

     assignin('base', 'TP', TP)
   assignin('base','ITNP',ITNP)
    assignin('base','XX',XX)
     
    assignin('base','AA',AA)
    assignin('base','pars',pars);
      
     error('Stopping criterion reached');
end


if runningPlot==1
    figure(1)
%     subplot(121),plot(1:nout,Inciddat,'k*-',1:nout,TP(end-50:end-10,1),'ro-');
%     subplot(122), plot(1:nout,ITNdat,'k*-',1:nout,ITNP(end-50:end-10,1),'g--');

subplot(121),plot(1:nout,Inciddat,'k*-',1:nout,TP(end-50:end-6,1),'ro-');
subplot(122), plot(1:nout,ITNdat,'k*-',1:nout,ITNP(end-50:end-6,1),'g--');

end 
     
end

% Stopping criterion function
function stop = stoppingCondition(XX)
    % Define your stopping criterion based on the intermediate values of x
    % Return true if the criterion is met, false otherwise
goodofit1=XX(1,1);
goodofit2=XX(1,2);

threshold1=0.4;
threshold2=0.5;

     if goodofit1 > threshold1 && goodofit2 > threshold2
%         disp('Stopping criterion reached');
        stop = true;
    else
        stop = false;
    end
   
end

