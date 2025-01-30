%% Bootstrap method for estimation of parameters and Confidence interval (Master file)
% This will call other *.m files to estimate the parameters and complete
% the bootstrap method for estimating CI of these parameters. 
% Written by Samit Bhattacharyya, University of Guelph, August 2009
%---------------------------------------------------------------------
function [TP, ITNP, XX, pars, AA]=optimizemastITNTogo2(beta0, r, cuppa, Incid, ITN)

% tic;


%%

 repos = 1;


IniPts = [1000 100 0 2500 1000 0.0];

%Start values of Initial condition
Init1S= IniPts(1);
Init2S= IniPts(2);
Init3S= IniPts(3);
Init4S= IniPts(4);
Init5S= IniPts(5);
Init6S= IniPts(6);


% lb = [ 999, 99,   0, 2499, 999, 0.0, 139,   0.27,  0.75,  0];    % for 6 month parameter regime
% ub = [1001,101, 0.1, 2510,1010, 0.0, 140,  0.28,   0.8,  1];

lb = [ 999, 99,   0, 2499, 999, 0.0, beta0-1,  r-0.01,  cuppa-0.05,  0];    % for 6 month parameter regime
ub = [1001,101, 0.1, 2510,1010, 0.0, beta0+1,  r+0.01,   cuppa+0.05,  1];



startvals = [Init1S, Init2S, Init3S, Init4S, Init5S, Init6S, beta0, r,cuppa, repos]; 


%%

numinit = 6;         

% Parameter names
parnames = {'S_h'; 'I_h'; 'R_h';'S_v'; 'I_v'; 'X_h'; 'beta0'; 'r'; 'cuppa'};          

% Number of running plot
runningPlot = 1;     


%-------------------------------------------------------------------------
% OPTIMIZATION for parameter estimation
 
 opts = optimoptions('fmincon','Display','iter','Algorithm','active-set','TolX',1e-100000,'TolCon',1e-10000,'TolFun',1e-100000,...
     'MaxIter',20, 'MaxFunEvals',6000);


%opts = optimset('Display','iter','Algorithm','SQP','MaxIter',40000); % Sets fminsearch options
% [pEstimates, ~, ~, output, ~, ~, HESSIAN] = fmincon(Objfcn,startvals,[],[],[],[],lb,...
%     ub,[],opts); %   as starting values
 


[x, fvalOpt, exitflag, output,~,~,HESSIAN]=fmincon(@(par)...        % Computes model parameter
    optimizeprogITNTogowrapper(par,Incid,ITN,runningPlot,...    %   estimates using startvals
    numinit),startvals,[],[],[],[],lb,...
    ub,[],opts) %   as starting values


    function mle = optimizeprogITNTogowrapper(startvals,Incid,ITN,runningPlot,...    %   estimates using startvals
    numinit)
    % This is a wrapper function that calls your actual objective function
    % It extracts the first output and returns it as the fval
    
    [mle, ~] = optimizeprogITNTogo2(startvals,Incid,ITN,runningPlot,numinit);
end


   
%--------------------------------------------------------------------------
% optMast = struct('TP',TP,'ITNP',ITNP,'GOF', goodofit,'estimates',pEstimates,'Hessian',HESSIAN);
% 
% 
%  save('TP',TP,'ITNP',ITNP,'XX',XX,'pars',pars,'AA',AA)
end
% toc;