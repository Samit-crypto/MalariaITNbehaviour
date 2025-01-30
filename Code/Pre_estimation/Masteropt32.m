% 
clear all
close all


[u11,V1]=Masteropt32a();

H1=u11(:,7);
upre=u11(end,:);
upre(:,7)=[];
upre(:,6)=[];
x0=0.057;     % ITN usage in 2000
uinitial=[upre x0];
% 
[u22,V2]=Masteropt32b(uinitial);

unew=[u11; u22];
vnew=[V1;V2];

load Togo_extractDeath.mat;
F3= Togo_extractDeath(:,2);

 load Togo_extractInci.mat;
  F1=Togo_extractInci(:,2);
 load Togo_extractITN.mat;
 F2= Togo_extractITN(:,2);

figure(1111)
plot(1:20,F3)
hold on
plot(H1(end-19:end))

 
 figure(102)

subplot(121), plot(unew(end-39:end,2))
hold on 
plot(1:40,F1)

subplot(122), plot(100*unew(end-39:end,6))
hold on
plot(1:40,F2)