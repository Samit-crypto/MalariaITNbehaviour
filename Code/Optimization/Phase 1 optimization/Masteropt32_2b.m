clear all
close all

[A1,B1]=Masteropt32a_2();

upre=A1(end,:);
upre(:,6)=[];
x0=0.057;   % ITN usage in 2000
uinitial=[upre x0];

load('Togo_2new.mat')
pars1(1,2)
pars1(1,3);
n=100;

A=lhsdesign(n,2);
A(:,1)=pars1(1,3)+(5-(pars1(1,3)))*A(:,1);   % for cuppa [0.55 10]
A(:,2)=0 + (pars1(1,2))*A(:,2);  % for r    [0 0.1]
A(101,:) = [pars1(1,3), pars1(1,2)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% 
for j=1:n+1
 [uj,Vj] = Masteropt32b_2b(uinitial, A(j,1), A(j,2));

M1(:,j)=Vj;
M2(:,j)=100*uj(:,end);

figure(121)
subplot(121)
if j == n+1
    plot(Vj,'ro')
    hold on
else
plot(Vj)
hold on
end

subplot(122)
if j == n+1
 plot(M2(:,j),'go')
 hold on
else
plot(M2(:,j))
hold on
end

J1(:,j)=M1(45:51,j);
J2(:,j)=M2(45:51,j);

for k=1:100
F1(:,j,k)=poissrnd(pars1(1,4)*J1(:,j));%prob put from Spread sheet 
F2(:,j,k)=poissrnd(0.72*J2(:,j));%take from optimizeITN Co.m
end   

 DataMean1_1(:,j)=mean(F1(:,j,k),3);

eta1=0.2;
eta2=0.1;
eta3=0.2;
Q1(1,j)=trapz(eta1*DataMean1_1(:,j)+eta2*abs( A(j,1)-pars1(1,3))+eta3*abs(A(j,2)-pars1(1,2)));

end

[val1, idx1] = min(Q1);

Malcas=352.5;
Emalcas2020=Malcas-(Malcas*(.40))
 EMalcas2025 = Malcas-(Malcas*(.75))
AchieveMalcas2020=236.2;							
   AchieveITNUse2020=80.8;

 W(1,:)=J1(end,:);
 ind = interp1(W,1:length(W),EMalcas2025,'nearest','extrap');   % find index position in (cuppa, r) less than and closet to Expected malaria incidence cases 
 Q1_1=Q1(1,idx1)
 Q1_2=Q1(1,ind)
 Q1_3=Q1(1,101)
 AAA=A(idx1,:)
 BBB=A(ind,:)
 CCC=A(101,:)

 Y11_1=TP(161:211);
 Y22_2=ITNP(161:211);
 all_pointsInci =Togo_extractInci(1:2:end,2);
  all_pointsITN = Togo_extractITN(1:2:end,2);
 NewY11_1=Y11_1(1:2:end);
  NewY22_2=Y22_2(1:2:end);


 TP(211)
EMalcas2025

 figure(222)
%  subplot(121)
 plot(M1(:,101))
 subplot(122)
 plot(M2(:,101))

AAA=A(idx1,:);
  J11(:,1)=M1(45:51,idx1);
   J22(:,1)=M2(45:51,idx1);
for k=1:100
F11(:,k)=poissrnd(pars1(1,4)*J11);
F22(:,k)=poissrnd(0.72*J22);
end   
% % % 
 O=cumsum(F11,1);
 O1=cumsum(J11,1);
 MeanInc= O1(end);
 O2=O(end,:);
 
figure(987)
histogram(O2,10)
hold on 
xline(MeanInc);
% % % 
X1=1:7;
 X2=1:7;
 Data1=F11;
 Data2=F22;
N1=size(Data1,1);
N2=size(Data2,1);
 DataMean1=mean(Data1,2);
 DataMean2=mean(Data2,2);
 DataSEM1=std(Data1,0,2)/sqrt(N1);
 DataSEM2=std(Data2,0,2)/sqrt(N2);
 CI951=tinv([0.025 0.975],N1-1);
 CI952=tinv([0.025 0.975],N2-1);
% % DataCI95=bsxfun(@times,DataSEM,CI95(:));
% % CI95=tinv(0.95,N-1);
  DataCI951=bsxfun(@times,DataSEM1,CI951.*ones(7,1));
  DataCI952=bsxfun(@times,DataSEM2,CI952.*ones(7,1));
 CIdist1=abs(DataCI951-DataMean1);
 CIdist2=abs(DataCI952-DataMean2);
% 
Y1=CIdist1(:,1);
Z1=CIdist1(:,2);
Y2=CIdist2(:,1);
Z2=CIdist2(:,2);

OptMalcas2035=DataMean1(end)

% % 
figure(888)
subplot(121)
plot(X1, DataMean1','k','LineWidth', 2) 
hold on
plot(X1, Y1', 'b', 'LineWidth', 2);
hold on;
plot(X1, Z1', 'b', 'LineWidth', 2);
X11 = [X1, fliplr(X1)];
inBetween = [Y1', fliplr(Z1')];
fill(X11, inBetween, 'g');
subplot(122) 
plot(X2, DataMean2','k','LineWidth', 2) 
hold on
plot(X2, Y2', 'b', 'LineWidth', 2);
hold on;
plot(X2, Z2', 'b', 'LineWidth', 2);
X22 = [X2, fliplr(X2)];
inBetween = [Y2', fliplr(Z2')];
fill(X22, inBetween, 'g');


 load CI_BootStrap_Togo2.mat;


max_values_Inci=max(subsetInci, [], 2);
min_values_Inci=min(subsetInci, [], 2);
max_values_ITN=max(subsetITN, [], 2);
min_values_ITN=min(subsetITN, [], 2);


 figure(111)
subplot(121)
plot(1:45, max_values_Inci','Color', [0.65 0.65 0.65])
hold on
plot(1:45, min_values_Inci','Color', [0.65 0.65 0.65])
hold on

X_A = [1:45, fliplr(1:45)];
inBetween = [min_values_Inci', fliplr(max_values_Inci')];
c=[0.9 0.9 0.9];
fill(X_A, inBetween, c);

plot(1:2:length(Y11_1(:,1)), NewY11_1, 'k');
hold on
plot(45:51,Y1)
hold on
plot(45:51,Z1)

scatter(1:2:length(Togo_extractInci(:,2)), all_pointsInci, 'k', 'filled', 'Marker', 'o');
%plot(all_pointsInci)
hold on 
xline(41)
X111 = [45:51, fliplr(45:51)];
inBetween = [Y1', fliplr(Z1')];
fill(X111, inBetween, 'g');
hold on
plot(45:51,DataMean1)

P1=31;
P2=Malcas;
plot(P1,P2,'r*')
text(P1,P2,['(', num2str(P1), ', ', num2str(P2), ')'])
hold on
P3=Emalcas2020;
P4=41;
text(P4,P3,['(', num2str(P4), ', ', num2str(P3), ')'])
hold on
plot(P4,P3,'b*')
text(P4,P3,['(', num2str(P4), ', ', num2str(P3), ')'])
hold on
scatter(P4, AchieveMalcas2020, 'k', 'filled', 'Marker', 'o')
hold on
scatter(43, AchieveMalcas2021, 'k', 'filled', 'Marker', 'o')
hold on
scatter(45, AchieveMalcas2022, 'k', 'filled', 'Marker', 'o')
hold on

P5=51;
P6=EMalcas2025;
plot(P5,P6,'b*')
text(P5,P6,['(', num2str(P5), ', ', num2str(P6), ')'])


subplot(122)
plot(1:45, max_values_ITN','Color', [0.65 0.65 0.65])
hold on
plot(1:45, min_values_ITN','Color', [0.65 0.65 0.65])
hold on

X_B = [1:45, fliplr(1:45)];
inBetween = [min_values_ITN', fliplr(max_values_ITN')];
c=[0.9 0.9 0.9];
fill(X_B, inBetween, c);


plot(1:2:length(Y22_2(:,1)), NewY22_2, 'k');
 hold on
 plot(45:51,Y2)
hold on
plot(45:51,Z2)
hold on
scatter(1:2:length(Togo_extractITN(:,2)), all_pointsITN, 'r', 'filled', 'Marker', 'o');

hold on 
scatter(41,  AchieveITNUse2020 ,'r', 'filled', 'Marker', 'o')
hold on 
scatter(43,  AchieveITNUse2021 ,'r', 'filled', 'Marker', 'o')
hold on 
scatter(45,  AchieveITNUse2022 ,'r', 'filled', 'Marker', 'o')
xline(41)
X222 = [45:51, fliplr(45:51)];
inBetween = [Y2', fliplr(Z2')];
fill(X222, inBetween, 'g');
hold on
plot(45:51, DataMean2)
