clear all
close all

[C1,D1]=Masteropt32a_2();

upre1=C1(end,:);
upre1(:,6)=[];
x0=0.057;     % ITN usage in 2000
uinitial1=[upre1 x0];
% 
load("optpara1new.mat");
load('optphase1new.mat')
AA1=AAA(1,1);
AA2=AAA(1,2);
n=10;

B=lhsdesign(n,3);
B(:,1)=AA1+(5-AA1)*B(:,1);     % for cuppa [AA1 10]
B(:,2)=0+(AA2)*B(:,2);     % for r    [0 AA2]
B(:,3)=0.5+(0.1)*B(:,3);    % for b_0  [0.5  0.6]
B(11,:) = [AA1,AA2,0.5];

for j=1:n+1
[u1j,V1j] = Masteropt32bb_2(uinitial1, B(j,1), B(j,2), B(j,3));

M11(:,j)=V1j;
M22(:,j)=100*u1j(:,end);

figure(121)
subplot(121)
if j == n+1
    plot(V1j,'ro')
    hold on
else
plot(V1j)
hold on
end

subplot(122)
if j == n+1
 plot(M22(:,j),'go')
 hold on
else
plot(M22(:,j))
hold on
end

J111(:,j)=M11(51:61,j);
 J222(:,j)=M22(51:61,j);

for k=1:100
F111(:,j,k)=poissrnd(pars1(1,4)*J111(:,j));%prob put from Spread sheet 
F222(:,j,k)=poissrnd(0.72*J222(:,j));%take from optimizeITN Co.m
end  


   DataMean2_2(:,j)=mean(F111(:,j,k),3);

   eta1=0.2;
eta2=0.1;
eta3=0.2;
eta4=0.01;
 Q2(1,j)=trapz(eta1*DataMean2_2(:,j)+eta2*abs(B(j,1)-AA1)+eta3*abs(B(j,2)-AA2)+eta4*abs(B(j,3)-0.5));


end

[val2, idx2] = min(Q2);

Malcas1=352.5;
 EMalcas2020=Malcas1-(Malcas1*(.40))
 EMalcas2025 = Malcas1-(Malcas1*(.75))
 EMalcas2030 = Malcas1-(Malcas1*(.90))

 W1(1,:)=J111(end,:);
 ind1 = interp1(W1,1:length(W1),EMalcas2030,'nearest','extrap');   % find index position in (cuppa, r) less than and closet to Expected malaria incidence cases 
 Q2_1=Q2(1,idx2)
 Q2_2=Q2(1,ind1)
 Q2_3=Q2(1,11)
 AAA2=B(idx2,:)
 BBB2=B(ind1,:)
 CCC2=B(11,:)

figure(222)
 subplot(121)
 plot(M11(:,idx2))
 subplot(122)
 plot(M22(:,idx2))
% P1_1=B(idx2,:);
  J1111(:,1)=M11(51:61,idx2);
   J2222(:,1)=M22(51:61,idx2);

for k=1:100
F1111(:,k)=poissrnd(pars1(1,4)*J1111);
F2222(:,k)=poissrnd(0.72*J2222);
end  
% % % 
 H=cumsum(F1111,1);
 H1=cumsum(J1111,1);
 MeanInc1= H1(end);
 H2=H(end,:);
 
figure(987)
histogram(H2,10)
hold on 
xline(MeanInc1);
% % % 
 X11=1:11;
 X22=1:11;
 Data11=F1111;
 Data22=F2222;
N11=size(Data11,1);
N22=size(Data22,1);
 DataMean11=mean(Data11,2);
 DataMean22=mean(Data22,2);
 DataSEM11=std(Data11,0,2)/sqrt(N11);
 DataSEM22=std(Data22,0,2)/sqrt(N22);
 CI9511=tinv([0.025 0.975],N11-1);
 CI9522=tinv([0.025 0.975],N22-1);

  DataCI9511=bsxfun(@times,DataSEM11,CI9511.*ones(11,1));
  DataCI9522=bsxfun(@times,DataSEM22,CI9522.*ones(11,1));
 CIdist11=abs(DataCI9511-DataMean11);
 CIdist22=abs(DataCI9522-DataMean22);
% 
Y11=CIdist11(:,1);
Z11=CIdist11(:,2);
Y22=CIdist22(:,1);
Z22=CIdist22(:,2);

% 
% check whether incidence corresp to (kappa, r) that 
% minimizes the cost function is less than GTS 2030 
OptMalcas2030=DataMean11(end)
EMalcas2030

F11111(:,1)=poissrnd(pars1(1,4)*J1111);
F22222(:,1)=poissrnd(0.72*J2222);

 
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
hold on
 scatter(1:2:length(Togo_extractInci(:,2)), all_pointsInci, 'k', 'filled', 'Marker', 'o');
hold on
 xline(41,'LineWidth',1.2)
 X111 = [45:51, fliplr(45:51)];
inBetween = [Y1', fliplr(Z1')];
 
 fill(X111, inBetween, 'g');
 hold on
 plot(45:51,DataMean1)
 hold on
 xline(51)
 X1111 = [51:61, fliplr(51:61)];
 inBetween = [Y11', fliplr(Z11')];
 fill(X1111, inBetween, 'r');
 hold on 
 plot(51:61,DataMean11)

hold on
P1=31;
P2=Malcas1;
plot(P1,P2,'r*','LineWidth',3,'MarkerSize',10)
text(P1,P2,['(', num2str(P1), ', ', num2str(P2), ')'])
P3= EMalcas2020;
P4=41;
text(P4,P3,['(', num2str(P4), ', ', num2str(P3), ')'])
hold on
plot(P4,P3,'b*','LineWidth',3,'MarkerSize',10)
text(P4,P3,['(', num2str(P4), ', ', num2str(P3), ')'])
hold on
scatter(P4,AchieveMalcas2020, 'k', 'filled', 'Marker', 'o')
hold on
scatter(43, AchieveMalcas2021, 'k', 'filled', 'Marker', 'o')
hold on
scatter(45, AchieveMalcas2022, 'k', 'filled', 'Marker', 'o')
hold on
P5=51;
P6= EMalcas2025;
plot(P5,P6,'b*','LineWidth',3,'MarkerSize',10)
text(P5,P6,['(', num2str(P5), ', ', num2str(P6), ')'])
hold on
% plot(P5,p0,'b*','LineWidth',3,'MarkerSize',10)

P7=61;
P8= EMalcas2030;
plot(P7,P8,'b*','LineWidth',3,'MarkerSize',10)
text(P7,P8,['(', num2str(P7), ', ', num2str(P8), ')'])

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

scatter(41,  AchieveITNUse2020 ,'r', 'filled', 'Marker', 'o')
hold on 
scatter(43,  AchieveITNUse2021 ,'r', 'filled', 'Marker', 'o')
hold on 
scatter(45,  AchieveITNUse2022 ,'r', 'filled', 'Marker', 'o')
hold on
xline(41)
 hold on
 X222 = [45:51, fliplr(45:51)];
inBetween = [Y2', fliplr(Z2')];
fill(X222, inBetween, 'g');
hold on
plot(45:51, DataMean2)
hold on
xline(51)
X2222 = [51:61, fliplr(51:61)];
 inBetween = [Y22', fliplr(Z22')];
 fill(X2222, inBetween, 'r');
 hold on 
 plot(51:61,DataMean22)

figure(333) 
subplot(121)
 plot(1:45,Y11_1(1:45))
 hold on
 scatter(1:2:length(Togo_extractInci(:,2)), all_pointsInci, 'k', 'filled', 'Marker', 'o','LineWidth',8);
 hold on
scatter(41,AchieveMalcas2020, 'k', 'filled', 'Marker', 'o','LineWidth',8)
 hold on
scatter(43,AchieveMalcas2021, 'k', 'filled', 'Marker', 'o','LineWidth',8)
 hold on
scatter(45,AchieveMalcas2022, 'k', 'filled', 'Marker', 'o','LineWidth',8)
hold on
subplot(122)
plot(1:45,Y22_2(1:45))
  hold on
scatter(1:2:length(Togo_extractITN(:,2)), all_pointsITN, 'red', 'filled', 'Marker', 'o','LineWidth',8);
hold on
scatter(41,AchieveITNUse2020, 'r', 'filled', 'Marker', 'o','LineWidth',8)
hold on
scatter(43,AchieveITNUse2021, 'r', 'filled', 'Marker', 'o','LineWidth',8)
hold on
scatter(45,AchieveITNUse2022, 'r', 'filled', 'Marker', 'o','LineWidth',8)

