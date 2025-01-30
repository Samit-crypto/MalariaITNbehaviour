%% Evaluate estimated parameters (beta_0, kappa_0, r)

load Togo_extractInci.mat;
  Incid1=Togo_extractInci(:,2);
  load Togo_extractITN.mat;
 ITN1= Togo_extractITN(:,2);

AchieveMalcas2020=236.2;
AchieveMalcas2021=231.8;
AchieveMalcas2022=230.5;
AchieveMalcas2020_2021=(AchieveMalcas2020+AchieveMalcas2021)/2;
AchieveMalcas2021_2022=(AchieveMalcas2021+AchieveMalcas2022)/2;

AchieveITNUse2020= 80.8;  % manipulated
AchieveITNUse2021=76.9;
AchieveITNUse2022= 64.3 ;
AchieveITNUse2020_2021=(AchieveITNUse2020+AchieveITNUse2021)/2;
AchieveITNUse2021_2022=(AchieveITNUse2021+AchieveITNUse2022)/2;


Incid1=[Incid1;AchieveMalcas2020;AchieveMalcas2020_2021;...
        AchieveMalcas2021;AchieveMalcas2021_2022;AchieveMalcas2022];
ITN1=[ITN1;AchieveITNUse2020;AchieveITNUse2020_2021;...
    AchieveITNUse2021;AchieveITNUse2021_2022;AchieveITNUse2022];


beta0= 215.6; %beta
r=0.1;      %relative risk
cuppa=0.55;  % imitation

[TP, ITNP, XX, pars, AA]=optimizemastITNTogo2(beta0, r, cuppa, Incid1, ITN1);

 save('TP',TP,'ITNP1',ITNP,'XX1',XX,'pars1',pars,'AA1',AA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Perform sensitivity analysis using Bootstrap method

Y11_1=TP(161:211);
 Y22_2=ITNP(161:211);
  NewY11_1=Y11_1(1:2:end);
  NewY22_2=Y22_2(1:2:end);

modelInci = TP(161:205,1);
modelITN = ITNP(161:205,1);
ErrInci =Incid1-modelInci;
ErrITN = ITN1-modelITN;
numArrays = 20;

% Initialize a cell array to store the generated arrays
generatedArraysInci = cell(1, numArrays);
generatedArraysITN = cell(1, numArrays);


% Generate the arrays using permutations with replacement
for i = 1:numArrays
%     perturbationsInci = lowerBound + (upperBound - lowerBound) * rand(1, numel(ErrInci));
%     generatedArraysInci{i}=ErrInci+perturbationsInci;
     generatedArraysInci{i} = datasample(ErrInci, numel(ErrInci), 'Replace', true);

% perturbationsITN = lowerBound + (upperBound - lowerBound) * rand(1, numel(ErrITN));
%     generatedArraysITN{i}=ErrITN+perturbationsITN;
    generatedArraysITN{i} = datasample(ErrITN, numel(ErrITN), 'Replace', true);
end

SyndataInci=generatedArraysInci{1,1};
SyndataITN=generatedArraysITN{1,1};

for i = 2:numArrays
SyndataInci = horzcat(SyndataInci, generatedArraysInci{1,i});
SyndataITN = horzcat(SyndataITN, generatedArraysITN{1,i});
end


SyndataInciNew= modelInci+SyndataInci;
SyndataITNNew= modelITN+SyndataITN;


TPdata=cell(1,numArrays);
ITNPdata=cell(1,numArrays);
XXdata=cell(1,numArrays);
parsdata=cell(1,numArrays);
AAdata=cell(1,numArrays);

for i=1:numArrays
i
 try
Incid2=SyndataInciNew(:,i);
ITN2=SyndataITNNew(:,i);

beta0= 215 + ( 217-215) * rand; %beta
r=0.05 + (0.15 - 0.05) * rand;     %relative risk
cuppa=0.54 + (0.56 - 0.54) * rand;  % imitation

[TP, ITNP, XX, pars, AA]=optimizemastITNTogo2(beta0, r, cuppa, Incid2, ITN2);

 catch
            % Handle the error (e.g., display a warning or handle exceptional cases)
            warning('An error occurred during fmincon optimization. Continuing with next iteration.');
  end
TPdata{1, i} = TP;
ITNPdata{1, i} =ITNP;
XXdata{1,i}=XX;
parsdata{1,i}=pars;
AAdata{1,i}=AA;

end

Estparams=parsdata{1,1};
MalInc=TPdata{1,1};
ITNuse=ITNPdata{1,1};

for i=2:20
 Estparams = vertcat(Estparams, parsdata{1,i});
MalInc = horzcat(MalInc, TPdata{1,i});
 ITNuse=horzcat(ITNuse,ITNPdata{1,i});
end
Estparams(:,4)=[];

% Assuming your matrix is stored in the variable 'data'
confidenceLevel = 0.95; % Specify the desired confidence level

% Get the number of columns in the matrix
numColumns = size(Estparams, 2);

% Initialize an array to store the confidence intervals
confidenceIntervals = zeros(2, numColumns);

% Iterate over each column
for col = 1:numColumns
    % Extract the data from the current column
    columnData = Estparams(:, col);
    
    % Compute the mean and standard deviation of the column data
    meanValue = mean(columnData);
    stdValue = std(columnData);
    
    % Compute the sample size
    sampleSize = numel(columnData);
    
    % Compute the standard error
    standardError = stdValue / sqrt(sampleSize);
    
    % Compute the critical value (t-score) based on the confidence level
    df = sampleSize - 1; % Degrees of freedom
    alpha = 1 - confidenceLevel;
    criticalValue = tinv(1 - alpha / 2, df);
    
    % Compute the lower and upper bounds of the confidence interval for the column
    lowerBound = meanValue - criticalValue * standardError;
    upperBound = meanValue + criticalValue * standardError;
    
    % Store the confidence interval in the array
    confidenceIntervals(:, col) = [lowerBound; upperBound];
end

% Display the confidence intervals for each column
for col = 1:numColumns
    fprintf('Column %d Confidence Interval: %.2f - %.2f\n', col, confidenceIntervals(1, col), confidenceIntervals(2, col));
end


% Assuming your matrix is stored in the variable 'data'
startRow = 161; % Starting row index
endRow = 205; % Ending row index
numColumns = 20; % Number of columns to plot

% Extract the desired subset of the matrix
subsetInci = MalInc(startRow:endRow, 1:numColumns);
subsetITN = ITNuse(startRow:endRow, 1:numColumns);


% Create a figure and axis for the plot
figure(121)
hold on

% Iterate over each column to plot
for col = 1:numColumns
    % Get the data for the current column
    columnDataInci(:,col) = subsetInci(:, col);
    columnDataITN(:,col) = subsetITN(:, col);
    
    % Plot the data
   subplot(121) 
   plot(startRow:endRow, columnDataInci,'Color', [0.65 0.65 0.65])
   hold on
   plot(startRow:endRow, TP(startRow:endRow),'k')
   subplot(122)
    plot(startRow:endRow, columnDataITN,'Color', [0.65 0.65 0.65])
    hold on
    plot(startRow:endRow, ITNP(startRow:endRow),'k')
end
