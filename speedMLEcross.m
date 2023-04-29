set(groot,'defaultAxesFontSize',18)
load('data_onlineConf/JX/JX_MaxSpeed.mat','maxSpeed')
load('JX_rightMaxSpeed.mat','rightward')
load('JX_LeftMaxSpeed.mat','leftward')
sTotal = copy(:,22);
sRight = copy(copy(:,8)<1000,22);
sLeft = copy(copy(:,8)>1000,22);

warning('off')
sampleRatio = 0.5;
test_n = 100;


sRightMSE = NaN(1,test_n);
sLeftMSE = NaN(1,test_n);
s1LateralMSE = NaN(1,test_n);
s2LateralMSE = NaN(1,test_n);
s3LateralMSE = NaN(1,test_n);
s4TotalMSE = NaN(1,test_n);
s5LateralMSE = NaN(1,test_n);
s6LateralMSE = NaN(1,test_n);
s7LateralMSE = NaN(1,test_n);
s8LateralMSE = NaN(1,test_n);
s9LateralMSE = NaN(1,test_n);
for i = 1:test_n

    sTotalSampleBin = zeros(1,length(sTotal));
    sTotalSampleBin(randsample(length(sTotal),round(length(sTotal)*sampleRatio))) = 1;
    sTotalSampleBin = logical(sTotalSampleBin);
    sTotalSample = sTotal(sTotalSampleBin);
    sTotalTest = sTotal(~sTotalSampleBin);
    maxSpeedSample = maxSpeed(sTotalSampleBin);
    maxSpeedTest = maxSpeed(~sTotalSampleBin);

    sRightSampleBin = zeros(1,length(sRight));
    sRightSampleBin(randsample(length(sRight),round(length(sRight)*sampleRatio))) = 1;
    sRightSampleBin = logical(sRightSampleBin);
    sRightSample = sRight(sRightSampleBin);
    sRightTest = sRight(~sRightSampleBin);
    rightwardSample = rightward(sRightSampleBin);
    rightwardTest = rightward(~sRightSampleBin);

    sLeftSampleBin = zeros(1,length(sLeft));
    sLeftSampleBin(randsample(length(sLeft),round(length(sLeft)*sampleRatio))) = 1;
    sLeftSampleBin = logical(sLeftSampleBin);
    sLeftSample = sLeft(sLeftSampleBin);
    sLeftTest = sLeft(~sLeftSampleBin);
    leftwardSample = leftward(sLeftSampleBin);
    leftwardTest = leftward(~sLeftSampleBin);

    sParams = NaN(3,2);
    
    mld = fitlm(maxSpeedSample,sTotalSample);
    sParams(1,:) = table2array(mld.Coefficients(:,1));
    mld = fitlm(rightwardSample,sRightSample);
    sParams(2,:) = table2array(mld.Coefficients(:,1));
    mld = fitlm(leftwardSample,sLeftSample);
    sParams(3,:) = table2array(mld.Coefficients(:,1));

    s4TotalMSE(i) = sum((sTotalTest - (sParams(1,1) .* maxSpeedTest + sParams(1,2))).^2) ./ length(maxSpeedTest);
    sRightMSE(i) = sum((sRightTest - (sParams(2,1) .* rightwardTest + sParams(2,2))).^2) ./ length(rightwardTest);
    sLeftMSE(i) = sum((sLeftTest - (sParams(3,1) .* leftwardTest + sParams(3,2))).^2) ./ length(leftwardTest);
    s1LateralMSE(i) = (sum((sRightTest - (sParams(2,1) .* rightwardTest + sParams(2,2))).^2) + sum((sLeftTest - (sParams(3,1) .* leftwardTest + sParams(3,2))).^2))./ (length(rightwardTest)+length(leftwardTest));
    s2LateralMSE(i) = (sum((sRightTest - (sParams(2,1) .* rightwardTest + sParams(1,2))).^2) + sum((sLeftTest - (sParams(3,1) .* leftwardTest + sParams(1,2))).^2))./ (length(rightwardTest)+length(leftwardTest));
    s3LateralMSE(i) = (sum((sRightTest - (sParams(1,1) .* rightwardTest + sParams(2,2))).^2) + sum((sLeftTest - (sParams(1,1) .* leftwardTest + sParams(3,2))).^2))./ (length(rightwardTest)+length(leftwardTest));
    s5LateralMSE(i) = (sum((sRightTest - (sParams(2,1) .* rightwardTest + 0)).^2) + sum((sLeftTest - (sParams(3,1) .* leftwardTest + 0)).^2))./ (length(rightwardTest)+length(leftwardTest));
    s6LateralMSE(i) = (sum((sRightTest - (0 .* rightwardTest + sParams(2,2))).^2) + sum((sLeftTest - (0 .* leftwardTest + sParams(3,2))).^2))./ (length(rightwardTest)+length(leftwardTest));
    s7LateralMSE(i) = (sum((sRightTest - (sParams(1,1) .* rightwardTest + 0)).^2) + sum((sLeftTest - (sParams(1,1) .* leftwardTest + 0)).^2))./ (length(rightwardTest)+length(leftwardTest));
    s8LateralMSE(i) = (sum((sRightTest - (0 .* rightwardTest + sParams(1,2))).^2) + sum((sLeftTest - (0 .* leftwardTest + sParams(1,2))).^2))./ (length(rightwardTest)+length(leftwardTest));
    s9LateralMSE(i) = (sum((sRightTest - (0 .* rightwardTest + 0)).^2) + sum((sLeftTest - (0 .* leftwardTest + 0)).^2))./ (length(rightwardTest)+length(leftwardTest));
    
    i
end

warning('on')
sModel = NaN(test_n,9);
sModel(:,1) = s1LateralMSE;
sModel(:,2) = s2LateralMSE;
sModel(:,3) = s3LateralMSE;
sModel(:,4) = s4TotalMSE;
sModel(:,5) = s5LateralMSE;
sModel(:,6) = s6LateralMSE;
sModel(:,7) = s7LateralMSE;
sModel(:,8) = s8LateralMSE;
sModel(:,9) = s9LateralMSE;
%%
figure(1)

sModel_n = size(sModel,2);
sSemMSE = NaN(1,sModel_n);
sMeanMSE = NaN(1,sModel_n);
for i = 1:model_n
    sSemMSE(i) = std(sModel(:,i)) ./ test_n;
    sMeanMSE(i) = mean(sModel(:,i));
end
errorbar(1:9,sMeanMSE,sSemMSE,'o','MarkerSize',8,'LineWidth',2)
xlim([0,10])
xlabel('Model #')
ylabel('SEM of MSE')
xticks(0:10)
title('Model Comparison for Avg/Max Speed Fit')
%%
% | Model | Intercept | Slope |
% |-------|-----------|-------|
% |   1   |   indiv   | indiv |
% |   2   |   total   | indiv |
% |   3   |   indiv   | total |
% |   4   |   total   | total |
% |   5   |     0     | indiv |
% |   6   |   indiv   |   0   |
% |   7   |     0     | total |
% |   8   |   total   |   0   |
% |   9   |     0     |   0   |