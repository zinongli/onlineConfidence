%% Parameter Generation
% Because this is a cross validation on the linear fit of BIAS, the null
% hypothesis is zero for both parameters: slope and intercept
set(groot,'defaultAxesFontSize',18)
load('data_onlineConf/JX/JX_MaxSpeed.mat','maxSpeed')
load('JX_rightMaxSpeed.mat','rightward')
load('JX_LeftMaxSpeed.mat','leftward')
xTotal = copy(:,23);
xRight = copy(copy(:,8)<1000,23);
xLeft = copy(copy(:,8)>1000,23);

warning('off')
sampleRatio = 0.5;
test_n = 100;


xRightMSE = NaN(1,test_n);
xLeftMSE = NaN(1,test_n);
x1LateralMSE = NaN(1,test_n);
x2LateralMSE = NaN(1,test_n);
x3LateralMSE = NaN(1,test_n);
x4TotalMSE = NaN(1,test_n);
x5LateralMSE = NaN(1,test_n);
x6LateralMSE = NaN(1,test_n);
x7LateralMSE = NaN(1,test_n);
x8LateralMSE = NaN(1,test_n);
x9LateralMSE = NaN(1,test_n);
for i = 1:test_n

    xTotalSampleBin = zeros(1,length(xTotal));
    xTotalSampleBin(randsample(length(xTotal),round(length(xTotal)*sampleRatio))) = 1;
    xTotalSampleBin = logical(xTotalSampleBin);
    xTotalSample = xTotal(xTotalSampleBin);
    xTotalTest = xTotal(~xTotalSampleBin);
    maxSpeedSample = maxSpeed(xTotalSampleBin);
    maxSpeedTest = maxSpeed(~xTotalSampleBin);
    save('xTotalSampleSpeed,mat','maxSpeedSample')

    xRightSampleBin = zeros(1,length(xRight));
    xRightSampleBin(randsample(length(xRight),round(length(xRight)*sampleRatio))) = 1;
    xRightSampleBin = logical(xRightSampleBin);
    xRightSample = xRight(xRightSampleBin);
    xRightTest = xRight(~xRightSampleBin);
    rightwardSample = rightward(xRightSampleBin);
    rightwardTest = rightward(~xRightSampleBin);
    save('xRightSampleSpeed,mat','rightwardSample')

    xLeftSampleBin = zeros(1,length(xLeft));
    xLeftSampleBin(randsample(length(xLeft),round(length(xLeft)*sampleRatio))) = 1;
    xLeftSampleBin = logical(xLeftSampleBin);
    xLeftSample = xLeft(xLeftSampleBin);
    xLeftTest = xLeft(~xLeftSampleBin);
    leftwardSample = leftward(xLeftSampleBin);
    leftwardTest = leftward(~xLeftSampleBin);
    save('xLeftSampleSpeed,mat','leftwardSample')

    xParams = NaN(3,2);
    [phat,~] = mle(xTotalSample,'pdf',@(xTotalSample,a,b,c,d) mlPDFfittingTotal(xTotalSample,a,b,c,d),'Start',[0,0,0,5]);
    xParams(1,:) = phat(1:2);
    [phat,~] = mle(xRightSample,'pdf',@(xRightSample,a,b,c,d) mlPDFfittingRight(xRightSample,a,b,c,d),'Start',[0,0,0,5]);
    xParams(2,:) = phat(1:2);
    [phat,~] = mle(xLeftSample,'pdf',@(xLeftSample,a,b,c,d) mlPDFfittingLeft(xLeftSample,a,b,c,d),'Start',[0,0,0,5]);
    xParams(3,:) = phat(1:2);

    x4TotalMSE(i) = sum((xTotalTest - (xParams(1,1) .* maxSpeedTest + xParams(1,2))).^2) ./ length(maxSpeedTest);
    xRightMSE(i) = sum((xRightTest - (xParams(2,1) .* rightwardTest + xParams(2,2))).^2) ./ length(rightwardTest);
    xLeftMSE(i) = sum((xLeftTest - (xParams(3,1) .* leftwardTest + xParams(3,2))).^2) ./ length(leftwardTest);
    x1LateralMSE(i) = (sum((xRightTest - (xParams(2,1) .* rightwardTest + xParams(2,2))).^2) + sum((xLeftTest - (xParams(3,1) .* leftwardTest + xParams(3,2))).^2))./ (length(rightwardTest)+length(leftwardTest));
    x2LateralMSE(i) = (sum((xRightTest - (xParams(2,1) .* rightwardTest + xParams(1,2))).^2) + sum((xLeftTest - (xParams(3,1) .* leftwardTest + xParams(1,2))).^2))./ (length(rightwardTest)+length(leftwardTest));
    x3LateralMSE(i) = (sum((xRightTest - (xParams(1,1) .* rightwardTest + xParams(2,2))).^2) + sum((xLeftTest - (xParams(1,1) .* leftwardTest + xParams(3,2))).^2))./ (length(rightwardTest)+length(leftwardTest));
    x5LateralMSE(i) = (sum((xRightTest - (xParams(2,1) .* rightwardTest + 0)).^2) + sum((xLeftTest - (xParams(3,1) .* leftwardTest + 0)).^2))./ (length(rightwardTest)+length(leftwardTest));
    x6LateralMSE(i) = (sum((xRightTest - (0 .* rightwardTest + xParams(2,2))).^2) + sum((xLeftTest - (0 .* leftwardTest + xParams(3,2))).^2))./ (length(rightwardTest)+length(leftwardTest));
    x7LateralMSE(i) = (sum((xRightTest - (xParams(1,1) .* rightwardTest + 0)).^2) + sum((xLeftTest - (xParams(1,1) .* leftwardTest + 0)).^2))./ (length(rightwardTest)+length(leftwardTest));
    x8LateralMSE(i) = (sum((xRightTest - (0 .* rightwardTest + xParams(1,2))).^2) + sum((xLeftTest - (0 .* leftwardTest + xParams(1,2))).^2))./ (length(rightwardTest)+length(leftwardTest));
    x9LateralMSE(i) = (sum((xRightTest - (0 .* rightwardTest + 0)).^2) + sum((xLeftTest - (0 .* leftwardTest + 0)).^2))./ (length(rightwardTest)+length(leftwardTest));
    
    i
end

warning('on')
xModel = NaN(test_n,9);
xModel(:,1) = x1LateralMSE;
xModel(:,2) = x2LateralMSE;
xModel(:,3) = x3LateralMSE;
xModel(:,4) = x4TotalMSE;
xModel(:,5) = x5LateralMSE;
xModel(:,6) = x6LateralMSE;
xModel(:,7) = x7LateralMSE;
xModel(:,8) = x8LateralMSE;
xModel(:,9) = x9LateralMSE;


yTotal = copy(:,29);
yRight = copy(copy(:,8)<1000,29);
yLeft = copy(copy(:,8)>1000,29);

warning('off')

yRightMSE = NaN(1,test_n);
yLeftMSE = NaN(1,test_n);
y1LateralMSE = NaN(1,test_n);
y2LateralMSE = NaN(1,test_n);
y3LateralMSE = NaN(1,test_n);
y4TotalMSE = NaN(1,test_n);
y5LateralMSE = NaN(1,test_n);
y6LateralMSE = NaN(1,test_n);
y7LateralMSE = NaN(1,test_n);
y8LateralMSE = NaN(1,test_n);
y9LateralMSE = NaN(1,test_n);
fprintf('Orthogonal MSE Progress:    ');
for i = 1:test_n
    yTotalSampleBin = zeros(1,length(yTotal));
    yTotalSampleBin(randsample(length(yTotal),round(length(yTotal)*sampleRatio))) = 1;
    yTotalSampleBin = logical(yTotalSampleBin);
    yTotalSample = yTotal(yTotalSampleBin);
    yTotalTest = yTotal(~yTotalSampleBin);
    maxSpeedSample = maxSpeed(yTotalSampleBin);
    maxSpeedTest = maxSpeed(~yTotalSampleBin);
    save('yTotalSampleSpeed,mat','maxSpeedSample')

    yRightSampleBin = zeros(1,length(yRight));
    yRightSampleBin(randsample(length(yRight),round(length(yRight)*sampleRatio))) = 1;
    yRightSampleBin = logical(yRightSampleBin);
    yRightSample = yRight(yRightSampleBin);
    yRightTest = yRight(~yRightSampleBin);
    rightwardSample = rightward(yRightSampleBin);
    rightwardTest = rightward(~yRightSampleBin);
    save('yRightSampleSpeed,mat','rightwardSample')

    yLeftSampleBin = zeros(1,length(yLeft));
    yLeftSampleBin(randsample(length(yLeft),round(length(yLeft)*sampleRatio))) = 1;
    yLeftSampleBin = logical(yLeftSampleBin);
    yLeftSample = yLeft(yLeftSampleBin);
    yLeftTest = yLeft(~yLeftSampleBin);
    leftwardSample = leftward(yLeftSampleBin);
    leftwardTest = leftward(~yLeftSampleBin);
    save('yLeftSampleSpeed,mat','leftwardSample')

    yParams = NaN(3,2);
    [phat,~] = mle(yTotalSample,'pdf',@(yTotalSample,a,b,c,d) mlPDFfittingTotal(yTotalSample,a,b,c,d),'Start',[0,0,0,5]);
    yParams(1,:) = phat(1:2);
    [phat,~] = mle(yRightSample,'pdf',@(yRightSample,a,b,c,d) mlPDFfittingRight(yRightSample,a,b,c,d),'Start',[0,0,0,5]);
    yParams(2,:) = phat(1:2);
    [phat,~] = mle(yLeftSample,'pdf',@(yLeftSample,a,b,c,d) mlPDFfittingLeft(yLeftSample,a,b,c,d),'Start',[0,0,0,5]);
    yParams(3,:) = phat(1:2);

    y4TotalMSE(i) = sum((yTotalTest - (yParams(1,1) .* maxSpeedTest + yParams(1,2))).^2) ./ length(maxSpeedTest);
    yRightMSE(i) = sum((yRightTest - (yParams(2,1) .* rightwardTest + yParams(2,2))).^2) ./ length(rightwardTest);
    yLeftMSE(i) = sum((yLeftTest - (yParams(3,1) .* leftwardTest + yParams(3,2))).^2) ./ length(leftwardTest);
    y1LateralMSE(i) = (sum((yRightTest - (yParams(2,1) .* rightwardTest + yParams(2,2))).^2) + sum((yLeftTest - (yParams(3,1) .* leftwardTest + yParams(3,2))).^2))./ (length(rightwardTest)+length(leftwardTest));
    y2LateralMSE(i) = (sum((yRightTest - (yParams(2,1) .* rightwardTest + yParams(1,2))).^2) + sum((yLeftTest - (yParams(3,1) .* leftwardTest + yParams(1,2))).^2))./ (length(rightwardTest)+length(leftwardTest));
    y3LateralMSE(i) = (sum((yRightTest - (yParams(1,1) .* rightwardTest + yParams(2,2))).^2) + sum((yLeftTest - (yParams(1,1) .* leftwardTest + yParams(3,2))).^2))./ (length(rightwardTest)+length(leftwardTest));
    y5LateralMSE(i) = (sum((yRightTest - (yParams(2,1) .* rightwardTest + 0)).^2) + sum((yLeftTest - (yParams(3,1) .* leftwardTest + 0)).^2))./ (length(rightwardTest)+length(leftwardTest));
    y6LateralMSE(i) = (sum((yRightTest - (0 .* rightwardTest + yParams(2,2))).^2) + sum((yLeftTest - (0 .* leftwardTest + yParams(3,2))).^2))./ (length(rightwardTest)+length(leftwardTest));
    y7LateralMSE(i) = (sum((yRightTest - (yParams(1,1) .* rightwardTest + 0)).^2) + sum((yLeftTest - (yParams(1,1) .* leftwardTest + 0)).^2))./ (length(rightwardTest)+length(leftwardTest));
    y8LateralMSE(i) = (sum((yRightTest - (0 .* rightwardTest + yParams(1,2))).^2) + sum((yLeftTest - (0 .* leftwardTest + yParams(1,2))).^2))./ (length(rightwardTest)+length(leftwardTest));
    y9LateralMSE(i) = (sum((yRightTest - (0 .* rightwardTest + 0)).^2) + sum((yLeftTest - (0 .* leftwardTest + 0)).^2))./ (length(rightwardTest)+length(leftwardTest));

    perc = floor(100*i/test_n);
    fprintf('\b\b\b\b%4d%%',perc);
end

warning('on')
yModel = NaN(test_n,9);
yModel(:,1) = y1LateralMSE;
yModel(:,2) = y2LateralMSE;
yModel(:,3) = y3LateralMSE;
yModel(:,4) = y4TotalMSE;
yModel(:,5) = y5LateralMSE;
yModel(:,6) = y6LateralMSE;
yModel(:,7) = y7LateralMSE;
yModel(:,8) = y8LateralMSE;
yModel(:,9) = y9LateralMSE;
%% Model Comparison
figure(1)
subplot(1,2,1)
xModel_n = size(xModel,2);
xSemMSE = NaN(1,xModel_n);
xMeanMSE = NaN(1,xModel_n);
for i = 1:model_n
    xSemMSE(i) = std(xModel(:,i)) ./ test_n;
    xMeanMSE(i) = mean(xModel(:,i));
end
errorbar(1:9,xMeanMSE,xSemMSE,'o','MarkerSize',8,'LineWidth',2)
xlim([0,10])
xlabel('Model #')
ylabel('SEM of MSE')
xticks(0:10)
title('Model Comparison for Along Bias')

subplot(1,2,2)
yModel_n = size(yModel,2);
ySemMSE = NaN(1,yModel_n);
yMeanMSE = NaN(1,yModel_n);
for i = 1:model_n
    ySemMSE(i) = std(yModel(:,i)) ./ test_n;
    yMeanMSE(i) = mean(yModel(:,i));
end
subplot(1,2,1)
errorbar(1:9,yMeanMSE,ySemMSE,'o','MarkerSize',8,'LineWidth',2)
xlim([0,10])
xlabel('Model #')
ylabel('SEM of MSE')
xticks(0:10)
title('Model Comparison for Orthogonal Bias')
%%
% | Model | Intercept | Slope |
% |   1   |   indiv   | indiv |
% |   2   |   total   | indiv |
% |   3   |   indiv   | total |
% |   4   |   total   | total |
% |   5   |     0     | indiv |
% |   6   |   indiv   |   0   |
% |   7   |     0     | total |
% |   8   |   total   |   0   |
% |   9   |     0     |   0   |

