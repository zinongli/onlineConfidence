%% Parameter Estimation
xTotal = copy(:,23);
xRight = copy(copy(:,8)<1000,23);
xLeft = copy(copy(:,8)>1000,23);
xParams = NaN(3,2);
[phat,~] = mle(xTotal,'pdf',@(xTotal,a,b,c,d) mlPDFfittingTotal(xTotal,a,b,c,d),'Start',[0,0,0,5]);
xParams(1,:) = phat(1:2);
[phat,~] = mle(xRight,'pdf',@(xRight,a,b,c,d) mlPDFfittingRight(xRight,a,b,c,d),'Start',[0,0,0,5]);
xParams(2,:) = phat(1:2);
[phat,~] = mle(xLeft,'pdf',@(xLeft,a,b,c,d) mlPDFfittingLeft(xLeft,a,b,c,d),'Start',[0,0,0,5]);
xParams(3,:) = phat(1:2);

%% Parameter Generation
sampleRatio = 0.5;
yTotal = copy(:,29);
yRight = copy(copy(:,8)<1000,29);
yLeft = copy(copy(:,8)>1000,29);
load('data_onlineConf/JX/JX_MaxSpeed.mat','maxSpeed')
load('JX_rightMaxSpeed.mat','rightward')
load('JX_LeftMaxSpeed.mat','leftward')

yTotalSampleBin = zeros(1,length(yTotal));
yTotalSampleBin(randsample(length(yTotal),round(length(yTotal)*sampleRatio))) = 1;
yTotalSampleBin = logical(yTotalSampleBin);
yTotalSample = yTotal(yTotalSampleBin);
yTotalTest = yTotal(~yTotalSampleBin);
maxSpeed = maxSpeed(yTotalSampleBin);
save('yTotalSampleSpeed,mat','maxSpeed')

yRightSampleBin = zeros(1,length(yRight));
yRightSampleBin(randsample(length(yRight),round(length(yRight)*sampleRatio))) = 1;
yRightSampleBin = logical(yRightSampleBin);
yRightSample = yRight(yRightSampleBin);
yRightTest = yRight(~yRightSampleBin);
rightward = rightward(yRightSampleBin);
save('yRightSampleSpeed,mat','rightward')

yLeftSampleBin = zeros(1,length(yLeft));
yLeftSampleBin(randsample(length(yLeft),round(length(yLeft)*sampleRatio))) = 1;
yLeftSampleBin = logical(yLeftSampleBin);
yLeftSample = yLeft(yLeftSampleBin);
yLeftTest = yLeft(~yLeftSampleBin);
leftward = leftward(yLeftSampleBin);
save('yLeftSampleSpeed,mat','leftward')

yParams = NaN(3,2);
[phat,~] = mle(yTotalSample,'pdf',@(yTotalSample,a,b,c,d) mlPDFfittingTotal(yTotalSample,a,b,c,d),'Start',[0,0,0,5]);
yParams(1,:) = phat(1:2);
[phat,~] = mle(yRightSample,'pdf',@(yRightSample,a,b,c,d) mlPDFfittingRight(yRightSample,a,b,c,d),'Start',[0,0,0,5]);
yParams(2,:) = phat(1:2);
[phat,~] = mle(yLeftSample,'pdf',@(yLeftSample,a,b,c,d) mlPDFfittingLeft(yLeftSample,a,b,c,d),'Start',[0,0,0,5]);
yParams(3,:) = phat(1:2);

%%

