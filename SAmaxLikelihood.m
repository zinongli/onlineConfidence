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

xRightSample = xRight;
rightwardSample = rightward;
save('xRightSampleSpeed.mat','rightwardSample')

xLeftSample = xLeft;
leftwardSample = leftward;
save('xLeftSampleSpeed.mat','leftwardSample')

xParams = NaN(4,2);

[phat,~] = mle(xRightSample,'pdf',@(xRightSample,a,b,c,d) xMLEfittingRight(xRightSample,a,b,c,d),'Start',[0,0,0,5]);
xParams(:,1) = phat(1:4);
[phat,~] = mle(xLeftSample,'pdf',@(xLeftSample,a,b,c,d) xMLEfittingLeft(xLeftSample,a,b,c,d),'Start',[0,0,0,5]);
xParams(:,2) = phat(1:4);

warning('on')


yTotal = copy(:,29);
yRight = copy(copy(:,8)<1000,29);
yLeft = copy(copy(:,8)>1000,29);

warning('off')

yRightSample = yRight;
rightwardSample = rightward;
save('yRightSampleSpeed.mat','rightwardSample')
yLeftSample = yLeft;
leftwardSample = leftward;
save('yLeftSampleSpeed.mat','leftwardSample')

yParams = NaN(4,2);

[phat,~] = mle(yRightSample,'pdf',@(yRightSample,a,b,c,d) yMLEfittingRight(yRightSample,a,b,c,d),'Start',[0,0,0,5]);
yParams(:,1) = phat(1:4);
[phat,~] = mle(yLeftSample,'pdf',@(yLeftSample,a,b,c,d) yMLEfittingLeft(yLeftSample,a,b,c,d),'Start',[0,0,0,5]);
yParams(:,2) = phat(1:4);

warning('on')

%%
