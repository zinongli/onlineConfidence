s = maxSpeed;
error = copy(:,29);
f = @(theta,s,e) -sum(-log(theta(3).*s + theta(4)) - (((error - (theta(1).*s + theta(2))).^2)./ (2.*(theta(3).*s+theta(4)))));
fun = @(theta) f(theta,s,error);
theta0 = [0,0,1,10];
UB = [1,30,1,30];
LB = [-1,-30,0,1];
theta = bads(fun,theta0,LB,UB);

%%
load('data_onlineConf/JX/JX_MaxSpeed.mat','maxSpeed')
load('JX_rightMaxSpeed.mat','rightward')
load('JX_LeftMaxSpeed.mat','leftward')
xTotal = copy(:,23);
xRight = copy(copy(:,8)<1000,23);
xLeft = copy(copy(:,8)>1000,23);

sampleRatio = 0.5;
test_n = 50;
theta0 = [0,0,1,5];
UB = [1,2,1,10];
LB = [-1,-2,0,1];
f = @(theta,s,error) -sum(-log(theta(3).*s + theta(4)) - (((error - (theta(1).*s + theta(2))).^2)./ (2.*(theta(3).*s+theta(4)))));


xRightLLH = NaN(1,test_n);
xLeftLLH = NaN(1,test_n);
xLateralLLH = NaN(1,test_n);
xTotalLLH = NaN(1,test_n);




for i = 1:test_n
    
    tic
    
    xTotalSampleBin = zeros(1,length(xTotal));
    xTotalSampleBin(randsample(length(xTotal),round(length(xTotal)*sampleRatio))) = 1;
    xTotalSampleBin = logical(xTotalSampleBin);
    xTotalSample = xTotal(xTotalSampleBin);
    xTotalTest = xTotal(~xTotalSampleBin);
    maxSpeedSample = maxSpeed(xTotalSampleBin);
    maxSpeedTest = maxSpeed(~xTotalSampleBin);

    xRightSampleBin = zeros(1,length(xRight));
    xRightSampleBin(randsample(length(xRight),round(length(xRight)*sampleRatio))) = 1;
    xRightSampleBin = logical(xRightSampleBin);
    xRightSample = xRight(xRightSampleBin);
    xRightTest = xRight(~xRightSampleBin);
    rightwardSample = rightward(xRightSampleBin);
    rightwardTest = rightward(~xRightSampleBin);

    xLeftSampleBin = zeros(1,length(xLeft));
    xLeftSampleBin(randsample(length(xLeft),round(length(xLeft)*sampleRatio))) = 1;
    xLeftSampleBin = logical(xLeftSampleBin);
    xLeftSample = xLeft(xLeftSampleBin);
    xLeftTest = xLeft(~xLeftSampleBin);
    leftwardSample = leftward(xLeftSampleBin);
    leftwardTest = leftward(~xLeftSampleBin);

    xParams = NaN(3,4);
    fun = @(theta) f(theta,maxSpeedSample,xTotalSample);
    [xParams(1,:),~,~,~] = bads(fun,theta0,LB,UB);
    fun = @(theta) f(theta,rightwardSample,xRightSample);
    [xParams(2,:),~,~,~] = bads(fun,theta0,LB,UB);
    fun = @(theta) f(theta,leftwardSample,xLeftSample);
    [xParams(3,:),~,~,~] = bads(fun,theta0,LB,UB);
    
    xTotalLLH(i) = -sum(-log(xParams(1,3).*maxSpeedTest + xParams(1,4)) - (((xTotalTest - (xParams(1,1).*maxSpeedTest + xParams(1,2))).^2)./ (2.*(xParams(1,3).*maxSpeedTest+xParams(1,4)))));
    xRightLLH(i) = -sum(-log(xParams(2,3).*rightwardTest + xParams(2,4)) - (((xRightTest - (xParams(2,1).*rightwardTest + xParams(2,2))).^2)./ (2.*(xParams(2,3).*rightwardTest+xParams(2,4)))));
    xLeftLLH(i) = -sum(-log(xParams(3,3).*leftwardTest + xParams(3,4)) - (((xLeftTest - (xParams(3,1).*leftwardTest + xParams(3,2))).^2)./ (2.*(xParams(3,3).*leftwardTest+xParams(3,4)))));
    xLateralLLH(i) = xRightLLH(i) + xLeftLLH(i);
    
    toc
    i
end

xModel = NaN(test_n,2);
xModel(:,1) = xTotalLLH;
xModel(:,2) = xLateralLLH;

yRightLLH = NaN(1,test_n);
yLeftLLH = NaN(1,test_n);
yLateralLLH = NaN(1,test_n);
yTotalLLH = NaN(1,test_n);




for i = 1:test_n
    
    tic
    
    yTotalSampleBin = zeros(1,length(yTotal));
    yTotalSampleBin(randsample(length(yTotal),round(length(yTotal)*sampleRatio))) = 1;
    yTotalSampleBin = logical(yTotalSampleBin);
    yTotalSample = yTotal(yTotalSampleBin);
    yTotalTest = yTotal(~yTotalSampleBin);
    maxSpeedSample = maxSpeed(yTotalSampleBin);
    maxSpeedTest = maxSpeed(~yTotalSampleBin);

    yRightSampleBin = zeros(1,length(yRight));
    yRightSampleBin(randsample(length(yRight),round(length(yRight)*sampleRatio))) = 1;
    yRightSampleBin = logical(yRightSampleBin);
    yRightSample = yRight(yRightSampleBin);
    yRightTest = yRight(~yRightSampleBin);
    rightwardSample = rightward(yRightSampleBin);
    rightwardTest = rightward(~yRightSampleBin);

    yLeftSampleBin = zeros(1,length(yLeft));
    yLeftSampleBin(randsample(length(yLeft),round(length(yLeft)*sampleRatio))) = 1;
    yLeftSampleBin = logical(yLeftSampleBin);
    yLeftSample = yLeft(yLeftSampleBin);
    yLeftTest = yLeft(~yLeftSampleBin);
    leftwardSample = leftward(yLeftSampleBin);
    leftwardTest = leftward(~yLeftSampleBin);

    yParams = NaN(3,4);
    fun = @(theta) f(theta,maxSpeedSample,yTotalSample);
    [yParams(1,:),~,~,~] = bads(fun,theta0,LB,UB);
    fun = @(theta) f(theta,rightwardSample,yRightSample);
    [yParams(2,:),~,~,~] = bads(fun,theta0,LB,UB);
    fun = @(theta) f(theta,leftwardSample,yLeftSample);
    [yParams(3,:),~,~,~] = bads(fun,theta0,LB,UB);
    
    yTotalLLH(i) = -sum(-log(yParams(1,3).*maxSpeedTest + yParams(1,4)) - (((yTotalTest - (yParams(1,1).*maxSpeedTest + yParams(1,2))).^2)./ (2.*(yParams(1,3).*maxSpeedTest+yParams(1,4)))));
    yRightLLH(i) = -sum(-log(yParams(2,3).*rightwardTest + yParams(2,4)) - (((yRightTest - (yParams(2,1).*rightwardTest + yParams(2,2))).^2)./ (2.*(yParams(2,3).*rightwardTest+yParams(2,4)))));
    yLeftLLH(i) = -sum(-log(yParams(3,3).*leftwardTest + yParams(3,4)) - (((yLeftTest - (yParams(3,1).*leftwardTest + yParams(3,2))).^2)./ (2.*(yParams(3,3).*leftwardTest+yParams(3,4)))));
    yLateralLLH(i) = yRightLLH(i) + yLeftLLH(i);
    
    toc
    i+50
end

yModel = NaN(test_n,2);
yModel(:,1) = yTotalLLH;
yModel(:,2) = yLateralLLH;


%% Model Comparison
set(groot,'defaultAxesFontSize',18)
figure(1)
subplot(1,2,1)
xModel_n = size(xModel,2);
xSemMSE = NaN(1,xModel_n);
xMeanMSE = NaN(1,xModel_n);
for i = 1:xModel_n
    xSemMSE(i) = std(xModel(:,i));
    xMeanMSE(i) = -mean(xModel(:,i));
end
errorbar(1:2,xMeanMSE,xSemMSE,'o','MarkerSize',18,'LineWidth',2)
xlim([0,3])
xlabel('Model')
ylabel('Log Likelihood')
xticks(0:3)
xticklabels(["", "Total", "Individal", ""])
title('Model Comparison for Along Error')
grid on
subplot(1,2,2)
yModel_n = size(yModel,2);
ySemMSE = NaN(1,yModel_n);
yMeanMSE = NaN(1,yModel_n);
for i = 1:yModel_n
    ySemMSE(i) = std(yModel(:,i));
    yMeanMSE(i) = -mean(yModel(:,i));
end
errorbar(1:2,yMeanMSE,ySemMSE,'o','MarkerSize',18,'LineWidth',2)
xlim([0,3])
xlabel('Model')
ylabel('Log Likelihood')
xticks(0:3)
xticklabels(["", "Total", "Individal", ""])
title('Model Comparison for Orthogonal Error')
grid on
sgtitle('Test N = 50, Errorbar for STD','FontSize',18)
%%
load('data_onlineConf/JX/JX_MaxSpeed.mat','maxSpeed')
load('JX_rightMaxSpeed.mat','rightward')
load('JX_LeftMaxSpeed.mat','leftward')
xTotal = copy(:,23);
xRight = copy(copy(:,8)<1000,23);
xLeft = copy(copy(:,8)>1000,23);


theta0 = [0,0,0,5];
UB = [1,25,1,50];
LB = [-1,-1,0,1];

f = @(theta,speed,error) log(1/sqrt(2*pi)) * length(speed) - sum(-log(theta(3).*speed + theta(4)) - (((error - (theta(1).*speed + theta(2))).^2)./ (2.*(theta(3).*speed+theta(4)).^2)));


xParams = NaN(3,4);
fun = @(theta) f(theta,maxSpeed,xTotal);
xParams(1,:) = bads(fun,theta0,LB,UB);
fun = @(theta) f(theta,rightward,xRight);
xParams(2,:) = bads(fun,theta0,LB,UB);
fun = @(theta) f(theta,leftward,xLeft);
xParams(3,:) = bads(fun,theta0,LB,UB);

%%
yTotal = copy(:,29);
yRight = copy(copy(:,8)<1000,29);
yLeft = copy(copy(:,8)>1000,29);


theta0 = [0,0,0,5];
UB = [1,25,1,50];
LB = [-1,-1,0,1];

f = @(theta,speed,error) log(1/sqrt(2*pi)) * length(speed) - sum(-log(theta(3).*speed + theta(4)) - (((error - (theta(1).*speed + theta(2))).^2)./ (2.*(theta(3).*speed+theta(4)).^2)));


yParams = NaN(3,4);
fun = @(theta) f(theta,maxSpeed,yTotal);
yParams(1,:) = bads(fun,theta0,LB,UB);
fun = @(theta) f(theta,rightward,yRight);
yParams(2,:) = bads(fun,theta0,LB,UB);
fun = @(theta) f(theta,leftward,yLeft);
yParams(3,:) = bads(fun,theta0,LB,UB);
%%
xTotalLLH = log(1/sqrt(2*pi)) * length(maxSpeed) + sum(-log(xParams(1,3).*maxSpeed + xParams(1,4)) - (((xTotal - (xParams(1,1).*maxSpeed + xParams(1,2))).^2)./ (2.*(xParams(1,3).*maxSpeed+xParams(1,4)))));
xRightLLH = log(1/sqrt(2*pi)) * length(rightward) + sum(-log(xParams(2,3).*rightward + xParams(2,4)) - (((xRight - (xParams(2,1).*rightward + xParams(2,2))).^2)./ (2.*(xParams(2,3).*rightward+xParams(2,4)))));
xLeftLLH = log(1/sqrt(2*pi)) * length(leftward) + sum(-log(xParams(3,3).*leftward + xParams(3,4)) - (((xLeft - (xParams(3,1).*leftward + xParams(3,2))).^2)./ (2.*(xParams(3,3).*leftward+xParams(3,4)))));
xLateralLLH = xRightLLH + xLeftLLH;

aicTotal = 2 * 4 - 2 * xTotalLLH
aicIndividual = 2 * 8 - 2 * xLateralLLH
bicTotal = 4 * log(length(maxSpeed)) - 2 * xTotalLLH
bicIndividual = 8 * log(length(maxSpeed)) - 2 * xLateralLLH
%%
figure(1)
indices = 1:1600;
subplot(1,3,1)

plot(maxSpeed,xTotal,'o')
hold on
plot(indices,(xParams(1,1) .* indices) + xParams(1,2),'--r')
plot(indices,(xParams(1,1) .* indices) + xParams(1,2) + (xParams(1,3) * indices + xParams(1,4)),'--r')
plot(indices,(xParams(1,1) .* indices) + xParams(1,2) - (xParams(1,3) * indices + xParams(1,4)),'--r')
plot(indices,(xParams(1,1) .* indices) + xParams(1,2) + 2*(xParams(1,3) * indices + xParams(1,4)),'--r')
plot(indices,(xParams(1,1) .* indices) + xParams(1,2) - 2*(xParams(1,3) * indices + xParams(1,4)),'--r')
hold off
xlabel('Speed mm/s')
ylabel('Error mm')

subplot(1,3,2)
plot(rightward,xRight,'o')
hold on
plot(indices,(xParams(2,1) .* indices) + xParams(2,2),'--r')
plot(indices,(xParams(2,1) .* indices) + xParams(2,2) + (xParams(2,3) * indices + xParams(2,4)),'--r')
plot(indices,(xParams(2,1) .* indices) + xParams(2,2) - (xParams(2,3) * indices + xParams(2,4)),'--r')
plot(indices,(xParams(2,1) .* indices) + xParams(2,2) + 2*(xParams(2,3) * indices + xParams(2,4)),'--r')
plot(indices,(xParams(2,1) .* indices) + xParams(2,2) - 2*(xParams(2,3) * indices + xParams(2,4)),'--r')
hold off
xlabel('Speed mm/s')
ylabel('Error mm')

subplot(1,3,3)
plot(leftward,xLeft,'o')
hold on
plot(indices,(xParams(3,1) .* indices) + xParams(3,2),'--r')
plot(indices,(xParams(3,1) .* indices) + xParams(3,2) + (xParams(3,3) * indices + xParams(3,4)),'--r')
plot(indices,(xParams(3,1) .* indices) + xParams(3,2) - (xParams(3,3) * indices + xParams(3,4)),'--r')
plot(indices,(xParams(3,1) .* indices) + xParams(3,2) + 2*(xParams(3,3) * indices + xParams(3,4)),'--r')
plot(indices,(xParams(3,1) .* indices) + xParams(3,2) - 2*(xParams(3,3) * indices + xParams(3,4)),'--r')
hold off
xlabel('Speed mm/s')
ylabel('Error mm')

sgtitle('Error Along the Reach Direction','FontSize',18)

%%
figure(2)
indices = 1:1600;
subplot(1,3,1)

plot(maxSpeed,yTotal,'o')
hold on
plot(indices,(yParams(1,1) .* indices) + yParams(1,2),'--r')
plot(indices,(yParams(1,1) .* indices) + yParams(1,2) + (yParams(1,3) * indices + yParams(1,4)),'--r')
plot(indices,(yParams(1,1) .* indices) + yParams(1,2) - (yParams(1,3) * indices + yParams(1,4)),'--r')
plot(indices,(yParams(1,1) .* indices) + yParams(1,2) + 2*(yParams(1,3) * indices + yParams(1,4)),'--r')
plot(indices,(yParams(1,1) .* indices) + yParams(1,2) - 2*(yParams(1,3) * indices + yParams(1,4)),'--r')
hold off
xlabel('Speed mm/s')
ylabel('Error mm')

subplot(1,3,2)
plot(rightward,yRight,'o')
hold on
plot(indices,(yParams(2,1) .* indices) + yParams(2,2),'--r')
plot(indices,(yParams(2,1) .* indices) + yParams(2,2) + (yParams(2,3) * indices + yParams(2,4)),'--r')
plot(indices,(yParams(2,1) .* indices) + yParams(2,2) - (yParams(2,3) * indices + yParams(2,4)),'--r')
plot(indices,(yParams(2,1) .* indices) + yParams(2,2) + 2*(yParams(2,3) * indices + yParams(2,4)),'--r')
plot(indices,(yParams(2,1) .* indices) + yParams(2,2) - 2*(yParams(2,3) * indices + yParams(2,4)),'--r')
hold off
xlabel('Speed mm/s')
ylabel('Error mm')

subplot(1,3,3)
plot(leftward,yLeft,'o')
hold on
plot(indices,(yParams(3,1) .* indices) + yParams(3,2),'--r')
plot(indices,(yParams(3,1) .* indices) + yParams(3,2) + (yParams(3,3) * indices + yParams(3,4)),'--r')
plot(indices,(yParams(3,1) .* indices) + yParams(3,2) - (yParams(3,3) * indices + yParams(3,4)),'--r')
plot(indices,(yParams(3,1) .* indices) + yParams(3,2) + 2*(yParams(3,3) * indices + yParams(3,4)),'--r')
plot(indices,(yParams(3,1) .* indices) + yParams(3,2) - 2*(yParams(3,3) * indices + yParams(3,4)),'--r')
hold off
xlabel('Speed mm/s')
ylabel('Error mm')

sgtitle('Error Orthogonal the Reach Direction','FontSize',18)