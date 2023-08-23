% copy column contents:
% 1,2: target x and y in wac pixels 
% 3: time out 
% 4,5: the onset and end time of the reach
% 6,7: endpoint x and y in wac pixels
% 8,9: start position in wac pixels
% 10: target distance in mm
% 11,12: target x and y in mm
% 13,14: endpoint x and y in mm
% 15: target size in mm
% 16: actual duration of the reach
% 17: error size in mm
% 18: switch logical array
% 19,20: start position in mm
% 21: actual reach distance
% 22: average speed
% 23: error along the reach direction (vector projection) in mm
% 24,25: relative target position in mm
% 26: score
% 27: trial order
% 28: hit rate
% 29: error orthogonal to the reach direction (vector rejection) in mm
%%
subj = 'JX';
load('data_onlineConf/JX/JM_practice_traj_S1_06-Apr-2023_tform.mat')
load('data_onlineConf/JX/JM_practice_traj_S1c_06-Apr-2023_rawtotal.mat')
load('data_onlineConf/JX/JM_practice_traj_S1c_06-Apr-2023_traXtotal.mat')
load('data_onlineConf/JX/JM_practice_traj_S1c_06-Apr-2023_traYtotal.mat')
load('data_onlineConf/JX/JX_MaxSpeed.mat')
load('data_onlineConf/JX/JX_rightMaxSpeed.mat')
load('data_onlineConf/JX/JX_leftMaxSpeed.mat')
load('data_onlineConf/JX/JX_SAparams8.mat')
load('data_onlineConf/JX/JX_OptSpeed8.mat')


index = NaN(length(data),1);
for i = 1:length(data)
    index(i) = ~isnan(sum(data(i,:)));
end
valid = data(index==true,:);
validTraX = traXtotal(index==true,:);
validTraY = traYtotal(index==true,:);

Affine2d =tform.T(1:2,1:2);
[~,s,~] = svd(Affine2d);
proj2tablet = 1./mean([s(1,1),s(2,2)]);
pixellength = 0.248;
copy = valid;
copy(:,[1,2]) = transformPointsInverse(tform,copy(:,[1,2]));
copy(:,[8,9]) = transformPointsInverse(tform,copy(:,[8,9]));
copy(:,10) = sqrt(sum((copy(:,1:2) - copy(:,8:9)).^2,2)) .* pixellength;
copy(:,[11,12]) = [copy(:,1)*pixellength (1080 - copy(:,2))*pixellength];
copy(:,[13,14]) = [copy(:,6)*pixellength (1080 - copy(:,7))*pixellength]; % 1080 = tablet pixel height
copy(:,15) = valid(:,10) .* pixellength .* proj2tablet ./ 2; % proj2tablet = projetor size to tablet size (physical size), /2 is diameter vs radius
copy(:,16) = copy(:,5) - copy(:,4);
copy(:,17) = sqrt( (copy(:,13)-copy(:,11)).^2 + (copy(:,14)-copy(:,12)).^2 );
copy(:,27) = 1:length(copy);
copy(:,19:20) = (copy(:,6:7) - copy(:,8:9)) .* pixellength;
copy(:,21) = sqrt(sum((copy(:,6:7) - copy(:,8:9)).^2,2)) .* pixellength;
copy(:,22) = copy(:,21) ./ copy(:,16);
copy(:,24:25) = (copy(:,1:2) - copy(:,8:9)) .* pixellength;% relative target coordinate
copy(:,23) = (abs(dot(copy(:,19:20),copy(:,24:25),2) ./ dot(copy(:,24:25),copy(:,24:25),2)) - 1) .*copy(:,10);
copy(:,26) = valid(:,11);
copy(:,28) = copy(:,26) ~= 0;
endPoints = (copy(:,6:7) - copy(:,8:9)) .* pixellength;% relative endpoint coordinate
projScale = (abs(dot(copy(:,19:20),copy(:,24:25),2) ./ dot(copy(:,24:25),copy(:,24:25),2)));
rejections = endPoints - projScale.* copy(:,24:25);
tooRight = NaN(1,length(copy));
for i = 1:length(copy)
    tooRight(i) = sum(cross([endPoints(i,:),0],[rejections(i,:),0])<0); 
end
tooRight(tooRight==0) = -1;
rejLength = sqrt(rejections(:,1).^2 + rejections(:,2).^2);
copy(:,29) = rejLength .* tooRight';

%% sigmoid fitting max speed from trajs
reCenteredTrajX = NaN(360,270);
reCenteredTrajY = NaN(360,270);

for i = 1:360
    reCenteredTrajX(i,:) = (validTraX(i,:) - copy(i,8)) .* pixellength;
    reCenteredTrajY(i,:) = (validTraY(i,:) - copy(i,9)) .* pixellength;
end

maxDuration = sum((sum(~isnan(reCenteredTrajX),1))~=0);
traProjection = NaN(360,maxDuration);
for i = 1:maxDuration
    traProjection(:,i) = (abs(dot([reCenteredTrajX(:,i),reCenteredTrajY(:,i)],copy(:,24:25),2) ./ dot(copy(:,24:25),copy(:,24:25),2)));
end
speedProjection = traProjection(:,2:end) - traProjection(:,1:end-1);
accProjection = speedProjection(:,2:end) - speedProjection(:,1:end-1);

endTime = sum(~isnan(reCenteredTrajX),2);
b0 = [1,-1,45];
sigmoidFit = NaN(3,360);
for i = 1:360
    x = 50:endTime(i);
    y = traProjection(i,x);
    x = x - 50;
    fun = @(b) (b(1)./(1+exp(b(2)*(x-b(3)))))-y;
    b = lsqnonlin(fun,b0);
    sigmoidFit(:,i) = b;
end
maxSpeed = -sigmoidFit(2,:)'.*0.25.*copy(:,10).*60;
%% testing with animation
% % feel free to comment this section as it's only for error checking
% for i = 1:360
%     x = 50:endTime(i);
%     y = traProjection(i,x);
%     x = x - 50;
%     b = sigmoidFit(:,i);
%     adjustedX = x - b(3);
%     plot(adjustedX,b(1)./(1+exp(b(2)*(x-b(3)))),'-')
%     hold on
%     plot(adjustedX,y,'o')
%     xline(0);
%     yline(1);
%     yline(0);
%     hold off
%     xlim([-45 45])
%     ylim([-0.1,1.2])
%     pause(0.1)
% end
%% Hit Rate 
% just wanna make sure the participants aren't hitting celings or floor all
% the time
tSizes = unique(copy(:,15));
SizeDistHitRate = NaN(5,3);
for i = 1:5
    for j = 1:3
        SizeDistHitRate(i,j) = mean(copy(copy(:,15)==tSizes(i) & round(copy(:,10),-2) == j*100,28));
    end
end
figure
for i = 1:3
    plot(tSizes,SizeDistHitRate(:,i),'-o')
    hold on
end

plot(tSizes,0.3:0.1:0.7,'r--')
hold off
xlabel('Target Size (mm)')
ylabel('Hit Rate')
legend('~100 mm','~200 mm','~300 mm','Expectation','Location','Northwest')
title('Hit Rates of Each Distance Range')

%% linear fit max & avg as param_0
plot(maxSpeed,copy(:,22),'o')

xlabel('Max Speed mm/s')
ylabel('Average Speed mm/s')
title('Max v.s. Average Speed, each point = one trial')
mld = fitlm(maxSpeed,copy(:,22));
%% MLE for max & avg speed relationship 
% this section is redundant as it will be done again later
% but it and the next section's graphs serve to validate the choice of using 
% MLE and that with a changing sigma
theta0 = [table2array(mld.Coefficients(2,1)),table2array(mld.Coefficients(1,1)),0.0001,30];
UB = [5,200,5,300];
LB = [0,-200,0.0001,1];

f = @(theta,speed,error) -(log(1/sqrt(2*pi)) * length(speed) + sum(-log(theta(3).*speed + theta(4)) - (((error - (theta(1).*speed + theta(2))).^2)./ (2.*(theta(3).*speed+theta(4)).^2))));
fun = @(theta) f(theta,maxSpeed,copy(:,22));
changingSig = bads(fun,theta0,LB,UB);

theta0 = [table2array(mld.Coefficients(2,1)),table2array(mld.Coefficients(1,1)),30];
UB = [5,200,300];
LB = [0,-200,1];
f2 = @(theta,speed,error) -(log(1/sqrt(2*pi)) * length(speed) + sum(-log(theta(3)) - (((error - (theta(1).*speed + theta(2))).^2)./ (2.*(theta(3)).^2))));
fun2 = @(theta) f2(theta,maxSpeed,copy(:,22));
fixSig = bads(fun2,theta0,LB,UB);

changingaic = 2 * 4 - 2 * -fun(changingSig);
changingbic = 4 * log(length(maxSpeed)) - 2 * -fun(changingSig);
fixaic = 2 * 5 - 2 * -fun2(fixSig);
fixbic = 5 * log(length(maxSpeed)) - 2 * -fun2(fixSig);
%% graphs for MLE max & avg 
indices = 1:1600;
plot(maxSpeed,copy(:,22),'o')
hold on
plot(indices,(changingSig(1) .* indices) + changingSig(2),'--r')
plot(indices,(changingSig(1) .* indices) + changingSig(2) + 1.5*(changingSig(3) * indices + changingSig(4)),'--r')
plot(indices,(changingSig(1) .* indices) + changingSig(2) - 1.5*(changingSig(3) * indices + changingSig(4)),'--r')
% plot(indices,(changingSig(1) .* indices) + changingSig(2) + 2*(changingSig(3) * indices + changingSig(4)),'--r')
% plot(indices,(changingSig(1) .* indices) + changingSig(2) - 2*(changingSig(3) * indices + changingSig(4)),'--r')
hold off
xlabel('Speed mm/s')
ylabel('Error mm')

%% Bias, Error, and Avg Speed fit in respect to Max Speed
set(groot,'defaultAxesFontSize',18)
load('data_onlineConf/JX/JX_MaxSpeed.mat','maxSpeed')

xTotal = copy(:,23);
xRight = copy(copy(:,8)<1000,23);
xLeft = copy(copy(:,8)>1000,23);

theta0 = [0,0,0.0001,5];
UB = [1,30,1,50];
LB = [-1,-30,0.0001,1];

f = @(theta,speed,error) -(log(1/sqrt(2*pi)) * length(speed) + sum(-log(theta(3).*speed + theta(4)) - (((error - (theta(1).*speed + theta(2))).^2)./ (2.*(theta(3).*speed+theta(4)).^2))));

xParams = NaN(3,4);
fun = @(theta) f(theta,maxSpeed,xTotal);
xParams(1,:) = bads(fun,theta0,LB,UB);
fun = @(theta) f(theta,maxSpeed(copy(:,8)<1000),xRight);
xParams(2,:) = bads(fun,theta0,LB,UB);
fun = @(theta) f(theta,maxSpeed(copy(:,8)>1000),xLeft);
xParams(3,:) = bads(fun,theta0,LB,UB);

yTotal = copy(:,29);
yRight = copy(copy(:,8)<1000,29);
yLeft = copy(copy(:,8)>1000,29);

theta0 = [0,0,0.0001,5];
UB = [1,30,1,50];
LB = [-1,-30,0.0001,1];

f = @(theta,speed,error) -(log(1/sqrt(2*pi)) * length(speed) + sum(-log(theta(3).*speed + theta(4)) - (((error - (theta(1).*speed + theta(2))).^2)./ (2.*(theta(3).*speed+theta(4)).^2))));

yParams = NaN(3,4);
fun = @(theta) f(theta,maxSpeed,yTotal);
yParams(1,:) = bads(fun,theta0,LB,UB);
fun = @(theta) f(theta,maxSpeed(copy(:,8)<1000),yRight);
yParams(2,:) = bads(fun,theta0,LB,UB);
fun = @(theta) f(theta,maxSpeed(copy(:,8)>1000),yLeft);
yParams(3,:) = bads(fun,theta0,LB,UB);

theta0 = [0,0,0.0001,30];
UB = [5,100,5,100];
LB = [0,-100,0.0001,1];

f = @(theta,speed,error) -(log(1/sqrt(2*pi)) * length(speed) + sum(-log(theta(3).*speed + theta(4)) - (((error - (theta(1).*speed + theta(2))).^2)./ (2.*(theta(3).*speed+theta(4)).^2))));
fun = @(theta) f(theta,maxSpeed,copy(:,22));
changingSig = bads(fun,theta0,LB,UB);

RightParams = [xParams(2,1:2);yParams(2,1:2);xParams(2,3:4);yParams(2,3:4);changingSig(1:2);changingSig(3:4)];
LeftParams = [xParams(3,1:2);yParams(3,1:2);xParams(3,3:4);yParams(3,3:4);changingSig(1:2);changingSig(3:4)];

cache_R = RightParams';
cache_L = LeftParams';
SAparams = [cache_R,cache_L];

save('data_onlineConf/JX/JX_SAparams.mat','SAparams');
