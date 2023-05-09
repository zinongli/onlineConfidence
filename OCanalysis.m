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

if subj == 'JX'
    load(['data_onlineConf/' subj '/' 'JM_practice_traj_S1_06-Apr-2023_tform.mat'])
    load(['data_onlineConf/' subj '/' 'JM_practice_traj_S1c_06-Apr-2023_rawtotal.mat'])
    load(['data_onlineConf/' subj '/' 'JM_practice_traj_S1c_06-Apr-2023_traXtotal.mat'])
    load(['data_onlineConf/' subj '/' 'JM_practice_traj_S1c_06-Apr-2023_traYtotal.mat'])
    load(['data_onlineConf/' subj '/' subj '_MaxSpeed.mat'])
    load(['data_onlineConf/' subj '/' subj '_rightMaxSpeed.mat'])
    load(['data_onlineConf/' subj '/' subj '_leftMaxSpeed.mat'])
    load(['data_onlineConf/' subj '/' subj '_SAparams8.mat'])
    load(['data_onlineConf/' subj '/' subj '_OptSpeed8.mat'])
elseif subj == 'ZL'
    load(['data_onlineConf/' subj '/' 'pilot_practice_traj_S1_31-Mar-2023_tform.mat'])
    load(['data_onlineConf/' subj '/' 'rawComplete.mat'])
    load(['data_onlineConf/' subj '/' 'xTrajComplete.mat'])
    load(['data_onlineConf/' subj '/' 'yTrajComplete.mat'])
    load(['data_onlineConf/' subj '/' subj '_MaxSpeed.mat'])
    load(['data_onlineConf/' subj '/' subj '_rightMaxSpeed.mat'])
    load(['data_onlineConf/' subj '/' subj '_leftMaxSpeed.mat'])
    load(['data_onlineConf/' subj '/' subj '_SAparams8.mat'])
    load(['data_onlineConf/' subj '/' subj '_OptSpeed8.mat'])

end

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
copy([1:60 181:240],3) = 3;
copy([61:120 241:300],3) = 2;
copy([121:180 301:360],3) = 1;
%%
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
%%
for i = 1:360
    x = 50:endTime(i);
    y = traProjection(i,x);
    x = x - 50;
    b = sigmoidFit(:,i);
    adjustedX = x - b(3);
    plot(adjustedX,b(1)./(1+exp(b(2)*(x-b(3)))),'-')
    hold on
    plot(adjustedX,y,'o')
    xline(0);
    yline(1);
    yline(0);
    hold off
    xlim([-45 45])
    ylim([-0.1,1.2])
    pause(0.1)
end


%% grouping by range (e.g. 50-100, 100-150, etc)
[sortspd,ind] = sort(copy(:,22));
group_n = 15;
stds = NaN(1,group_n);
hitrate = NaN(1,group_n);
means = NaN(1,group_n);
speedStd = NaN(1,group_n);
speedPoint = NaN(1,group_n);
for i = 1:group_n
    stds(i) = std(copy((i*50 < copy(:,22) & i*50+50 > copy(:,22)),23));
    hitrate(i) = sum(copy((i*50 < copy(:,22) & i*50+50 > copy(:,22)),28));
    means(i) = mean(copy((i*50 < copy(:,22) & i*50+50 > copy(:,22)),23));
    speedStd(i) = std(copy(i*50 < copy(:,22) & i*50+50 > copy(:,22),22));
    speedPoint(i) = mean(copy(i*50 < copy(:,22) & i*50+50 > copy(:,22),22));
    sum(i*50 < copy(:,22) & i*50+50 > copy(:,22))
end

mld = fitlm(speedPoint,stds)

errorstdX = std(copy(:,23));
bootMax = 1000;
numTestTrials = 50;

numConds = 13;
numCondd = 14;
bootXs = NaN(numTestTrials,numConds);
bootXd = NaN(numTestTrials,numCondd);
stdBoots = NaN(bootMax,numConds);
meanBoots = NaN(bootMax,numConds);

for ii = 1:bootMax
    
    for jj = 1:numConds
        i = jj+2;
        selectionLogic = i*50 < copy(:,22) & i*50+50 > copy(:,22);
        numTrials = sum(selectionLogic);
        bootInd = randi(numTrials,numTestTrials,1);
        subgroup = copy(selectionLogic,23);
        bootXs(:,jj) = subgroup(bootInd); 
%         bootY(:,jj) = epMatY(bootInd(:,jj),jj);
    end

    bootstdXs = std(bootXs);
    bootmeanXs = mean(bootXs);
%     bootstdY = std(bootY);
%     bootstdXY = sqrt(bootstdX.^2 + bootstdY.^2);
    stdBoots(ii,:) = bootstdXs;
    meanBoots(ii,:) = bootmeanXs;
end

errorsErrors = NaN(1,numConds);
for i = 1:numConds
    errorsErrors(i) = std(stdBoots(:,i));
end

figure
errorbar(speedPoint(3:15),stds(3:15),errorsErrors,'vertical','-bo')
hold on 
errorbar(speedPoint(3:15),stds(3:15),speedStd(3:15),'horizontal','-bo')
hold off
ylim([0,25])
ylabel('STD along the reach (mm)')
xlabel('Average speed of the reach (mm/s)')
title('Errorbar data from 1000 bootstrap samples of 50 trials each')
%% Hit Rate
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

%% grouping by order (e.g. first 30, second 30, etc.)
[~,ind] = sort(copy(:,22));
SpeedSortedTrials = copy(ind,:);
group_n = 20;
stds = NaN(1,group_n);
hitrate = NaN(1,group_n);
means = NaN(1,group_n);
speedStd = NaN(1,group_n);
speedPoint = NaN(1,group_n);
groupSize = length(copy)/group_n;
for i = 1:group_n
    stds(i) = std(SpeedSortedTrials(i*groupSize-groupSize+1:i*groupSize,29));
    speedPoint(i) = mean(SpeedSortedTrials(i*groupSize-groupSize+1:i*groupSize,22));
    hitrate(i) = sum(SpeedSortedTrials(i*groupSize-groupSize+1:i*groupSize,28));
    means(i) = mean(SpeedSortedTrials(i*groupSize-groupSize+1:i*groupSize,29));
    speedStd(i) = std(SpeedSortedTrials(i*groupSize-groupSize+1:i*groupSize,22));
end
figure
errorbar(speedPoint,stds,speedStd,'horizontal','o')
xlabel('Average Speed (mm/s)')
ylabel('Standard Error (mm)')
title('Average Speed vs Error')
mld = fitlm(speedPoint,stds)
%% max speed order, fixed group size, along
[sortedMaxSpeed,ind] = sort(maxSpeed);
maxSpeedSortedTrials = copy(ind,:);
group_ns = [6,8,9,10,12,15,18,20,24,30,36,40];
interceptFits = NaN(1,length(group_ns));
slopeFits = NaN(1,length(group_ns));
rSquareds = NaN(1,length(group_ns));
LLH = NaN(1,length(group_ns));
for j = 1:length(group_ns)
    
    group_n = group_ns(j);
    stds = NaN(1,group_n);
    hitrate = NaN(1,group_n);
    means = NaN(1,group_n);
    speedStd = NaN(1,group_n);
    speedPoint = NaN(1,group_n);
    groupSize = length(copy)/group_n;
    for i = 1:group_n
        stds(i) = std(maxSpeedSortedTrials(i*groupSize-groupSize+1:i*groupSize,23));
        speedPoint(i) = mean(sortedMaxSpeed(i*groupSize-groupSize+1:i*groupSize));
        hitrate(i) = sum(maxSpeedSortedTrials(i*groupSize-groupSize+1:i*groupSize,28));
        means(i) = mean(maxSpeedSortedTrials(i*groupSize-groupSize+1:i*groupSize,23));
        speedStd(i) = std(sortedMaxSpeed(i*groupSize-groupSize+1:i*groupSize));
    end
    mld = fitlm(speedPoint,stds);
    interceptFits(j) = table2array(mld.Coefficients(1,1));
    slopeFits(j) = table2array(mld.Coefficients(2,1));
    rSquareds(j) = mld.Rsquared.Adjusted;
    LLH(j) = mld.LogLikelihood;
end

figure
plot(slopeFits,'-o')
ylabel('Slope')
xlabel('Group Size')
xticks(1:length(group_ns))
xticklabels(group_ns)
title('Slope vs Group Size')
figure
plot(interceptFits,'-o')
ylabel('Intercept')
xlabel('Group Size')
xticks(1:length(group_ns))
xticklabels(group_ns)
title('Intercept vs Group Size')
figure
plot(rSquareds,'-o')
ylabel('R Sqaured')
xlabel('Group Size')
xticks(1:length(group_ns))
xticklabels(group_ns)
title('R Sqaured vs Group Size')

% figure
% errorbar(speedPoint,stds,speedStd,'horizontal','o')
% xlabel('Max Speed (mm/s)')
% ylabel('Standard Error (mm)')
% title('Max Speed vs Error')
% mld = fitlm(speedPoint,stds)
% hold on
% plot(200:1600,(200:1600)*table2array(mld.Coefficients(2,1))+table2array(mld.Coefficients(1,1)),'--')
% hold off
%% max speed order, fixed group size, ortho
[sortedMaxSpeed,ind] = sort(maxSpeed);
maxSpeedSortedTrials = copy(ind,:);
group_n = 20;
stds = NaN(1,group_n);
hitrate = NaN(1,group_n);
means = NaN(1,group_n);
speedStd = NaN(1,group_n);
speedPoint = NaN(1,group_n);
groupSize = length(copy)/group_n;
for i = 1:group_n
    stds(i) = std(maxSpeedSortedTrials(i*groupSize-groupSize+1:i*groupSize,29));
    speedPoint(i) = mean(sortedMaxSpeed(i*groupSize-groupSize+1:i*groupSize));
    hitrate(i) = sum(maxSpeedSortedTrials(i*groupSize-groupSize+1:i*groupSize,28));
    means(i) = mean(maxSpeedSortedTrials(i*groupSize-groupSize+1:i*groupSize,29));
    speedStd(i) = std(sortedMaxSpeed(i*groupSize-groupSize+1:i*groupSize));
end
figure
errorbar(speedPoint,stds,speedStd,'horizontal','o')
xlabel('Max Speed (mm/s)')
ylabel('Standard Error (mm)')
title('Max Speed vs Lateral Error')
mld = fitlm(speedPoint,stds)
hold on
plot(200:1600,(200:1600)*table2array(mld.Coefficients(2,1))+table2array(mld.Coefficients(1,1)),'--')
hold off
%%
plot(maxSpeed,copy(:,22),'o')

xlabel('Max Speed mm/s')
ylabel('Average Speed mm/s')
title('Max v.s. Average Speed, each point = one trial')
mld = fitlm(maxSpeed,copy(:,22))
%%
figure;
plot3(copy(:,15),copy(:,10),maxSpeed,'o')
grid on
xlabel('target size')
ylabel('distance')
zlabel('speed')