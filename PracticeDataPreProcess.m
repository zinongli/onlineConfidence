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
%% Preprocess
subj = 'JX';
try
    load('data_onlineConf/JX/JX_preProcessed.mat')
    disp('>>>>>>>>>Pre-processed data exists.<<<<<<<<<')
catch
    load('data_onlineConf/JX/JM_practice_traj_S1_06-Apr-2023_tform.mat')
    load('data_onlineConf/JX/JM_practice_traj_S1c_06-Apr-2023_rawtotal.mat')
    load('data_onlineConf/JX/JM_practice_traj_S1c_06-Apr-2023_traXtotal.mat')
    load('data_onlineConf/JX/JM_practice_traj_S1c_06-Apr-2023_traYtotal.mat')

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
    
    % sigmoid fitting max speed from trajs
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
    % feel free to comment this section as it's only for error checking
    % randomly select 10 trials to see how the fit went
    figure(1)
    for i = 1:10
        trial_i = randi(length(copy),1);
        x = 50:endTime(trial_i);
        y = traProjection(trial_i,x);
        x = x - 50;
        b = sigmoidFit(:,trial_i);
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
        xlabel('Recentered Time (s)')
        ylabel('Standardized Trajectory (mm)')
        title(['Trial #' num2str(trial_i)])
        pause(0.2)
    end
    copy(:,30) = -sigmoidFit(2,:)'.*0.25.*copy(:,10).*60;
    save('data_onlineConf/JX/JX_preProcessed.mat','copy');
end
%% linear fit max & avg as param_0
figure(2)
plot(copy(:,30),copy(:,22),'o')
xlabel('Fit Max Speed (mm/s)')
ylabel('Average Speed (mm/s)')
title('Max v.s. Average Speed, each point = one trial')

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

figure(3)
for i = 1:3
    plot(tSizes,SizeDistHitRate(:,i),'-o')
    hold on
end
plot(tSizes,0.3:0.1:0.7,'r--') % Thjs so called Expectation is ARBITARY
hold off
xlabel('Target Size (mm)')
ylabel('Hit Rate')
legend('~100 mm','~200 mm','~300 mm','Expectation','Location','Northwest')
title('Hit Rates of Each Distance Range')

