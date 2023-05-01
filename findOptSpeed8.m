function opt_speed = findOptSpeed8(distance,tSize,totalTime,speedRange,SAparams,Left)
speedResolution = 2000;
speeds = linspace(speedRange(1),speedRange(2),speedResolution);
penalty = 0.2;
tSizeScale = 1.2;

points1 = NaN(1,length(speeds));
points2 = NaN(1,length(speeds));

params = SAparams(:,Left + 1); %therefore SAparams is right column | left column

for ss = 1:length(speeds) %loop over all speeds (incorporates errors)
    avgSpeed = speeds(ss) * params(9) + params(10);
    time = distance/avgSpeed;
    remainTime = totalTime - time;
    remainScore = (remainTime/totalTime)*10;
    remainScoreSwitch = ((remainTime - penalty)/totalTime)*10;
    
    errorx = speeds(ss) * params(1) + params(2);
    errory = speeds(ss) * params(3) + params(4);
    biasx = speeds(ss) * params(5) + params(6);
    biasy = speeds(ss) * params(7) + params(8);
    
    probHit1 = compute_phit8(tSize,errorx,errory, biasx, biasy);
    probHit2 = compute_phit8(tSize.*tSizeScale,errorx,errory, biasx, biasy);
    
    points1(ss) = remainScore * probHit1;
    points2(ss) = remainScoreSwitch * probHit2;
end
[~,ind] = max(points1);
opt_speed = speeds(ind);
end