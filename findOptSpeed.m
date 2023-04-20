function opt_speed = findOptSpeed(distance,tSize,totalTime,speedRange,SAparams)
speedResolution = 2000;
speeds = linspace(speedRange(1),speedRange(2),speedResolution);
penalty = 0.2;
tSizeScale = 1.2;
a = SAparams(1);
b = SAparams(2);
c = SAparams(3);
d = SAparams(4);
e = SAparams(5);
f = SAparams(6);
points1 = NaN(1,length(speeds));
points2 = NaN(1,length(speeds));

for ss = 1:length(speeds) %loop over all speeds (incorporates errors)
    avgSpeed = speeds(ss) * e + f;
    time = distance/avgSpeed;
    remainTime = totalTime - time;
    remainScore = (remainTime/totalTime)*10;
    remainScoreSwitch = ((remainTime - penalty)/totalTime)*10;
    errorx = speeds(ss) * a + b;
    errory = speeds(ss) * c + d;
    probHit1 = compute_phit2(tSize,errorx,errory);
    probHit2 = compute_phit2(tSize.*tSizeScale,errorx,errory);
    
    %this simple point calculation only works because using 5-seconds.
    %If another time is used a more complex points calculation is
    %needed.
    points1(ss) = remainScore * probHit1;
    points2(ss) = remainScoreSwitch * probHit2;
end
[~,ind] = max(points1);
opt_speed = speeds(ind);
end