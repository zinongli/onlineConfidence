totalTime = 3; % secs
speedRange = [100 1600];


participant = 0;
if participant == 1 % JX
    SAparams = [0.00744, 10.9, 0.00570, 3.41, 0.465, 72.52];
elseif participant == 0 % ZL
    SAparams = [0.00432, 6.61, 0.00423, 1.74, 0.480, 65.71];
end

speedResolution = 2000;
speeds = linspace(speedRange(1),speedRange(2),speedResolution);

tSizeScale = 2;
a = SAparams(1);
b = SAparams(2);
c = SAparams(3);
d = SAparams(4);
e = SAparams(5);
f = SAparams(6);

distance = 90:10:350; %mm
penalty = 0.2; % secs
target = 2:0.1:6; %mm
%%
points1 = NaN(length(target),length(distance),length(speeds));
points2 = NaN(length(target),length(distance),length(speeds));
alt_phit = NaN(length(target),length(distance));
alt_scale = NaN(length(target),length(distance));
opt_speed = NaN(length(target),length(distance));
for tt = 1:length(target) %loop over all target sizes
    target1 = target(tt);
    target2 = target1 * tSizeScale; %20% increase
    
    for dd = 1:length(distance) %loop over all distances
        
        for ss = 1:length(speeds) %loop over all speeds (incorporates errors)
            avgSpeed = speeds(ss) * e + f;
            time = distance(dd)/avgSpeed;
            remainTime = totalTime - time;
            remainScore = (remainTime/totalTime)*10;
            remainScoreSwitch = ((remainTime - penalty)/totalTime)*10;
            errorx = speeds(ss) * a + b;
            errory = speeds(ss) * c + d;
            probHit1 = compute_phit2(target1,errorx,errory);
            probHit2 = compute_phit2(target2,errorx,errory);
            
            points1(tt,dd,ss) = remainScore * probHit1;
            points2(tt,dd,ss) = remainScoreSwitch * probHit2;
        end
        [value,ind] = max(points1(tt,dd,:));
        time = distance(dd)/speeds(ind);
        errorx = speeds(ind) * a + b;
        alt_phit(tt,dd) = value / (totalTime - time - penalty); 
        % alt_phit is the alternative probability required for switching to be beneficial
        opt_speed(tt,dd) = speeds(ind);
        if alt_phit(tt,dd) < 0.99
            fun = @(x) abs(compute_phit2(x,errorx,errory) - alt_phit(tt,dd));
            x0 = target(tt);
            x = fminsearch(fun,x0);
            alt_scale(tt,dd) = x./target(tt);
        end
        dd
    end
    tt
end
% %%
% figure
% for ii = 1:length(speeds)
% hold on
% surf(target,distance,points1(:,:,ii)','FaceAlpha',0,'EdgeColor','b')
% surf(target,distance,points2(:,:,ii)','FaceAlpha',0,'EdgeColor','r')
% hold off
% xlabel('target size')
% ylabel('distance')
% zlabel('points')
% view(3)
% pause(.5)
% end
%%
figure
surf(target,distance,alt_phit')
xlabel('target size')
ylabel('distance')
zlabel('P(Hit|Switch)')
title('The P(Hit|Switch) required for switching to be equally desirable')
hold on 
surf(target,distance,repmat([1],length(target),length(distance))','EdgeColor','r')
hold off
%%
figure
surf(target,distance,alt_scale')
xlabel('target size')
ylabel('distance')
zlabel('scale')
title('The magnification scale required for switching to be equally desirable')
%%
figure
surf(target,distance,opt_speed')
xlabel('target size')
ylabel('distance')
zlabel('speed')
title('Optimal Max Speed')