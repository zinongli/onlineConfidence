totalTime = 3; % secs
speedRange = [100 1600];


participant = 1;
if participant == 1 % JX
    SAparams = [0.007439540690441, 10.919040445152977, 0.006120633676816, 4.625618271214958, 0.465, 72.52];
elseif participant == 0 % ZL
    SAparams = [0.00432, 6.61, 0.00423, 1.74, 0.480, 65.71]; %group_n = 20
end

speedResolution = 100;
speeds = linspace(speedRange(1),speedRange(2),speedResolution);

tSizeScale = 0.5;
a = SAparams(1);
b = SAparams(2);
c = SAparams(3);
d = SAparams(4);
e = SAparams(5);
f = SAparams(6);

distance = 80:10:370; %mm
penalty = 0.2; % secs
target = 2:1:20; %mm
%%
points1 = NaN(length(target),length(distance),length(speeds));
points2 = NaN(length(target),length(distance),length(speeds));
alt_phit = NaN(length(target),length(distance));
alt_scale = NaN(length(target),length(distance));
opt_speed = NaN(length(target),length(distance));
for tt = 1:length(target) %loop over all target sizes
    target1 = target(tt);
    target2 = target1 * tSizeScale; 
    
    for dd = 1:length(distance) %loop over all distances
        
        for ss = 1:length(speeds) %loop over all speeds (incorporates errors)
            avgSpeed = speeds(ss) * e + f;
            time = distance(dd)/avgSpeed;
            remainTime = totalTime - time;
            remainScore = (remainTime/totalTime)*10;
            remainScoreSwitch = ((remainTime + penalty)/totalTime)*10;
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
        errory = speeds(ind) * c + d;
        alt_phit(tt,dd) = value / (totalTime - time + penalty); 
        % alt_phit is the alternative probability required for switching to be beneficial
        opt_speed(tt,dd) = speeds(ind);
        if alt_phit(tt,dd) < 0.99
            fun = @(x) abs(compute_phit2(x,errorx,errory) - alt_phit(tt,dd));
            x0 = target(tt);
            x = fminsearch(fun,x0);
            alt_scale(tt,dd) = x./target(tt);
        end
        
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
surf(target,distance,alt_phit', 'EdgeAlpha',0)
hold on
% surf(target,distance,alt_phit2', 'EdgeAlpha',0)
surf(target,distance,repmat([1],length(target),length(distance))','EdgeAlpha',0, 'FaceColor',[0.8,0.8,0.8])
hold off
grid on
xlabel('target size')
ylabel('distance')
zlabel('P(Hit|Switch)')
title('The P(Hit|Switch) required for switching to be equally desirable')

%%
figure
surf(target,distance,alt_scale', 'EdgeAlpha',0)
% hold on
% surf(target,distance,alt_scale2', 'EdgeAlpha',0)
% hold off
grid on
xlabel('target size')
ylabel('distance')
zlabel('scale')
title('The magnification scale required for switching to be equally desirable')
%%
figure

surf(target,distance,opt_speed', 'EdgeAlpha',0)
% hold on3
% surf(target,distance,opt_speed2', 'EdgeAlpha',0)
% hold off
grid on
xlabel('target size')
ylabel('distance')
zlabel('speed')
title('Optimal Max Speed')
