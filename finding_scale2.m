totalTime = 3; % secs
speedRange = [100 1600];






speedResolution = 200;
speeds = linspace(speedRange(1),speedRange(2),speedResolution);
tSizeScale = 1.5;


Left = 0;
UniLatParams = SAparams(:,Left + 1); %therefore SAparams is right column | left column

distance = 80:20:370; %mm
penalty = 0.2; % secs
target = 2:2:20; %mm
maxScore = 10;
%%
points1 = NaN(length(target),length(distance),length(speeds));
points2 = NaN(length(target),length(distance),length(speeds));
alt_phit = NaN(length(target),length(distance));
alt_scale = NaN(length(target),length(distance));
opt_speed = NaN(length(target),length(distance));
for tt = 1:length(target) %loop over all target sizes
    target1 = target(tt);
    target2 = target1.*tSizeScale;
    
    for dd = 1:length(distance) %loop over all distances

        target_i = target1;
        
        optS_i = findOptSpeed8(distance(dd),target_i,totalTime,speedRange,SAparams,Left);
        time = distance(dd)/optS_i;
        remainTime = totalTime - time;
        remainScore = (remainTime/totalTime)*maxScore;
        biasx = optS_i * UniLatParams(1) + UniLatParams(2);
        biasy = optS_i * UniLatParams(3) + UniLatParams(4);
        errorx = optS_i * UniLatParams(5) + UniLatParams(6);
        errory = optS_i * UniLatParams(7) + UniLatParams(8);
        eGainReg = remainScore * compute_phit8(target_i,errorx,errory, biasx, biasy);
        alt_phit(tt,dd) = eGainReg * totalTime / ((totalTime - time - penalty) * maxScore);
        % alt_phit is the alternative probability required for switching to be beneficial
        opt_speed(tt,dd) = optS_i;
        if alt_phit(tt,dd) < 0.99
            fun = @(x) abs(compute_phit8(x,errorx,errory, biasx, biasy) - alt_phit(tt,dd));
            x0 = target(tt);
            x = bads(fun,x0,0,50);
%             x = fminsearch(fun,x0);
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
