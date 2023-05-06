totalTime = 3; % secs
speedRange = [100 1600];

Left = 0;
UniLatParams = SAparams(:,Left + 1); %therefore SAparams is right column | left column

distance = 80:20:370; %mm
penalty = 0.5; % secs
target = [2:2:40] .* 1.2; %mm
maxScore = 10;
%%
points1 = NaN(length(target),length(distance),length(speeds));
points2 = NaN(length(target),length(distance),length(speeds));
reg_phit = NaN(length(target),length(distance));
alt_phit = NaN(length(target),length(distance));
alt_scale = NaN(length(target),length(distance));
opt_speed = NaN(length(target),length(distance));
for tt = 1:length(target) %loop over all target sizes
    for dd = 1:length(distance) %loop over all distances
        optS_i = findOptSpeed8(distance(dd),target(tt),totalTime,speedRange,SAparams,Left);
        time = distance(dd)/optS_i;
        remainTime = totalTime - time - penalty;
        remainScore = (remainTime/totalTime)*maxScore;
        biasx = optS_i * UniLatParams(1) + UniLatParams(2);
        biasy = optS_i * UniLatParams(3) + UniLatParams(4);
        errorx = optS_i * UniLatParams(5) + UniLatParams(6);
        errory = optS_i * UniLatParams(7) + UniLatParams(8);
        reg_phit(tt,dd) = compute_phit8(target(tt),errorx,errory, biasx, biasy);
        eGainReg = remainScore * reg_phit(tt,dd);
        alt_phit(tt,dd) = eGainReg * totalTime / ((totalTime - time) * maxScore);
        % alt_phit is the alternative probability required for switching to be beneficial
        opt_speed(tt,dd) = optS_i;
        if alt_phit(tt,dd) < 0.99
            fun = @(x) abs(compute_phit8(x,errorx,errory, biasx, biasy) - alt_phit(tt,dd));
            x0 = target(tt);
%             x = bads(fun,x0,0,50);
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
set(groot,'defaultAxesFontSize',18)
figure
subplot(1,3,1)
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


subplot(1,3,2)
surf(target,distance,alt_scale', 'EdgeAlpha',0)
% hold on
% surf(target,distance,alt_scale2', 'EdgeAlpha',0)
% hold off
grid on
xlabel('target size')
ylabel('distance')
zlabel('scale')
title('The magnification scale required for switching to be equally desirable')

subplot(1,3,3)
surf(target,distance,opt_speed', 'EdgeAlpha',0)
% hold on3
% surf(target,distance,opt_speed2', 'EdgeAlpha',0)
% hold off
grid on
xlabel('target size')
ylabel('distance')
zlabel('speed')
title('Optimal Max Speed')