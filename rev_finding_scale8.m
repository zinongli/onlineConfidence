totalTime = 3; % secs
speedRange = [100 1600];

Left = 0;
UniLatParams = SAparams(:,Left + 1); %therefore SAparams is right column | left column

distance = 80:20:370; %mm
penalty = 0.8; % secs
target = 2:2:40; %mm
maxScore = 10;
%%

rev_reg_phit = NaN(length(target),length(distance));
rev_alt_phit = NaN(length(target),length(distance));
rev_alt_scale = NaN(length(target),length(distance));
rev_opt_speed = NaN(length(target),length(distance));
for tt = 1:length(target) %loop over all target sizes
    for dd = 1:length(distance) %loop over all distances
        target_ij = target(tt);
        optS_i = findOptSpeed8(distance(dd),target_ij,totalTime,speedRange,SAparams,Left);
        time = distance(dd)/optS_i;
        remainTime = totalTime - time -penalty;
        remainScore = (remainTime/totalTime)*maxScore;
        biasx = optS_i * UniLatParams(1) + UniLatParams(2);
        biasy = optS_i * UniLatParams(3) + UniLatParams(4);
        errorx = optS_i * UniLatParams(5) + UniLatParams(6);
        errory = optS_i * UniLatParams(7) + UniLatParams(8);
        rev_reg_phit(tt,dd) = compute_phit0(target_ij,errorx,errory, biasx, biasy);
        eGainReg = remainScore * rev_reg_phit(tt,dd);
        rev_alt_phit(tt,dd) = eGainReg * totalTime / ((totalTime - time) * maxScore);
        % alt_phit is the alternative probability required for switching to be beneficial
        rev_opt_speed(tt,dd) = optS_i;
%         if alt_phit(tt,dd) < 0.99
            fun = @(x) abs(compute_phit0(x,errorx,errory, biasx, biasy) - rev_alt_phit(tt,dd));
            x0 = target_ij;
%             x = bads(fun,x0,0,50);
            x = fminsearch(fun,x0);
            rev_alt_scale(tt,dd) = x./target_ij;
%         end
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
subplot(1,4,1)

surf(target,distance,rev_reg_phit', 'EdgeAlpha',0)
hold on
% surf(target,distance,alt_phit2', 'EdgeAlpha',0)
surf(target,distance,repmat([1],length(target),length(distance))','EdgeAlpha',0, 'FaceColor',[0.8,0.8,0.8])
hold off
grid on
xlabel('target size')
ylabel('distance')
zlabel('P(Hit|Switch)')
title('The P(Hit|Switch) of Normal Target Size')

subplot(1,4,2)
surf(target,distance,rev_alt_phit', 'EdgeAlpha',0)
hold on
% surf(target,distance,alt_phit2', 'EdgeAlpha',0)
surf(target,distance,repmat([1],length(target),length(distance))','EdgeAlpha',0, 'FaceColor',[0.8,0.8,0.8])
hold off
grid on
xlabel('target size')
ylabel('distance')
zlabel('P(Hit|Switch)')
title('The P(Hit|Switch) Required for Switching to be Equally Desirable')

subplot(1,4,3)
surf(target,distance,rev_alt_scale', 'EdgeAlpha',0)
% hold on
% surf(target,distance,alt_scale2', 'EdgeAlpha',0)
% hold off
grid on
xlabel('target size')
ylabel('distance')
zlabel('scale')
title('The Magnification Scale Required for Switching to be Equally Desirable')


subplot(1,4,4)
surf(target,distance,rev_opt_speed', 'EdgeAlpha',0)
% hold on3
% surf(target,distance,opt_speed2', 'EdgeAlpha',0)
% hold off
grid on
xlabel('target size')
ylabel('distance')
zlabel('speed')
title('Optimal Max Speed')
