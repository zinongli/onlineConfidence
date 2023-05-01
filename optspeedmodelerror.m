participant = 1;
if participant == 1 % JX
%     SAparams = [0.007439540690441, 10.919040445152977, 0.006120633676816, 4.625618271214958, 0.465, 72.52];
    SAparams =  [xParams(3:4,:);yParams(3:4,:);xParams(1:2,:);yParams(1:2,:)];
    SAparams(9,:) = 0.465;
    SAparams(10,:) = 72.52;
elseif participant == 0 % ZL
    SAparams = [0.00432, 6.61, 0.00423, 1.74, 0.480, 65.71];
end

speedRange = [100 1600];
optSpeed = NaN(1,360);


for i = 1:360
    if copy(i,8)<1000
        Left = 0;
    else
        Left = 1;
    end
    optSpeed(i) = findOptSpeed8(copy(i,10),copy(i,15),copy(i,3),speedRange,SAparams,Left);
    i
end

%%

figure
plot(optSpeed,maxSpeed,'o')
hold on
plot(1:1800,1:1800,'--')
hold off
xlabel('Optimal Max Speed (mm/s)','FontSize',18)
ylabel('Actual Max Speed (mm/s)','FontSize',18)
if participant == 1
    title('JX Maximum Speed Data vs Prediction','FontSize',18)
elseif participant == 0
    title('ZL Maximum Speed Data vs Prediction','FontSize',18)
end
    %%
figure
plot3(copy(:,10),copy(:,15),optSpeed,'or')
hold on
plot3(copy(:,10),copy(:,15),maxSpeed,'ob')
hold off
xlabel('target size')
ylabel('distance')
zlabel('speed')