%% misc (arbitary values)
speedRange = [100 1800];
penalty = 0.8; % secs
maxScore = 10;
%% Task Env Generation
dists_n = 3;
proj2tablet = 1.886;
pixellength = 0.248;
proj2mm = proj2tablet .* pixellength;
UniRandRadius = 70 .* proj2mm;
edgesize = 50;
sizes_n = 6;
WindowWidth = 1024; % projector window width
yCenter = 384; % projector screen y center
rep = 10;
totalTime = 3;

distances = linspace(edgesize,WindowWidth-edgesize,dists_n+2)-edgesize;
distances = repmat(distances(2:end-1),1,sizes_n*rep) .* proj2mm;


target_sizes = [5,10,15,20,25,30] ./ pixellength ./ proj2tablet;
target_sizes = repmat(target_sizes,1,dists_n*rep);
target_sizes = target_sizes';
target_sizes = target_sizes(:)';

%% Stochastic Section
seeds = [randperm(size(distances,2)), randperm(size(distances,2))];
randdists = distances(seeds);
randdists = randdists(:);
randsizes = target_sizes(seeds);
randsizes = randsizes(:);
xbs = unifrnd(-0.03,0.03);
xbi = unifrnd(-25,25);
ybs = unifrnd(-0.02,0.02);
ybi = unifrnd(-8,8);
xes = unifrnd(0.002,0.02);
xei = unifrnd(4,12);
yes = unifrnd(0.002,0.02);
yei = unifrnd(4,12);
mas = unifrnd(0.4,0.6);
mai = unifrnd(50,80);
mes = unifrnd(0.0001,0.05);
mei = unifrnd(10,30);
sigma_s = unifrnd(10,45); % sigma of actual max speed deviation from opt max speed
UniLatParams = [xbs;xbi;ybs;ybi;xes;xei;yes;yei;mas;mai;mes;mei];
SAparams = [UniLatParams UniLatParams]; % placeholder for now
UniLatParams = [UniLatParams(1:2:11,:),UniLatParams(2:2:12,:)];
Left = 0;
StoSimData1 = NaN(length(randsizes),9);
StoSimData2 = NaN(length(randsizes),9);
mirrorTest = NaN(length(randsizes),1);
alongError1 = NaN(length(randsizes),1);
orthoError1 = NaN(length(randsizes),1);
alongError2 = NaN(length(randsizes),1);
orthoError2 = NaN(length(randsizes),1);
for i = 1:length(randsizes)
    if rem(i,2)
        theta = -pi/12 + (pi/6) * rand(1);
        rho = randdists(i)-UniRandRadius + UniRandRadius * 2 * rand(1);
        [offset(1),offset(2)] = pol2cart(theta,rho);
        StoSimData1(i,1:2) = - offset;
        StoSimData2(i,1:2) = - offset;
    else
        theta = -pi/12 + (pi/6) * rand(1);
        rho = randdists(i)-UniRandRadius + UniRandRadius * 2 * rand(1);
        [offset(1),offset(2)] = pol2cart(theta,rho);
        StoSimData1(i,1:2) = offset;
        StoSimData2(i,1:2) = offset;
    end
    StoSimData1(i,3) = norm(offset);
    StoSimData2(i,3) = norm(offset);
    StoSimData1(i,4) = randsizes(i); % size of the small and solid target
    
    optS_1 = findOptSpeed8(norm(offset),randsizes(i),totalTime,speedRange,SAparams,Left);
    time1 = norm(offset)/(UniLatParams(5,1) * optS_1 + UniLatParams(5,2));
    remainTime1 = totalTime - time1;
    remainScore = (remainTime1/totalTime)*maxScore;
    xyBiasError1 = optS_1 .* UniLatParams(:,1) + UniLatParams(:,2);
    reg_phit1 = compute_phit0(randsizes(i), xyBiasError1);
    eGainReg1 = remainScore * reg_phit1;
    alt_phit1 = eGainReg1 * totalTime / ((totalTime - time1 - penalty) * maxScore);
    % alt_phit is the alternative probability required for switching to be beneficial
    if alt_phit1 > 0.99 || alt_phit1 <= 0% impossible for a profitable alternative as switch
        StoSimData1(i,5) = 0;
    else
        % model 1 size generation
        fun = @(x) abs(compute_phit0(x,xyBiasError1) - alt_phit1);
        x0 = randsizes(i);
        StoSimData1(i,5) = fmincon(fun,x0,[],[],[],[],5,80); % size of the big and hollow target
        StoSimData2(i,4) = StoSimData1(i,5);
        alt_scale1 = StoSimData1(i,5)./randsizes(i); % the scale of magnification
        % model 2 size generation and error checking
        optS_2 = findOptSpeed8(norm(offset),StoSimData1(i,5),totalTime,speedRange,SAparams,Left);
        time2 = norm(offset)/(UniLatParams(5,1) * optS_2 + UniLatParams(5,2));
        remainTime2 = totalTime - time2 - penalty;
        remainScore2 = (remainTime2/totalTime)*maxScore;
        xyBiasError2 = optS_2 .* UniLatParams(:,1) + UniLatParams(:,2);
        reg_phit2 = compute_phit0(StoSimData1(i,5), xyBiasError2);
        eGainReg2 = remainScore2 * reg_phit2;
        alt_phit2 = eGainReg2 * totalTime / ((totalTime - time2) * maxScore);
        fun2 = @(x) abs(compute_phit0(x,xyBiasError2) - alt_phit2);
        x02 = randsizes(i);
        StoSimData2(i,5) = fmincon(fun2,x02,[],[],[],[],1,50); % for error checking
        alt_scale2 = StoSimData2(i,5)./StoSimData1(i,5); % the scale of reduction
        mirrorTest(i) = alt_scale1 * alt_scale2; % the result should be as close to one as possible and is (derived small size / given small size)
        % model 1 endpoint gen
        actS_1 = normrnd(optS_1,sigma_s);
        if actS_1 < speedRange(1)
            actS_1 = speedRange(1);
        end
        StoSimData1(i,9) = actS_1;
        StoSimData1(i,10) = normrnd(UniLatParams(5,1) * actS_1 + UniLatParams(5,2),UniLatParams(6,1) * actS_1 + UniLatParams(6,2));
        xyBiasErrorAct1 = actS_1 .* UniLatParams(:,1) + UniLatParams(:,2);
        alongError1 = normrnd(xyBiasErrorAct1(1),xyBiasErrorAct1(3));
        orthoError1 = normrnd(xyBiasErrorAct1(2),xyBiasErrorAct1(4));
%         endPointErrorX1 = offset .* alongError1 / norm(offset);
%         endPointErrorY1 = [offset(2),-offset(1)] .* orthoError1 / norm(offset);
%         StoSimData1(i,6:7) = offset + endPointErrorX2 + endPointErrorY2;
        StoSimData1(i,6:7) = [alongError1, orthoError1];
        if norm([alongError1 orthoError1]) > StoSimData1(i,5)
            % miss
            StoSimData1(i,8) = 0;
        elseif (StoSimData1(i,4) < norm([alongError1 orthoError1])) && ((norm([alongError1 orthoError1]) < StoSimData1(i,5)))
            % side hit
            StoSimData1(i,8) = 1;
        elseif norm([alongError1 orthoError1]) < StoSimData1(i,4)
            % bullseye
            StoSimData1(i,8) = 2;
        end
        % model 2 endpoint gen
        actS_2 = normrnd(optS_2,sigma_s);
        if actS_2 < speedRange(1)
            actS_2 = speedRange(1);
        end
        StoSimData2(i,9) = actS_2;
        StoSimData2(i,10) = normrnd(UniLatParams(5,1) * actS_2 + UniLatParams(5,2),UniLatParams(6,1) * actS_2 + UniLatParams(6,2));
        xyBiasErrorAct2 = actS_2 .* UniLatParams(:,1) + UniLatParams(:,2);
        alongError2 = normrnd(xyBiasErrorAct2(1),xyBiasErrorAct2(3));
        orthoError2 = normrnd(xyBiasErrorAct2(2),xyBiasErrorAct2(4));
%         endPointErrorX2 = offset .* alongError2 / norm(offset);
%         endPointErrorY2 = [offset(2),-offset(1)] .* orthoError2 / norm(offset);
%         StoSimData2(i,6:7) = offset + endPointErrorX2 + endPointErrorY2;
        StoSimData2(i,6:7) = [alongError2 , orthoError2];
        if norm([alongError2 orthoError2]) > StoSimData1(i,5) 
            % miss
            StoSimData2(i,8) = 0;
        elseif (StoSimData1(i,4) < norm([alongError2 orthoError2])) && ((norm([alongError2 orthoError2]) < StoSimData1(i,5)))
            % body shot
            StoSimData2(i,8) = 1;
        elseif norm([alongError2 orthoError2]) < StoSimData1(i,4)
            % bullseye
            StoSimData2(i,8) = 2;
        end
    end
    i
end
%% elimination
vStoSimData1 = StoSimData1(StoSimData1(:,5) ~= 0,:);
vStoSimData2 = StoSimData2(StoSimData1(:,5) ~= 0,:);

%% table of content
% 1,2 = subjective target location
% 3 = target distance
% 4 = target size
% 5 = alternative target size
% 6,7 = along & orthogonal error
% 8 = hit state ( 0 = miss, 1 = body shot, 2 = bullseye)
% 9 = max speed
% 10 = avg speed
%% retrieval

xTotal = [vStoSimData1(:,6);vStoSimData2(:,6)];
% xRight1 = vStoSimData1(vStoSimData1(:,1)>0,6);
% xLeft1 = vStoSimData1(vStoSimData1(:,1)<0,6);


theta0 = [0,0,0.0001,5];
UB = [1,30,1,50];
LB = [-1,-30,0.0001,1];

f = @(theta,speed,error) -(log(1/sqrt(2*pi)) * length(speed) + sum(-log(theta(3).*speed + theta(4)) - (((error - (theta(1).*speed + theta(2))).^2)./ (2.*(theta(3).*speed+theta(4)).^2))));


xParams = NaN(3,4);
fun = @(theta) f(theta,[vStoSimData1(:,9);vStoSimData2(:,9)],xTotal);
xParams(1,:) = bads(fun,theta0,LB,UB);
% fun = @(theta) f(theta,vStoSimData1(vStoSimData1(:,1)>0,9),xRight1);
% xParams1(2,:) = bads(fun,theta0,LB,UB);
% fun = @(theta) f(theta,vStoSimData1(vStoSimData1(:,1)<0,9),xLeft1);
% xParams1(3,:) = bads(fun,theta0,LB,UB);


% xTotal2 = vStoSimData2(:,6);
% xRight2 = vStoSimData2(vStoSimData2(:,1)>0,6);
% xLeft2 = vStoSimData2(vStoSimData2(:,1)<0,6);
% 
% 
% theta0 = [0,0,0.0001,5];
% UB = [1,30,1,50];
% LB = [-1,-30,0.0001,1];
% 
% f = @(theta,speed,error) -(log(1/sqrt(2*pi)) * length(speed) + sum(-log(theta(3).*speed + theta(4)) - (((error - (theta(1).*speed + theta(2))).^2)./ (2.*(theta(3).*speed+theta(4)).^2))));
% 

% xParams2 = NaN(3,4);
% fun = @(theta) f(theta,vStoSimData2(:,9),xTotal2);
% xParams2(1,:) = bads(fun,theta0,LB,UB);
% fun = @(theta) f(theta,vStoSimData2(vStoSimData2(:,1)>0,9),xRight2);
% xParams2(2,:) = bads(fun,theta0,LB,UB);
% fun = @(theta) f(theta,vStoSimData2(vStoSimData2(:,1)<0,9),xLeft2);
% xParams2(3,:) = bads(fun,theta0,LB,UB);


yTotal = [vStoSimData1(:,7);vStoSimData2(:,7)];
% yRight1 = vStoSimData1(vStoSimData1(:,1)>0,7);
% yLeft1 = vStoSimData1(vStoSimData1(:,1)<0,7);


theta0 = [0,0,0.0001,5];
UB = [1,30,1,50];
LB = [-1,-30,0.0001,1];

f = @(theta,speed,error) -(log(1/sqrt(2*pi)) * length(speed) + sum(-log(theta(3).*speed + theta(4)) - (((error - (theta(1).*speed + theta(2))).^2)./ (2.*(theta(3).*speed+theta(4)).^2))));


yParams = NaN(3,4);
fun = @(theta) f(theta,[vStoSimData1(:,9);vStoSimData2(:,9)],yTotal);
yParams(1,:) = bads(fun,theta0,LB,UB);
% fun = @(theta) f(theta,vStoSimData1(vStoSimData1(:,1)>0,9),yRight1);
% yParams1(2,:) = bads(fun,theta0,LB,UB);
% fun = @(theta) f(theta,vStoSimData1(vStoSimData1(:,1)<0,9),yLeft1);
% yParams1(3,:) = bads(fun,theta0,LB,UB);

% yTotal2 = vStoSimData2(:,7);
% yRight2 = vStoSimData2(vStoSimData2(:,1)>0,7);
% yLeft2 = vStoSimData2(vStoSimData2(:,1)<0,7);
% 
% 
% theta0 = [0,0,0.0001,5];
% UB = [1,30,1,50];
% LB = [-1,-30,0.0001,1];
% 
% f = @(theta,speed,error) -(log(1/sqrt(2*pi)) * length(speed) + sum(-log(theta(3).*speed + theta(4)) - (((error - (theta(1).*speed + theta(2))).^2)./ (2.*(theta(3).*speed+theta(4)).^2))));
% 
% 
% yParams2 = NaN(3,4);
% fun = @(theta) f(theta,vStoSimData2(:,9),yTotal2);
% yParams2(1,:) = bads(fun,theta0,LB,UB);
% fun = @(theta) f(theta,vStoSimData2(vStoSimData2(:,1)>0,9),yRight2);
% yParams2(2,:) = bads(fun,theta0,LB,UB);
% fun = @(theta) f(theta,vStoSimData2(vStoSimData2(:,1)<0,9),yLeft2);
% yParams2(3,:) = bads(fun,theta0,LB,UB);
mld1 = fitlm([vStoSimData1(:,9);vStoSimData2(:,9)],[vStoSimData1(:,10);vStoSimData2(:,10)]);

theta0 = [table2array(mld1.Coefficients(2,1)),table2array(mld1.Coefficients(1,1)),0.0001,30];
UB = [5,100,5,100];
LB = [0,-100,0.0001,1];

f = @(theta,speed,error) -(log(1/sqrt(2*pi)) * length(speed) + sum(-log(theta(3).*speed + theta(4)) - (((error - (theta(1).*speed + theta(2))).^2)./ (2.*(theta(3).*speed+theta(4)).^2))));
fun = @(theta) f(theta,[vStoSimData1(:,9);vStoSimData2(:,9)],[vStoSimData1(:,10);vStoSimData2(:,10)]);
changingSig = bads(fun,theta0,LB,UB);

ReUniLatParams = [xParams(1,1:2);yParams(1,1:2);xParams(1,3:4);yParams(1,3:4);changingSig(1:2);changingSig(3:4)]
UniLatParams
cache = ReUniLatParams';
ReSAparams = [cache(:),cache(:)];
%% opt speed prediction
ReConcat = [vStoSimData1;vStoSimData2];
ReOptSpeed8 = NaN(1,length(ReConcat));

for i = 1:length(ReConcat)
%     Left = ReConcat(i,1)<0;
    ReOptSpeed8(i) = findOptSpeed8(ReConcat(i,3),ReConcat(i,4),totalTime,speedRange,ReSAparams,0);

end
% %%
% sss = @(sigma,m,o) log(1/sqrt(2*pi)) * length(m) - (-log(sigma) * length(m) + sum(-((m'-o).^2)./(2 * sigma.^2)));
% sfun = @(sigma) sss(sigma,ReConcat(:,9),ReOptSpeed8);
% sigma_fit = bads(sfun,0.001,0.001,300);
% sss(sigma_fit,ReConcat(:,9),ReOptSpeed8);
%%
pb = makedist('Normal');
qqplot(ReConcat(:,9) - ReOptSpeed8',pb)
%%
[a,b,c] = swtest(ReConcat(:,9) - ReOptSpeed8')
%%
set(groot,'defaultAxesFontSize',18)
figure
plot(ReConcat(:,9),ReConcat(:,10),'o')
hold on
plot(100:1800,100:1800,'--')
hold off
xlabel('Optimal Max Speed (mm/s)','FontSize',18)
ylabel('Actual Max Speed (mm/s)','FontSize',18)
