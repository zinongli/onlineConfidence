%% misc
speedRange = [100 1800];
penalty = 0.8; % secs
maxScore = 10;
sigma_s = 260; % sigma of actual max speed deviation from opt max speed
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
rep = 20;
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
yei = unifrnd(4,10);
mas = unifrnd(0.04,0.06);
mai = unifrnd(50,80);
UniLatParams = [xbs;xbi;ybs;ybi;xes;xei;yes;yei;mas;mai];
SAparams = [UniLatParams UniLatParams]; % placeholder for now
UniLatParams = [UniLatParams(1:2:7,:),UniLatParams(2:2:8,:)];
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
        StoSimData2(i,1:2) = - offset;
    end
    StoSimData1(i,3) = norm(offset);
    StoSimData2(i,3) = norm(offset);
    StoSimData1(i,4) = randsizes(i); % size of the small and solid target
    StoSimData2(i,4) = randsizes(i);
    optS_1 = findOptSpeed8(norm(offset),randsizes(i),totalTime,speedRange,SAparams,Left);
    time1 = norm(offset)/optS_1;
    remainTime1 = totalTime - time1;
    remainScore = (remainTime1/totalTime)*maxScore;
    xyBiasError1 = optS_1 .* UniLatParams(:,1) + UniLatParams(:,2);
    reg_phit1 = compute_phit0(randsizes(i), xyBiasError1);
    eGainReg1 = remainScore * reg_phit1;
    alt_phit1 = eGainReg1 * totalTime / ((totalTime - time1 - penalty) * maxScore);
    % alt_phit is the alternative probability required for switching to be beneficial
    if alt_phit1 > 0.99 % impossible for a profitable alternative as switch
        StoSimData1(i,5) = 0;
    else
        % model 1 size generation
        fun = @(x) abs(compute_phit0(x,xyBiasError1) - alt_phit1);
        x0 = randsizes(i);
        StoSimData1(i,5) = fminsearch(fun,x0); % size of the big and hollow target
        alt_scale1 = StoSimData1(i,5)./randsizes(i); % the scale of magnification
        % model 2 size generation and error checking
        optS_2 = findOptSpeed8(norm(offset),StoSimData1(i,5),totalTime,speedRange,SAparams,Left);
        time2 = norm(offset)/optS_2;
        remainTime2 = totalTime - time2 - penalty;
        remainScore2 = (remainTime2/totalTime)*maxScore;
        xyBiasError2 = optS_2 .* UniLatParams(:,1) + UniLatParams(:,2);
        reg_phit2 = compute_phit0(StoSimData1(i,5), xyBiasError2);
        eGainReg2 = remainScore2 * reg_phit2;
        alt_phit2 = eGainReg2 * totalTime / ((totalTime - time2) * maxScore);
        fun2 = @(x) abs(compute_phit0(x,xyBiasError2) - alt_phit2);
        x02 = randsizes(i);
        StoSimData2(i,5) = fminsearch(fun2,x02); % for error checking
        alt_scale2 = StoSimData2(i,5)./StoSimData1(i,5); % the scale of reduction
        mirrorTest(i) = alt_scale1 * alt_scale2; % the result should be as close to one as possible and is (derived small size / given small size)
        % model 1 endpoint gen
        actS_1 = normrnd(optS_1,sigma_s);
        StoSimData1(i,9) = actS_1;
        xyBiasErrorAct1 = actS_1 .* UniLatParams(:,1) + UniLatParams(:,2);
        alongError1(i) = normrnd(xyBiasErrorAct1(1),xyBiasErrorAct1(3));
        orthoError1(i) = normrnd(xyBiasErrorAct1(2),xyBiasErrorAct1(4));
        endPointErrorX1 = offset .* alongError1(i) / norm(offset);
        endPointErrorY1 = [offset(2),-offset(1)] .* orthoError1(i) / norm(offset);
        StoSimData1(i,6:7) = offset + endPointErrorX1 + endPointErrorY1;
        if norm([alongError1(i) orthoError1(i)]) > StoSimData1(i,5)
            % miss
            StoSimData1(i,8) = 0;
        elseif (StoSimData1(i,4) < norm([alongError1(i) orthoError1(i)])) && ((norm([alongError1(i) orthoError1(i)]) < StoSimData1(i,5)))
            % side hit
            StoSimData1(i,8) = 1;
        elseif norm([alongError1(i) orthoError1(i)]) < StoSimData1(i,4)
            % bullseye
            StoSimData1(i,8) = 2;
        end
        % model 2 endpoint gen
        actS_2 = normrnd(optS_2,sigma_s);
        StoSimData2(i,9) = actS_2;
        xyBiasErrorAct2 = actS_2 .* UniLatParams(:,1) + UniLatParams(:,2);
        alongError2(i) = normrnd(xyBiasErrorAct2(1),xyBiasErrorAct2(3));
        orthoError2(i) = normrnd(xyBiasErrorAct2(2),xyBiasErrorAct2(4));
        endPointErrorX2 = offset .* alongError2(i) / norm(offset);
        endPointErrorY2 = [offset(2),-offset(1)] .* orthoError2(i) / norm(offset);
        StoSimData2(i,6:7) = offset + endPointErrorX2 + endPointErrorY2;
        if norm([alongError2(i) orthoError2(i)]) > StoSimData1(i,5) 
            % miss
            StoSimData2(i,8) = 0;
        elseif (StoSimData1(i,4) < norm([alongError2(i) orthoError2(i)])) && ((norm([alongError2(i) orthoError2(i)]) < StoSimData1(i,5)))
            % side hit
            StoSimData2(i,8) = 1;
        elseif norm([alongError2(i) orthoError2(i)]) < StoSimData1(i,4)
            % bullseye
            StoSimData2(i,8) = 2;
        end
    end
    i
end
%%
compare = [StoSimData1(:,8),StoSimData2(:,8)];
figure;bar(sort(compare(:,1)))
figure;bar(sort(compare(:,2)))