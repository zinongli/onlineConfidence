%% Participant Generation
abs = unifrnd(-0.03,0.03);
abi = unifrnd(-25,25);
obs = unifrnd(-0.02,0.02);
obi = unifrnd(-8,8);
aes = unifrnd(0.002,0.02);
aei = unifrnd(4,12);
oes = unifrnd(0.002,0.02);
oei = unifrnd(4,10);
mas = unifrnd(0.04,0.06);
mai = unifrnd(50,80);

UniLatParams = [abs;abi;obs;aes;aei;oes;mas;mai];

%% misc
SAparams = [UniLatParams UniLatParams]; % placeholder for now
Left = 0;
speedRange = [100 1800];
penalty = 0.8; % secs
maxScore = 10;
%% Task Env Generation
dists_n = 3;
proj2tablet = 1.886;
proj2mm = proj2tablet .* pixellength;
UniRandRadius = 70 .* proj2mm;
edgesize = 50;
sizes_n = 6;
WindowWidth = 1024; % projector window width
yCenter = 384; % projector screen y center
rep = 2;
totalTime = 3;

distances = linspace(edgesize,WindowWidth-edgesize,dists_n+2)-edgesize;
distances = repmat(distances(2:end-1),1,sizes_n*rep) .* proj2mm;


target_sizes = [5,10,15,20,25,30] ./ pixellength ./ proj2tablet;
target_sizes = repmat(target_sizes,1,dists_n*rep);
target_sizes = target_sizes';
target_sizes = target_sizes(:)';

seeds = [randperm(size(distances,2)), randperm(size(distances,2))];
randdists = distances(seeds);
randdists = randdists(:);
randsizes = target_sizes(seeds);
randsizes = randsizes(:);
%%
StoSimData = NaN(length(randsizes),1);
for i = 1:length(randsizes)
    if rem(i,2)
        theta = -pi/12 + (pi/6) * rand(1);
        rho = randdists(i)-UniRandRadius + UniRandRadius * 2 * rand(1);
        [offset(1),offset(2)] = pol2cart(theta,rho);
        StoSimData(i,1:2) = - offset;
    else
        theta = -pi/12 + (pi/6) * rand(1);
        rho = randdists(i)-UniRandRadius + UniRandRadius * 2 * rand(1);
        [offset(1),offset(2)] = pol2cart(theta,rho);
        StoSimData(i,1:2) = offset;
    end
    StoSimData(i,3) = randsizes(i); % size of the small and solid target
    
    optS_i = findOptSpeed8(norm(offset),randsizes(i),totalTime,speedRange,SAparams,Left);
    time = norm(offset)/optS_i;
    remainTime = totalTime - time;
    remainScore = (remainTime/totalTime)*maxScore;
    biasx = optS_i * UniLatParams(1) + UniLatParams(2);
    biasy = optS_i * UniLatParams(3) + UniLatParams(4);
    errorx = optS_i * UniLatParams(5) + UniLatParams(6);
    errory = optS_i * UniLatParams(7) + UniLatParams(8);
    reg_phit = compute_phit0(randsizes(i),errorx,errory, biasx, biasy);
    eGainReg = remainScore * reg_phit;
    alt_phit = eGainReg * totalTime / ((totalTime - time - penalty) * maxScore);
    % alt_phit is the alternative probability required for switching to be beneficial
    if alt_phit < 0.99
        fun = @(x) abs(compute_phit0(x,errorx,errory, biasx, biasy) - alt_phit);
        x0 = randsizes(i);
        StoSimData(i,4) = fminsearch(fun,x0); % size of the big and hollow target
        alt_scale = StoSimData(i,4)./randsizes(i); % the scale of magnification
    end
    
end