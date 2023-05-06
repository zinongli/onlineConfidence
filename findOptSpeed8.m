function opt_speed = findOptSpeed8(distance,tSize,totalTime,speedRange,SAparams,Left)

% SAparams: 
% |    R    |    L    |
% |--------------------
% |  X B S  |  X B S  |
% |  X B I  |  X B I  |
% |  Y B S  |  Y B S  |
% |  Y B I  |  Y B I  |
% |  X S S  |  X S S  |
% |  X S I  |  X S I  |
% |  Y S S  |  Y S S  |
% |  Y S I  |  Y S I  |

maxScore = 10;
S_0 = mean(speedRange);
params = SAparams(:,Left + 1); %therefore SAparams is right column | left column

f = @(speed,params,tSize,maxScore,distance,totalTime) -((totalTime - (distance/(speed*params(9)+params(10))))/totalTime)*maxScore * compute_phit8(tSize,speed*params(5)+params(6),speed*params(7)+params(8),speed*params(1)+params(2),speed*params(3)+params(4));
fun = @(speed) f(speed,params,tSize,maxScore,distance,totalTime);
opt_speed = bads(fun,S_0,speedRange(1),speedRange(2),speedRange(1),speedRange(2));

%     avgSpeed = speed * params(9) + params(10);
%     time = distance/avgSpeed;
%     remainTime = totalTime - time;
%     remainScore = (remainTime/totalTime)*maxScore;
%     
%     biasx = speed * params(1) + params(2);
%     biasy = speed * params(3) + params(4);
%     errorx = speed * params(5) + params(6);
%     errory = speed * params(7) + params(8);
%     
%     probHit1 = compute_phit0(tSize,errorx,errory, biasx, biasy);
%     E_Gain = remainScore * probHit1;

end