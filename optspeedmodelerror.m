
% SAparams = [xParams(2:3,1:2)';yParams(2:3,1:2)';xParams(2:3,3:4)';yParams(2:3,3:4)']
speedRange = [100 2000];
optSpeed8 = NaN(1,360);


for i = 1:360
    Left = copy(i,8)>1000;
    optSpeed8(i) = findOptSpeed8(copy(i,10),copy(i,15),copy(i,3),speedRange,SAparams,Left);
    i
end

%%
% sss = @(sigma,m,o) log(1/sqrt(2*pi)) * length(m) - (-log(sigma) * length(m) + sum(-((m'-o).^2)./(2 * sigma.^2)));
% sfun = @(sigma) sss(sigma,maxSpeed,optSpeed8);
% sigma_fit = bads(sfun,0.001,0.001,300);
% sss(sigma_fit,maxSpeed,optSpeed8);

UB = [5,500,5,500];
LB = [-5,-500,0.0001,1];
theta0 = [0,0,0.0001,1];
f = @(theta,maxS,optS) -(log(1/sqrt(2*pi)) * length(optS) + sum(-log(theta(3).*optS + theta(4)) - (((maxS - (theta(1).*optS + theta(2))).^2)./ (2.*(theta(3).*optS+theta(4)).^2))));
fun = @(theta) f(theta,maxSpeed,optSpeed8');
theta = bads(fun,theta0,LB,UB)

theta0 = [0,0,30];
UB = [5,500,500];
LB = [-5,-500,1];
f2 = @(theta,maxS,optS) -(log(1/sqrt(2*pi)) * length(optS) + sum(-log(theta(3)) - (((maxS - (theta(1).*optS + theta(2))).^2)./ (2.*(theta(3)).^2))));
fun2 = @(theta) f2(theta,maxSpeed,optSpeed8');
fixSig = bads(fun2,theta0,LB,UB)

fun2(fixSig)
fun(theta)

indices = 1:1600;
plot(optSpeed8,maxSpeed,'o')
hold on
plot(indices,(theta(1) .* indices) + theta(2),'--r')
plot(indices,(theta(1) .* indices) + theta(2) + ((theta(3) .* indices) + theta(4)),'--r')
plot(indices,(theta(1) .* indices) + theta(2) - ((theta(3) .* indices) + theta(4)),'--r')
plot(indices,(theta(1) .* indices) + theta(2) + 2*((theta(3) .* indices) + theta(4)),'--r')
plot(indices,(theta(1) .* indices) + theta(2) - 2*((theta(3) .* indices) + theta(4)),'--r')
hold off
xlabel('Speed mm/s')
ylabel('Error mm')
%%
pb = makedist('Normal');
qqplot(maxSpeed - optSpeed8',pb)
%%
[a,b,c] = swtest(maxSpeed - optSpeed8')
%%
set(groot,'defaultAxesFontSize',18)
figure
plot(optSpeed8,maxSpeed,'o')
hold on
plot(300:1400,300:1400,'--')
hold off
xlabel('Optimal Max Speed (mm/s)','FontSize',18)
ylabel('Actual Max Speed (mm/s)','FontSize',18)
% if participant == 1
%     title('JX Maximum Speed Data vs Prediction','FontSize',18)
% elseif participant == 0
%     title('ZL Maximum Speed Data vs Prediction','FontSize',18)
% end
%%

rSquared = sum((optSpeed8-mean(maxSpeed)).^2) / sum((maxSpeed-mean(maxSpeed)).^2)

%%
findOptSpeed8(80,10,3,speedRange,SAparams,Left)