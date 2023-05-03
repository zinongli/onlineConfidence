%%

xTotal = copy(:,23);
xRight = copy(copy(:,8)<1000,23);
xLeft = copy(copy(:,8)>1000,23);


theta0 = [0,0,0.0001,5];
UB = [1,10,1,50];
LB = [-1,-10,0.0001,1];

f = @(theta,speed,error) log(1/sqrt(2*pi)) * length(speed) - sum(-log(theta(3).*speed + theta(4)) - (((error - (theta(1).*speed + theta(2))).^2)./ (2.*(theta(3).*speed+theta(4)).^2)));


xParams = NaN(3,4);
fun = @(theta) f(theta,maxSpeed,xTotal);
xParams(1,:) = bads(fun,theta0,LB,UB);
fun = @(theta) f(theta,rightward,xRight);
xParams(2,:) = bads(fun,theta0,LB,UB);
fun = @(theta) f(theta,leftward,xLeft);
xParams(3,:) = bads(fun,theta0,LB,UB);


yTotal = copy(:,29);
yRight = copy(copy(:,8)<1000,29);
yLeft = copy(copy(:,8)>1000,29);


theta0 = [0,0,0.0001,5];
UB = [1,10,1,50];
LB = [-1,-10,0.0001,1];

f = @(theta,speed,error) log(1/sqrt(2*pi)) * length(speed) - sum(-log(theta(3).*speed + theta(4)) - (((error - (theta(1).*speed + theta(2))).^2)./ (2.*(theta(3).*speed+theta(4)).^2)));


yParams = NaN(3,4);
fun = @(theta) f(theta,maxSpeed,yTotal);
yParams(1,:) = bads(fun,theta0,LB,UB);
fun = @(theta) f(theta,rightward,yRight);
yParams(2,:) = bads(fun,theta0,LB,UB);
fun = @(theta) f(theta,leftward,yLeft);
yParams(3,:) = bads(fun,theta0,LB,UB);
%%
xTotalLLH = log(1/sqrt(2*pi)) * length(maxSpeed) + sum(-log(xParams(1,3).*maxSpeed + xParams(1,4)) - (((xTotal - (xParams(1,1).*maxSpeed + xParams(1,2))).^2)./ (2.*(xParams(1,3).*maxSpeed+xParams(1,4)))));
xRightLLH = log(1/sqrt(2*pi)) * length(rightward) + sum(-log(xParams(2,3).*rightward + xParams(2,4)) - (((xRight - (xParams(2,1).*rightward + xParams(2,2))).^2)./ (2.*(xParams(2,3).*rightward+xParams(2,4)))));
xLeftLLH = log(1/sqrt(2*pi)) * length(leftward) + sum(-log(xParams(3,3).*leftward + xParams(3,4)) - (((xLeft - (xParams(3,1).*leftward + xParams(3,2))).^2)./ (2.*(xParams(3,3).*leftward+xParams(3,4)))));
xLateralLLH = xRightLLH + xLeftLLH;

xaicTotal = 2 * 4 - 2 * xTotalLLH
xaicIndividual = 2 * 8 - 2 * xLateralLLH
xbicTotal = 4 * log(length(maxSpeed)) - 2 * xTotalLLH
xbicIndividual = 8 * log(length(maxSpeed)) - 2 * xLateralLLH

yTotalLLH = log(1/sqrt(2*pi)) * length(maxSpeed) + sum(-log(yParams(1,3).*maxSpeed + yParams(1,4)) - (((yTotal - (yParams(1,1).*maxSpeed + yParams(1,2))).^2)./ (2.*(yParams(1,3).*maxSpeed+yParams(1,4)))));
yRightLLH = log(1/sqrt(2*pi)) * length(rightward) + sum(-log(yParams(2,3).*rightward + yParams(2,4)) - (((yRight - (yParams(2,1).*rightward + yParams(2,2))).^2)./ (2.*(yParams(2,3).*rightward+yParams(2,4)))));
yLeftLLH = log(1/sqrt(2*pi)) * length(leftward) + sum(-log(yParams(3,3).*leftward + yParams(3,4)) - (((yLeft - (yParams(3,1).*leftward + yParams(3,2))).^2)./ (2.*(yParams(3,3).*leftward+yParams(3,4)))));
yLateralLLH = yRightLLH + yLeftLLH;

yaicTotal = 2 * 4 - 2 * yTotalLLH
yaicIndividual = 2 * 8 - 2 * yLateralLLH
ybicTotal = 4 * log(length(maxSpeed)) - 2 * yTotalLLH
ybicIndividual = 8 * log(length(maxSpeed)) - 2 * yLateralLLH
%%
figure(1)
indices = 1:1600;
subplot(1,3,1)

plot(maxSpeed,xTotal,'o')
hold on
plot(indices,(xParams(1,1) .* indices) + xParams(1,2),'--r')
plot(indices,(xParams(1,1) .* indices) + xParams(1,2) + (xParams(1,3) * indices + xParams(1,4)),'--r')
plot(indices,(xParams(1,1) .* indices) + xParams(1,2) - (xParams(1,3) * indices + xParams(1,4)),'--r')
plot(indices,(xParams(1,1) .* indices) + xParams(1,2) + 2*(xParams(1,3) * indices + xParams(1,4)),'--r')
plot(indices,(xParams(1,1) .* indices) + xParams(1,2) - 2*(xParams(1,3) * indices + xParams(1,4)),'--r')
hold off
xlabel('Speed mm/s')
ylabel('Error mm')

subplot(1,3,2)
plot(rightward,xRight,'o')
hold on
plot(indices,(xParams(2,1) .* indices) + xParams(2,2),'--r')
plot(indices,(xParams(2,1) .* indices) + xParams(2,2) + (xParams(2,3) * indices + xParams(2,4)),'--r')
plot(indices,(xParams(2,1) .* indices) + xParams(2,2) - (xParams(2,3) * indices + xParams(2,4)),'--r')
plot(indices,(xParams(2,1) .* indices) + xParams(2,2) + 2*(xParams(2,3) * indices + xParams(2,4)),'--r')
plot(indices,(xParams(2,1) .* indices) + xParams(2,2) - 2*(xParams(2,3) * indices + xParams(2,4)),'--r')
hold off
xlabel('Speed mm/s')
ylabel('Error mm')

subplot(1,3,3)
plot(leftward,xLeft,'o')
hold on
plot(indices,(xParams(3,1) .* indices) + xParams(3,2),'--r')
plot(indices,(xParams(3,1) .* indices) + xParams(3,2) + (xParams(3,3) * indices + xParams(3,4)),'--r')
plot(indices,(xParams(3,1) .* indices) + xParams(3,2) - (xParams(3,3) * indices + xParams(3,4)),'--r')
plot(indices,(xParams(3,1) .* indices) + xParams(3,2) + 2*(xParams(3,3) * indices + xParams(3,4)),'--r')
plot(indices,(xParams(3,1) .* indices) + xParams(3,2) - 2*(xParams(3,3) * indices + xParams(3,4)),'--r')
hold off
xlabel('Speed mm/s')
ylabel('Error mm')

sgtitle('Error Along the Reach Direction','FontSize',18)


figure(2)
indices = 1:1600;
subplot(1,3,1)

plot(maxSpeed,yTotal,'o')
hold on
plot(indices,(yParams(1,1) .* indices) + yParams(1,2),'--r')
plot(indices,(yParams(1,1) .* indices) + yParams(1,2) + (yParams(1,3) * indices + yParams(1,4)),'--r')
plot(indices,(yParams(1,1) .* indices) + yParams(1,2) - (yParams(1,3) * indices + yParams(1,4)),'--r')
plot(indices,(yParams(1,1) .* indices) + yParams(1,2) + 2*(yParams(1,3) * indices + yParams(1,4)),'--r')
plot(indices,(yParams(1,1) .* indices) + yParams(1,2) - 2*(yParams(1,3) * indices + yParams(1,4)),'--r')
hold off
xlabel('Speed mm/s')
ylabel('Error mm')

subplot(1,3,2)
plot(rightward,yRight,'o')
hold on
plot(indices,(yParams(2,1) .* indices) + yParams(2,2),'--r')
plot(indices,(yParams(2,1) .* indices) + yParams(2,2) + (yParams(2,3) * indices + yParams(2,4)),'--r')
plot(indices,(yParams(2,1) .* indices) + yParams(2,2) - (yParams(2,3) * indices + yParams(2,4)),'--r')
plot(indices,(yParams(2,1) .* indices) + yParams(2,2) + 2*(yParams(2,3) * indices + yParams(2,4)),'--r')
plot(indices,(yParams(2,1) .* indices) + yParams(2,2) - 2*(yParams(2,3) * indices + yParams(2,4)),'--r')
hold off
xlabel('Speed mm/s')
ylabel('Error mm')

subplot(1,3,3)
plot(leftward,yLeft,'o')
hold on
plot(indices,(yParams(3,1) .* indices) + yParams(3,2),'--r')
plot(indices,(yParams(3,1) .* indices) + yParams(3,2) + (yParams(3,3) * indices + yParams(3,4)),'--r')
plot(indices,(yParams(3,1) .* indices) + yParams(3,2) - (yParams(3,3) * indices + yParams(3,4)),'--r')
plot(indices,(yParams(3,1) .* indices) + yParams(3,2) + 2*(yParams(3,3) * indices + yParams(3,4)),'--r')
plot(indices,(yParams(3,1) .* indices) + yParams(3,2) - 2*(yParams(3,3) * indices + yParams(3,4)),'--r')
hold off
xlabel('Speed mm/s')
ylabel('Error mm')

sgtitle('Error Orthogonal the Reach Direction','FontSize',18)