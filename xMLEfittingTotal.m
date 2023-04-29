function pdf = xMLEfittingTotal(x,a,b,c,d)
load('xTotalSampleSpeed.mat','maxSpeedSample')
mu = a .* maxSpeedSample + b;
sigma = c .* maxSpeedSample + d;
pdf = normpdf(x,mu,sigma);
end