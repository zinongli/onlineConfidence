function pdf = mlPDFfittingTotal(x,a,b,c,d)
load('yTotalSampleSpeed,mat','maxSpeedSample')
mu = a .* maxSpeedSample + b;
sigma = c .* maxSpeedSample + d;
pdf = normpdf(x,mu,sigma);
end