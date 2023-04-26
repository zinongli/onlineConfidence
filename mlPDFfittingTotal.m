function pdf = mlPDFfittingTotal(x,a,b,c,d)
load('yTotalSampleSpeed,mat','maxSpeed')
mu = a .* maxSpeed + b;
sigma = c .* maxSpeed + d;
pdf = normpdf(x,mu,sigma);
end