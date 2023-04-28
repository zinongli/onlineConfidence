function pdf = mlPDFfittingLeft(x,a,b,c,d)
load('yLeftSampleSpeed,mat','leftwardSample')
mu = a .* leftwardSample + b;
sigma = c .* leftwardSample + d;
pdf = normpdf(x,mu,sigma);
end