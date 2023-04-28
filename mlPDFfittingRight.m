function pdf = mlPDFfittingRight(x,a,b,c,d)
load('yRightSampleSpeed,mat','rightwardSample')
mu = a .* rightwardSample + b;
sigma = c .* rightwardSample + d;
pdf = normpdf(x,mu,sigma);
end