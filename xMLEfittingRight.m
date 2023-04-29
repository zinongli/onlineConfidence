function pdf = xMLEfittingRight(x,a,b,c,d)
load('xRightSampleSpeed.mat','rightwardSample')
mu = a .* rightwardSample + b;
sigma = c .* rightwardSample + d;
pdf = normpdf(x,mu,sigma);
end