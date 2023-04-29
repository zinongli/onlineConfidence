function pdf = xMLEfittingLeft(x,a,b,c,d)
load('xLeftSampleSpeed.mat','leftwardSample')
mu = a .* leftwardSample + b;
sigma = c .* leftwardSample + d;
pdf = normpdf(x,mu,sigma);
end