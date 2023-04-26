function pdf = mlPDFfittingRight(x,a,b,c,d)
load('yRightSampleSpeed,mat','rightward')
mu = a .* rightward + b;
sigma = c .* rightward + d;
pdf = normpdf(x,mu,sigma);
end