function pdf = mlPDFfittingLeft(x,a,b,c,d)
load('yLeftSampleSpeed,mat','leftward')
mu = a .* leftward + b;
sigma = c .* leftward + d;
pdf = normpdf(x,mu,sigma);
end