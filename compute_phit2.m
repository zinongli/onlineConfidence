function phit2 = compute_phit2(radius, sdx, sdy)
% radius = 45;    % target/penalty radius in pixels
% dsd = 1;       % resolution of SD in pixels
% SDs = 10:dsd:70;
% dx = .1;         % resolution of shifts in pixels

ninterval = 10000;
dy = radius/ninterval;
yList = ((1:ninterval) - .5)*dy;
xBound = sqrt(radius.^2 - yList.^2);
yGauss = (1/(sqrt(2*pi)*sdy))*exp(-yList.^2/(2*sdy^2));
phit2 = dy*sum(yGauss.*(erf((xBound)/(sdx*sqrt(2)))-erf((-xBound)/(sdx*sqrt(2)))));
end
