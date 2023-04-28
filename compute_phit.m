function phit = compute_phit(radius, sd, percent_dist)
% radius = 45;    % target/penalty radius in pixels
% dsd = 1;       % resolution of SD in pixels
% SDs = 10:dsd:70;
% dx = .1;         % resolution of shifts in pixels
maxshift = 400;
aims = maxshift .* (1-percent_dist);
% numsd = length(SDs);
numaim = length(aims);
% phit = zeros(numsd,numaim);
phit = NaN(1,numaim);

r2 = radius^2;
ninterval = 1000;
dy = radius/ninterval;
midy = ((1:ninterval) - .5)*dy;
midy2 = midy.^2;
maxx = sqrt(r2 - midy.^2);

c1 = 1/(sqrt(2*pi)*sd);
c2 = 2*sd^2;
c3 = sd*sqrt(2);
yg = c1*exp(-midy2/c2);
for x = 1:numaim
    aim = aims(x);
    phit(x) = dy*sum(yg.*(erf((maxx-aim)/c3)-erf((-maxx-aim)/c3)));
end

