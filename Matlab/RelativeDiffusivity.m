function [r2, kappa_r] = RelativeDiffusivity(t,x,y,timeIndices)
nDrifters = size(x,2);
r2 = zeros(nDrifters,1);
kappa_r = zeros(nDrifters,1);
stride = 1;

iteration = 1;
for iDrifter=1:stride:nDrifters
    for jDrifter = (iDrifter+1):stride:nDrifters
        dq = x(timeIndices,iDrifter) - x(timeIndices,jDrifter);
        dr = y(timeIndices,iDrifter) - y(timeIndices,jDrifter);
        
        r2(iteration) = mean((dq).^2 + (dr).^2,1);
        kappa_r(iteration) = ((dq(end))^2 + (dr(end)).^2 - (dq(1))^2 - (dr(1)).^2)/(4*max(t(timeIndices))-min(t(timeIndices)));
        
        iteration = iteration+1;
    end
end

end