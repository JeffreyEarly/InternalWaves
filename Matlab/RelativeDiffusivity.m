function [r2, kappa_r,correlation] = RelativeDiffusivity(t,x,y,timeIndices, method)
nDrifters = size(x,2);
nReps = nDrifters*(nDrifters-1)/2;
r2 = zeros(nReps,1);
kappa_r = zeros(nReps,1);
correlation = zeros(nReps,1);
velocities = zeros(length(timeIndices)-1,nReps);

stride = 1;
iteration = 1;
for iDrifter=1:stride:nDrifters
    for jDrifter = (iDrifter+1):stride:nDrifters
        dq = x(timeIndices,iDrifter) - x(timeIndices,jDrifter);
        dr = y(timeIndices,iDrifter) - y(timeIndices,jDrifter);
        
        % squared-separation distance as a function of time.
        r2t = (dq).^2 + (dr).^2;
               
        if strcmp(method,'endpoint')
            kappa_r(iteration) = (r2t(max(timeIndices)) - r2t(min(timeIndices)))/(4*max(t(timeIndices))-min(t(timeIndices)));
        elseif strcmp(method,'slope')
            [p,~,mu]=polyfit(t(timeIndices),r2t(timeIndices),1);
            kappa_r(iteration) = (p(1)/mu(2))/4;
        elseif strcmp(method,'powspec')
            velocities(:,iteration) = diff(dq+sqrt(-1)*dr)/(t(2)-t(1));
        end
        
        % mean-squared-separation distance over the requested time-interval
        r2(iteration) = mean(r2t(timeIndices),1);
        
%         [~,~,q,r] = CenterOfMass( x(:,[iDrifter jDrifter]), y(:,[iDrifter jDrifter]) );
        dt = t(2)-t(1);
        u = diff(x(:,[iDrifter jDrifter]))/dt;
        v = diff(y(:,[iDrifter jDrifter]))/dt;
        u = u - mean(u,1);
        v = v - mean(v,1);
        f = @(u1,u2) mean(u1.*u2)/(std(u1,1)*std(u2,1));
        correlation(iteration) = (f(u(:,1),u(:,2)) + f(v(:,1),v(:,2)))/2;
        
        iteration = iteration+1;
    end
end

if strcmp(method,'powspec')
    averaging_bandwidth = 1;
    taper_bandwidth = 0;
    kappa_r = DiffusivityFromZeroFrequency(t(2)-t(1),velocities,averaging_bandwidth,taper_bandwidth)';
end

end