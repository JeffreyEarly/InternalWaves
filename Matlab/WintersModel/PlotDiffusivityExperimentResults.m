load('particles_LIN.mat');

tIndices = 1:length(t);

diffusivityMethod = 'powspec';

n = floatsPerLevel;
nLevels = size(x,2)/floatsPerLevel;

r2_r = cell(nLevels,1);
kappa_r = cell(nLevels,1);
kappa_r_corr = cell(nLevels,1);
r2_a = cell(nLevels,1);
kappa_a = cell(nLevels,1);
theZlabels = cell(nLevels,1);
thelabels = cell(nLevels,1);
for zLevel=1:nLevels
    zLevelIndices = (zLevel-1)*floatsPerLevel + (1:floatsPerLevel);
    zLevelIndices = reshape(zLevelIndices,round(sqrt(floatsPerLevel)),[]);
%     zLevelIndices = zLevelIndices(1:3:end,1:3:end);
    zLevelIndices = reshape(zLevelIndices,1,[]);
    
    x_float = x(tIndices,zLevelIndices);
    y_float = y(tIndices,zLevelIndices);
    z_float = z(tIndices,zLevelIndices);
    
    thelabels{zLevel} = sprintf('%d meters',round(mean(z_float(1,:))));
    [r2_r{zLevel}, kappa_r{zLevel}, kappa_r_corr{zLevel}] = PairwiseRelativeDiffusivity(t(tIndices),x_float,y_float,diffusivityMethod);
    [r2_a{zLevel}, kappa_a{zLevel}] = SingleParticleDiffusivity(t(tIndices),x_float,y_float,diffusivityMethod);
end

dx = abs(x(1,2)-x(1,1));
theBins = dx*(1:sqrt(floatsPerLevel)) + dx/2;
theBins(1) = 0;
figure
subplot(2,1,1)
for zLevel=1:nLevels
    a = sqrt(r2_r{zLevel});
    b = kappa_r{zLevel};
    
    [N,edges,bin] = histcounts(a,theBins);
    yMean = zeros(length(edges)-1,1);
    yStdErr = zeros(length(edges)-1,1);
    for i=1:max(bin)
        yMean(i) = mean(b(bin==i));
        yStdErr(i) = std(b(bin==i))/sqrt(sum(bin==i));
    end
    xMean = ((edges(1:end-1)+edges(2:end))/2)';

    errorbar(xMean, yMean,2*yStdErr), hold on
end

legend(thelabels,'Location', 'northwest')
xlabel('distance (m)')
ylabel('diffusivity (m^2/s)')
title(sprintf('horizontal *relative* diffusivity at different depths (%s method)',diffusivityMethod))

subplot(2,1,2)
for zLevel=1:nLevels
    a = sqrt(r2_a{zLevel});
    b = kappa_a{zLevel};
    
    [N,edges,bin] = histcounts(a,theBins);
    yMean = zeros(length(edges)-1,1);
    yStdErr = zeros(length(edges)-1,1);
    for i=1:max(bin)
        yMean(i) = mean(b(bin==i));
        yStdErr(i) = std(b(bin==i))/sqrt(sum(bin==i));
    end
    xMean = ((edges(1:end-1)+edges(2:end))/2)';

    errorbar(xMean, yMean,2*yStdErr), hold on
end

legend(thelabels,'Location', 'northeast')
xlabel('distance (m)')
ylabel('diffusivity (m^2/s)')
title(sprintf('horizontal *absolute* diffusivity at different depths (%s method)',diffusivityMethod))

