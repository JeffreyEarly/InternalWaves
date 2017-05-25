% 50km domain, full wave spectrum
file = '/Volumes/OceanTransfer/DiffusivityExperiment_2017-05-21T081340_128x128x129.nc';

% 500km domain, full wave spectrum
file = '/Volumes/OceanTransfer/DiffusivityExperiment_2017-05-23T224848_128x128x129.nc';

t = ncread(file, 't');

Nx = length(ncread(file, 'x'));
Ny = length(ncread(file, 'y'));
Nz = length(ncread(file, 'z'));

nLevels = ncreadatt(file, '/', 'nFloatLevels');
interpolationMethod = ncreadatt(file, '/', 'interpolation-method');
N0 = ncreadatt(file, '/', 'N0');
rho0 = 1025;
dz_drho = 9.81/(N0*N0*rho0);

x = ncread(file, 'x-position');
y = ncread(file, 'y-position');
z = ncread(file, 'z-position');
rho = ncread(file, 'density');

x = ncread(file, 'x-position-diffusive');
y = ncread(file, 'y-position-diffusive');
z = ncread(file, 'z-position-diffusive');
rho = ncread(file, 'density-diffusive');

% x = ncread(file, 'x-position-drifter');
% y = ncread(file, 'y-position-drifter');
% z = ncread(file, 'z-position-drifter');
% rho = ncread(file, 'density-drifter');

floatsPerLevel = size(x,1)/nLevels;

tIndices = (1:(length(t)-1));

diffusivityMethod = 'endpoint';

n = floatsPerLevel;
r2 = cell(nLevels,1);
kappa_r = cell(nLevels,1);
kappa_r_corr = cell(nLevels,1);
kappa_z = zeros(nLevels,1);
D2z = zeros(length(tIndices),nLevels);
theZlabels = cell(nLevels,1);
thelabels = cell(nLevels,1);
for zLevel=1:nLevels
    zLevelIndices = (zLevel-1)*floatsPerLevel + (1:floatsPerLevel);
    zLevelIndices = reshape(zLevelIndices,round(sqrt(floatsPerLevel)),[]);
    zLevelIndices = zLevelIndices(1:3:end,1:3:end);
    zLevelIndices = reshape(zLevelIndices,1,[]);
    
    x_float = x(zLevelIndices,tIndices)';
    y_float = y(zLevelIndices,tIndices)';
    z_float = z(zLevelIndices,tIndices)';
    rho_float = rho(zLevelIndices,tIndices)';
    dz = (dz_drho * (rho(zLevelIndices,tIndices)-mean(rho(zLevelIndices,tIndices),1)))';
      
    D2z(:,zLevel) = mean(dz.*dz,2);
    [p,~,mu]=polyfit(t(tIndices),D2z(:,zLevel),1);
    kappa_z(zLevel) = (p(1)/mu(2))/2;
    theZlabels{zLevel} = sprintf('%d meters (kappa=%.2g)',round(mean(z_float(1,:))),kappa_z(zLevel));
    
    thelabels{zLevel} = sprintf('%d meters',round(mean(z_float(1,:))));
    [r2{zLevel}, kappa_r{zLevel}, kappa_r_corr{zLevel}] = RelativeDiffusivity(t(tIndices),x_float,y_float,diffusivityMethod);
    
%     [a, b] = PatchDiffusivity(t(tIndices),x_float,y_float,1,sqrt(floatsPerLevel));
%     [r2(1:length(a),zLevel), kappa_r(1:length(a),zLevel)] = PatchDiffusivity(t(tIndices),x_float,y_float,1,sqrt(floatsPerLevel));
%     
%     [ACx, DOFx] = Autocorrelation(x_float, length(dt)-1);
%     [ACy, DOFy] = Autocorrelation(y_float, length(dt)-1);
end

% u0 = real(velocities(:,1));
% for n=1:16
%     un = real(velocities(:,n));
%     r(n) = mean((u0 - mean(u0)).*(un - mean(un)))/(std(u0,1)*std(un,1));
% end

figure
plot(t(tIndices),D2z)
legend(theZlabels,'Location', 'northwest')

xGrid = ncread(file,'x');
theBins = xGrid(1:sqrt(floatsPerLevel)+1) + (xGrid(2)-xGrid(1))/2;
theBins(1) = 0;
figure
subplot(2,1,1)
for zLevel=1:nLevels
    a = sqrt(r2{zLevel});
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
% vlines(4*(xGrid(2)-xGrid(1)),'k2--')
legend(thelabels,'Location', 'northwest')
xlabel('distance (m)')
ylabel('diffusivity (m^2/s)')
title(sprintf('horizontal diffusivity at different depths (%s method)',diffusivityMethod))

return

subplot(2,1,2)
for zLevel=1:nFloatLevels
    a = sqrt(r2{zLevel});
    b = kappa_r_corr{zLevel};
    
    [N,edges,bin] = histcounts(a,theBins);
    yMean = zeros(max(bin),1);
    yStdErr = zeros(max(bin),1);
    for i=1:max(bin)
        yMean(i) = mean(b(bin==i));
        yStdErr(i) = std(b(bin==i))/sqrt(sum(bin==i));
    end
    xMean = ((edges(1:end-1)+edges(2:end))/2)';

    errorbar(xMean, yMean,2*yStdErr), hold on
end
% vlines(4*(xGrid(2)-xGrid(1)),'k2--')
legend(thelabels,'Location', 'northeast')
xlabel('distance (m)')
ylabel('correlation coefficient')
title('velocity correlation vs distance')

