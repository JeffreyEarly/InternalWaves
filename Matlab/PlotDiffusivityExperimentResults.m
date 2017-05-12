file = '/Volumes/OceanTransfer/DiffusivityExperiment_2017-05-12T142805_128x16x33.nc';

t = ncread(file, 't');

Nx = length(ncread(file, 'x'));
Ny = length(ncread(file, 'y'));
Nz = length(ncread(file, 'z'));

nFloatLevels = ncreadatt(file, '/', 'nFloatLevels');
N0 = ncreadatt(file, '/', 'N0');
rho0 = 1025;
dz_drho = 9.81/(N0*N0*rho0);

x = ncread(file, 'x-position');
y = ncread(file, 'y-position');
z = ncread(file, 'z-position');
rho = ncread(file, 'density');

% x = ncread(file, 'x-position-diffusive');
% y = ncread(file, 'y-position-diffusive');
% z = ncread(file, 'z-position-diffusive');
% rho = ncread(file, 'density-diffusive');

% x = ncread(file, 'x-position-drifter');
% y = ncread(file, 'y-position-drifter');
% z = ncread(file, 'z-position-drifter');
% rho = ncread(file, 'density-drifter');

floatsPerLevel = size(x,1)/nFloatLevels;

r2 = zeros(floatsPerLevel*(floatsPerLevel-1)/2,nFloatLevels);
kappa_r = zeros(floatsPerLevel*(floatsPerLevel-1)/2,nFloatLevels);
D2z = zeros(length(t),nFloatLevels);
thelabels = cell(nFloatLevels,1);
for zLevel=1:nFloatLevels
    zLevelIndices = (zLevel-1)*floatsPerLevel + (1:floatsPerLevel);
    
    x_float = x(zLevelIndices,:)';
    y_float = y(zLevelIndices,:)';
    z_float = z(zLevelIndices,:)';
    rho_float = rho(zLevelIndices,:)';
    dz = (dz_drho * (rho(zLevelIndices,:)-mean(rho(zLevelIndices,:),1)))';
    
    thelabels{zLevel} = sprintf('%d meters',round(mean(z_float(1,:))));
    
   
    D2z(:,zLevel) = mean(dz.*dz,2);
    
    [r2(:,zLevel), kappa_r(:,zLevel)] = RelativeDiffusivity(t,x_float,y_float,1:round(length(t)/1.0));
    
end


figure
plot(t,D2z)
legend(thelabels,'Location', 'northwest')

figure
for zLevel=1:nFloatLevels
    x = sqrt(r2(:,zLevel));
    y = kappa_r(:,zLevel);
    
    [N,edges,bin] = histcounts(x,15);
    yMean = zeros(max(bin),1);
    yStdErr = zeros(max(bin),1);
    for i=1:max(bin)
        yMean(i) = mean(y(bin==i));
        yStdErr(i) = std(y(bin==i))/sqrt(sum(bin==i));
    end
    xMean = ((edges(1:end-1)+edges(2:end))/2)';

    errorbar(xMean, yMean,2*yStdErr), hold on
end
legend(thelabels,'Location', 'northwest')

