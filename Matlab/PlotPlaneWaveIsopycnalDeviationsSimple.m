file = '/Volumes/OceanTransfer/PlaneWaveIsopycnalExperiment_2017-06-08T131852_64x64x65.nc';

t = ncread(file, 't');

Nx = length(ncread(file, 'x'));
Ny = length(ncread(file, 'y'));
Nz = length(ncread(file, 'z'));

nFloatLevels = ncreadatt(file, '/', 'nFloatLevels');
N0 = ncreadatt(file, '/', 'N0');
rho0 = 1025;
dz_drho = 9.81/(N0*N0*rho0);

x = ncread(file, 'x-position-exact');
y = ncread(file, 'y-position-exact');
z = ncread(file, 'z-position-exact');
rho = ncread(file, 'density-exact');

floatsPerLevel = size(x,1)/nFloatLevels;
dz = zeros(size(z));
for zLevel=1:nFloatLevels
    zLevelIndices = (zLevel-1)*floatsPerLevel + (1:floatsPerLevel);
    initialRho = mean(rho(zLevelIndices,1));
    
    % Distance from the initial isopycnal.
    dz(zLevelIndices,:) = dz_drho * (rho(zLevelIndices,:)-initialRho);
end

figure, plot(t(1:end-1),dz(2,1:end-1))
