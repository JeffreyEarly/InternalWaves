
file = '/Users/jearly/Desktop/LagrangianErrorExperiment_2017-05-05T113306_128x16x33.nc';
file = '/Users/jearly/Desktop/LagrangianErrorExperiment_2017-05-05T130219_128x16x17.nc';
file = '/Users/jearly/Desktop/LagrangianErrorExperiment_2017-05-05T131106_128x16x65.nc';
file = '/Volumes/OceanTransfer/LagrangianErrorExperiment_2017-05-08T143435_128x16x17.nc';
file = '/Volumes/OceanTransfer/LagrangianErrorExperiment_2017-05-08T151314_128x16x17.nc';
file = '/Users/jearly/Desktop/LagrangianErrorExperiment_2017-05-08T193423_16x16x17.nc';
file = '/Users/jearly/Desktop/LagrangianErrorExperiment_2017-05-08T195011_16x16x17.nc';

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

xLinear = ncread(file, 'x-position-linear');
yLinear = ncread(file, 'y-position-linear');
zLinear = ncread(file, 'z-position-linear');
rhoLinear = ncread(file, 'density-linear');

xSpline = ncread(file, 'x-position-spline');
ySpline = ncread(file, 'y-position-spline');
zSpline = ncread(file, 'z-position-spline');
rhoSpline = ncread(file, 'density-spline');

floatsPerLevel = size(x,1)/nFloatLevels;
dz = zeros(size(z));
dzLinear = zeros(size(z));
dzSpline = zeros(size(z));
for zLevel=1:nFloatLevels
    zLevelIndices = (zLevel-1)*floatsPerLevel + (1:floatsPerLevel);
    initialRho = mean(rho(zLevelIndices,1));
    dz(zLevelIndices,:) = dz_drho * (rho(zLevelIndices,:)-initialRho);
    dzLinear(zLevelIndices,:) = dz_drho * (rhoLinear(zLevelIndices,:)-initialRho);
    dzSpline(zLevelIndices,:) = dz_drho * (rhoSpline(zLevelIndices,:)-initialRho);
end

% z_plane  = (z(101:200,:))';
% zeta_plane = (zeta(101:200,:))';
% 
% z2 = MeanSquareSeparation( z_plane, zeros(size(z_plane)) );
% iso2 = MeanSquareSeparation( z_plane-zeta_plane, zeros(size(z_plane)) );

% xSplineDiffusion = mean((xSpline(:,end)-x(:,end)).^2)/t(end)
% ySplineDiffusion = mean((ySpline(:,end)-y(:,end)).^2)/t(end)
% zSplineDiffusion = mean((zSpline(:,end)-z(:,end)).^2)/t(end)

% drifterIndices = 101;
% figure
% plot( x(drifterIndices,:), y(drifterIndices,:) ), hold on
% plot( xSpline(drifterIndices,:), ySpline(drifterIndices,:) )


xSplineDiffusion = mean((xSpline-x).^2,1)'./t;
ySplineDiffusion = mean((ySpline-y).^2,1)'./t;
zSplineDiffusion = mean((zSpline(101:300,:)-z(101:300,:)).^2,1)'./t;
zIsoSplineDiffusion = mean((dz(101:300,:)-dzSpline(101:300,:)).^2,1)'./t;
zDiffusion = mean(dz(101:300,:).^2,1)'./t;

kappa_h = (xSplineDiffusion + xSplineDiffusion)/4;
kappa_z = zIsoSplineDiffusion/2;
kappa_z_abs = zDiffusion/2;

kappa_z0000 = mean((dzSpline(1:100,:)-dz(1:100,:)).^2,1)'./t;
kappa_z1250 = mean((dzSpline(101:200,:)-dz(101:200,:)).^2,1)'./t;
kappa_z2500 = mean((dzSpline(201:300,:)-dz(201:300,:)).^2,1)'./t;

kappa_z1250_abs = mean(dz(101:200,:).^2,1)'./t;
kappa_z2500_abs = mean(dz(201:300,:).^2,1)'./t;

figure
subplot(1,2,1)
plot(t,[kappa_h kappa_z kappa_z_abs]), ylog
title(sprintf('%dx%dx%d',Nx,Ny,Nz))
xlabel('time (s)')
ylabel('diffusivity (m^2/s)')
legend('horizontal diffusivity', 'vertical diffusivity', 'abs vertical diffusivity')

subplot(1,2,2)
plot(t, [kappa_z1250 kappa_z2500]), ylog, hold on
plot(t, [kappa_z1250_abs kappa_z2500_abs],'--')
xlabel('time (s)')
ylabel('diffusivity (m^2/s)')
legend('1250 m', '2500 m', '1250 m (abs)', '2500 m (abs)')
