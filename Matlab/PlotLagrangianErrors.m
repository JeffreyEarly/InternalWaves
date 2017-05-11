
file = '/Volumes/OceanTransfer/LagrangianErrorExperiment_2017-05-10T144505_128x16x33.nc';

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
xSplineInterpDiffusivity = zeros(length(t),nFloatLevels);
ySplineInterpDiffusivity = zeros(length(t),nFloatLevels);
zSplineInterpDiffusivity = zeros(length(t),nFloatLevels);
zIsoSplineInterpDiffusivity = zeros(length(t),nFloatLevels);
zIsoLinearInterpDiffusivity = zeros(length(t),nFloatLevels);
kappa_z = zeros(length(t),nFloatLevels);
for zLevel=1:nFloatLevels
    zLevelIndices = (zLevel-1)*floatsPerLevel + (1:floatsPerLevel);
    initialRho = mean(rho(zLevelIndices,1));
    
    % Distance from the initial isopycnal.
    dz(zLevelIndices,:) = dz_drho * (rho(zLevelIndices,:)-initialRho);
    dzLinear(zLevelIndices,:) = dz_drho * (rhoLinear(zLevelIndices,:)-initialRho);
    dzSpline(zLevelIndices,:) = dz_drho * (rhoSpline(zLevelIndices,:)-initialRho);
    
    % Diffusivity associated with the interpolation method, at each level
    zIsoLinearInterpDiffusivity(:,zLevel) = mean((dzLinear(zLevelIndices,:)-dz(zLevelIndices,:)).^2,1)'./t;
    
    xSplineInterpDiffusivity(:,zLevel) = mean((xSpline(zLevelIndices,:)-x(zLevelIndices,:)).^2,1)'./t;
    ySplineInterpDiffusivity(:,zLevel) = mean((ySpline(zLevelIndices,:)-y(zLevelIndices,:)).^2,1)'./t;
    zSplineInterpDiffusivity(:,zLevel) = mean((zSpline(zLevelIndices,:)-z(zLevelIndices,:)).^2,1)'./t;
    zIsoSplineInterpDiffusivity(:,zLevel) = mean((dzSpline(zLevelIndices,:)-dz(zLevelIndices,:)).^2,1)'./t;
    
    % Physical diffusivity at each level
    kappa_z(:,zLevel) = (mean(dz(zLevelIndices,:).^2,1)'./t)/2;
    
end

kappa_h_interp = (xSplineInterpDiffusivity + ySplineInterpDiffusivity)/4;
kappa_z_interp = zIsoSplineInterpDiffusivity/2;
kappa_z_interp_linear = zIsoLinearInterpDiffusivity/2;

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

figure
plot(t,kappa_z,'--'), ylog, hold on
plot(t,kappa_z_interp)
plot(t,kappa_z_interp_linear)
title(sprintf('%dx%dx%d',Nx,Ny,Nz))
xlabel('time (s)')
ylabel('diffusivity (m^2/s)')

return

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
