% plane wave
file = '/Volumes/OceanTransfer/LagrangianErrorExperiment_2018-02-20T143117_64x64x65.nc';

% plane wave, set to external
% file = '/Volumes/OceanTransfer/LagrangianErrorExperiment_2018-02-20T142341_64x64x65.nc';

% GM spectrum
% file = '/Volumes/OceanTransfer/LagrangianErrorExperiment_2018-02-20T140104_64x64x65.nc';

% 2 plane waves
file = '/Volumes/OceanTransfer/LagrangianErrorExperiment_2018-02-20T143439_64x64x65.nc';

% single, high mode plane wave.
file = '/Volumes/OceanTransfer/LagrangianErrorExperiment_2018-02-20T145027_64x64x65.nc';

% two plane waves
% file = '/Volumes/OceanTransfer/LagrangianErrorExperiment_2018-02-20T145216_64x64x65.nc';

t = ncread(file, 't');

Nx = length(ncread(file, 'x'));
Ny = length(ncread(file, 'y'));
Nz = length(ncread(file, 'z'));

nFloatLevels = ncreadatt(file, '/', 'nFloatLevels');
N0 = ncreadatt(file, '/', 'N0');
rho0 = 1025;
dz_drho = 9.81/(N0*N0*rho0);

x = ncread(file, 'x-position-exact')';
y = ncread(file, 'y-position-exact')';
z = ncread(file, 'z-position-exact')';
rho = ncread(file, 'density-exact')';

xLinear = ncread(file, 'x-position-linear')';
yLinear = ncread(file, 'y-position-linear')';
zLinear = ncread(file, 'z-position-linear')';
rhoLinear = ncread(file, 'density-linear')';

xSpline = ncread(file, 'x-position-spline')';
ySpline = ncread(file, 'y-position-spline')';
zSpline = ncread(file, 'z-position-spline')';
rhoSpline = ncread(file, 'density-spline')';

kappa_z_exact =  ((dz_drho*(rho(2:end)-rho(1))).^2)./t(2:end);
kappa_z_linear =  ((dz_drho*(rhoLinear(2:end)-rhoLinear(1))).^2)./t(2:end);
kappa_z_spline =  ((dz_drho*(rhoSpline(2:end)-rhoSpline(1))).^2)./t(2:end);

figure
plot(t(2:end),kappa_z_exact), ylog, hold on
plot(t(2:end),kappa_z_linear)
plot(t(2:end),kappa_z_spline)
legend('exact', 'linear', 'spline')
ylabel('diffusivity (m^2/s)')
xlabel('time (s)')
title('Numerical diffusivity in the vertical')