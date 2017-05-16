
% file = '/Volumes/OceanTransfer/LagrangianErrorExperiment_2017-05-09T224927_256x32x65.nc';
file = '/Volumes/OceanTransfer/LagrangianErrorExperiment_2017-05-10T144505_128x16x33.nc';
file = '/Volumes/OceanTransfer/LagrangianErrorExperiment_2017-05-11T103240_128x16x65.nc';

% This file zeros everything outside of the k=nyquist/4 band
file = '/Volumes/OceanTransfer/LagrangianErrorExperiment_2017-05-13T153340_64x64x65.nc';

% This file zeros everything outside of the k=nyquist/2 band
file = '/Volumes/OceanTransfer/LagrangianErrorExperiment_2017-05-14T165631_64x64x65.nc';

% This file zeros everything outside of the k=3*nyquist/4 band
file = '/Volumes/OceanTransfer/LagrangianErrorExperiment_2017-05-14T170830_64x64x65.nc';
% Cubic spline does really well, largely matching even at 3/4th of nyquist.

% This file zeros everything outside of the k=nyquist-2*dk band
file = '/Volumes/OceanTransfer/LagrangianErrorExperiment_2017-05-14T172441_64x64x65.nc';

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
dzTracer = zeros(size(z));
dzLinear = zeros(size(z));
dzSpline = zeros(size(z));
xSplineInterpDiffusivity = zeros(length(t),nFloatLevels);
ySplineInterpDiffusivity = zeros(length(t),nFloatLevels);
zSplineInterpDiffusivity = zeros(length(t),nFloatLevels);
zIsoSplineInterpDiffusivity = zeros(length(t),nFloatLevels);
xLinearInterpDiffusivity = zeros(length(t),nFloatLevels);
yLinearInterpDiffusivity = zeros(length(t),nFloatLevels);
zLinearInterpDiffusivity = zeros(length(t),nFloatLevels);
zIsoLinearInterpDiffusivity = zeros(length(t),nFloatLevels);
kappa_z = zeros(length(t),nFloatLevels);
kappa_z_centered = zeros(length(t),nFloatLevels);
kappa_h = zeros(length(t),nFloatLevels);
for zLevel=1:nFloatLevels
    zLevelIndices = (zLevel-1)*floatsPerLevel + (1:floatsPerLevel);
    initialRho = mean(rho(zLevelIndices,1));
    
    % Distance from the initial isopycnal.
    dz(zLevelIndices,:) = dz_drho * (rho(zLevelIndices,:)-initialRho);
    dzTracer(zLevelIndices,:) = dz_drho * (rho(zLevelIndices,:)-mean(rho(zLevelIndices,:),1));
    dzLinear(zLevelIndices,:) = dz_drho * (rhoLinear(zLevelIndices,:)-initialRho);
    dzSpline(zLevelIndices,:) = dz_drho * (rhoSpline(zLevelIndices,:)-initialRho);
    dzTracerSpline(zLevelIndices,:) = dz_drho * (rhoSpline(zLevelIndices,:)-mean(rhoSpline(zLevelIndices,:),1));
    dzTracerLinear(zLevelIndices,:) = dz_drho * (rhoLinear(zLevelIndices,:)-mean(rhoLinear(zLevelIndices,:),1));
    
    zD2Tracer(:,zLevel) = mean((dzTracer(zLevelIndices,:)).^2,1)';
    zD2TracerSpline(:,zLevel) = mean((dzTracerSpline(zLevelIndices,:)).^2,1)';
    
    zD2Linear(:,zLevel) = mean((dzLinear(zLevelIndices,:)-dz(zLevelIndices,:)).^2,1)';
    xyD2Linear(:,zLevel) = mean((xLinear(zLevelIndices,:)-x(zLevelIndices,:)).^2+(yLinear(zLevelIndices,:)-y(zLevelIndices,:)).^2,1)';
    
    zD2Spline(:,zLevel) = mean((dzSpline(zLevelIndices,:)-dz(zLevelIndices,:)).^2,1)';
    xyD2Spline(:,zLevel) = mean((xSpline(zLevelIndices,:)-x(zLevelIndices,:)).^2+(ySpline(zLevelIndices,:)-y(zLevelIndices,:)).^2,1)';
    
    % Diffusivity associated with the interpolation method, at each level
    % This is measured by looking at the two-particle separation
    % statistics.
    xLinearInterpDiffusivity(:,zLevel) = mean((xLinear(zLevelIndices,:)-x(zLevelIndices,:)).^2,1)'./t;
    yLinearInterpDiffusivity(:,zLevel) = mean((yLinear(zLevelIndices,:)-y(zLevelIndices,:)).^2,1)'./t;
    zLinearInterpDiffusivity(:,zLevel) = mean((zLinear(zLevelIndices,:)-z(zLevelIndices,:)).^2,1)'./t;
    zIsoLinearInterpDiffusivity(:,zLevel) = mean((dzLinear(zLevelIndices,:)-dz(zLevelIndices,:)).^2,1)'./t;
    
    xSplineInterpDiffusivity(:,zLevel) = mean((xSpline(zLevelIndices,:)-x(zLevelIndices,:)).^2,1)'./t;
    ySplineInterpDiffusivity(:,zLevel) = mean((ySpline(zLevelIndices,:)-y(zLevelIndices,:)).^2,1)'./t;
    zSplineInterpDiffusivity(:,zLevel) = mean((zSpline(zLevelIndices,:)-z(zLevelIndices,:)).^2,1)'./t;
    zIsoSplineInterpDiffusivity(:,zLevel) = mean((dzSpline(zLevelIndices,:)-dz(zLevelIndices,:)).^2,1)'./t;
    
    % Physical diffusivity at each level
    kappa_z(:,zLevel) = (mean(dz(zLevelIndices,:).^2,1)'./t)/2;
    kappa_z_centered(:,zLevel) = (mean(dzTracer(zLevelIndices,:).^2,1)'./t)/2;
    kappa_z_centered_spline(:,zLevel) = (mean(dzTracerSpline(zLevelIndices,:).^2,1)'./t)/2;
    kappa_z_centered_linear(:,zLevel) = (mean(dzTracerLinear(zLevelIndices,:).^2,1)'./t)/2;
    [minD, maxD] = SecondMomentMatrix( x(zLevelIndices,:)', y(zLevelIndices,:)','eigen');
    D2 = (minD + maxD)/2 - (minD(1)+maxD(1))/2;
    kappa_h(:,zLevel) = 0.5 * D2./t;
    [Mxx, Myy] = SecondMomentMatrix( z(zLevelIndices,:)', zeros(size(z(zLevelIndices,:)')) );
end

kappa_h_interp = (xSplineInterpDiffusivity + ySplineInterpDiffusivity)/4;
kappa_h_interp_linear = (xLinearInterpDiffusivity + yLinearInterpDiffusivity)/4;
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
plot(t,kappa_z,'k--','LineWidth',2), ylog, hold on
plot(t,kappa_z_centered,'k','LineWidth',2)
plot(t,kappa_z_centered_spline)
plot(t,kappa_z_centered_linear)
% plot(t,kappa_z_interp)
% plot(t,kappa_z_interp_linear)
title(sprintf('%dx%dx%d',Nx,Ny,Nz))
xlabel('time (s)')
ylabel('diffusivity (m^2/s)')


figure
plot(t,kappa_h,'k','LineWidth',2), hold on
plot(t,kappa_h_interp), ylog
plot(t,kappa_h_interp_linear)

return

% figure
% zLevel = nFloatLevels;
% zLevelIndices = (zLevel-1)*floatsPerLevel + (1:floatsPerLevel);
% edges = -1:0.1:1;
% subplot(2,2,1)
% dzHist = dzTracer;
% histogram(dzHist(zLevelIndices,round(end/4)), edges)
% subplot(2,2,2)
% histogram(dzHist(zLevelIndices,round(2*end/4)), edges)
% subplot(2,2,3)
% histogram(dzHist(zLevelIndices,round(3*end/4)), edges)
% subplot(2,2,4)
% histogram(dzHist(zLevelIndices,end), edges)
% 
% figure
% zLevel = nFloatLevels;
% zLevelIndices = (zLevel-1)*floatsPerLevel + (1:floatsPerLevel);
% edges = -1:0.1:1;
% subplot(2,2,1)
% dzHist = dz;
% histogram(dzHist(zLevelIndices,round(end/4)), edges)
% subplot(2,2,2)
% histogram(dzHist(zLevelIndices,round(2*end/4)), edges)
% subplot(2,2,3)
% histogram(dzHist(zLevelIndices,round(3*end/4)), edges)
% subplot(2,2,4)
% histogram(dzHist(zLevelIndices,end), edges)
% 
% figure
% scatter(reshape(sqrt(zD2Spline),[],1),reshape(kappa_z_interp,[],1)), xlog, ylog, hold on
% scatter(reshape(sqrt(xyD2Spline),[],1),reshape(kappa_h_interp,[],1))
% scatter(reshape(sqrt(zD2Tracer),[],1),reshape(kappa_z,[],1))
% title('Diffusivity of spline interpolation, total time')
% xlabel('rms separation (m)')
% ylabel('diffusivity (m^2/s)')
% legend('vertical', 'horizontal', 'vertical physical')
% 
% figure
% scatter(reshape(sqrt(zD2Linear),[],1),reshape(kappa_z_interp_linear,[],1)), xlog, ylog, hold on
% scatter(reshape(sqrt(xyD2Linear),[],1),reshape(kappa_h_interp_linear,[],1))
% scatter(reshape(sqrt(zD2Tracer),[],1),reshape(kappa_z,[],1))
% title('Diffusivity of linear interpolation, total time')
% xlabel('rms separation (m)')
% ylabel('diffusivity (m^2/s)')
% legend('vertical', 'horizontal', 'vertical physical')

%%%%%%%%%%%%%%%%%%
% Notes:
% Interpolation diffusivity is a measure of how much an interpolated
% particle separates from the true particle.
%
% Physical/true diffusivity is a measure of how much a particle separates
% from its neighbor.
%
% 

% timeIndices = 1:10:length(t);
% t2 = t(timeIndices);
% D2 = zD2TracerSpline(timeIndices,:);
% 
% x = reshape(sqrt(D2(1:end-1,:)),[],1);
% y = reshape(diff(D2)./diff(t2),[],1);
% 
% [N,edges,bin] = histcounts(log10(x),15);
% yMean = zeros(max(bin),1);
% yStdErr = zeros(max(bin),1);
% for i=1:max(bin)
%     yMean(i) = mean(y(bin==i));
%     yStdErr(i) = std(y(bin==i))/sqrt(sum(bin==i));
% end
% xMean = ((edges(1:end-1)+edges(2:end))/2)';
% 
% figure
% errorbar(xMean, yMean,2*yStdErr)


figure
scatter(reshape(sqrt(zD2Spline(1:end-1,:)),[],1),reshape(diff(zD2Spline)./diff(t),[],1)), xlog, ylog, hold on
scatter(reshape(sqrt(xyD2Spline(1:end-1,:)),[],1),reshape(diff(xyD2Spline)./diff(t),[],1))
scatter(reshape(sqrt(zD2Tracer(1:end-1,:)),[],1),reshape(diff(zD2Tracer)./diff(t),[],1))
scatter(reshape(sqrt(zD2TracerSpline(1:end-1,:)),[],1),reshape(diff(zD2TracerSpline)./diff(t),[],1))
title('Diffusivity of spline interpolation, local time')
xlabel('rms separation (m)')
ylabel('diffusivity (m^2/s)')
legend('vertical', 'horizontal', 'vertical physical', 'vertical physical spline')

figure
scatter(reshape(sqrt(zD2Linear(1:end-1,:)),[],1),reshape(diff(zD2Linear)./diff(t),[],1)), xlog, ylog, hold on
scatter(reshape(sqrt(xyD2Linear(1:end-1,:)),[],1),reshape(diff(xyD2Linear)./diff(t),[],1))
scatter(reshape(sqrt(zD2Tracer(1:end-1,:)),[],1),reshape(diff(zD2Tracer)./diff(t),[],1))
title('Diffusivity of linear interpolation, local time')
xlabel('rms separation (m)')
ylabel('diffusivity (m^2/s)')
legend('vertical', 'horizontal', 'vertical physical')

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
