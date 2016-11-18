%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% InternalWaveModelGMSpectrumUnitTest
%
% This script uses the InternalWaveModel to create, and validate, a
% Garrett-Munk spectrum in a linear internal wave field with constant
% stratification.
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com
%
% March 25th, 2016      Version 1.0
% March 30th, 2016      Version 1.1
% November 17th, 2016   Version 1.2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx = 30e3;
Ly = 15e3;
Lz = 5000;

Nx = 256;
Ny = 128;
Nz = 64;

latitude = 31;
N0 = 5.2e-3/2; % Choose your stratification 7.6001e-04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavemodel = InternalWaveModel([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);
wavemodel.InitializeWithGMSpectrum(1.0);

t = 0;
[u,v]=wavemodel.VelocityFieldAtTime(t);
[w,zeta] = wavemodel.VerticalFieldsAtTime(t);

depth = 1000;
depthIndex = find(wavemodel.z-Lz > -depth,1,'first');
stride = 4;
t = 0:15*60:4*86400;
xIndices = 1:stride:Nx;
yIndices = 1:stride:Ny;
cv_mooring = zeros([length(t) length(xIndices)*length(yIndices)]);
for iTime=1:length(t)
    fprintf('time @ %d of %d\n', iTime, length(t));
    [u,v]=wavemodel.VelocityFieldAtTime(t(iTime));
    cv_mooring(iTime,:) = reshape(u(xIndices,yIndices,depthIndex),1,[]) + sqrt(-1)*reshape(v(xIndices,yIndices,depthIndex),1,[]);
end

taper_bandwidth = 2;
psi=[];
%  [psi,lambda]=sleptap(size(cv_mooring,1),taper_bandwidth);
[omega_p, Spp, Snn, Spn] = mspec(t(2)-t(1),cv_mooring,psi);

omega = [ -flipud(omega_p(2:end)); omega_p];
[S_gm] = GarrettMunkHorizontalKineticEnergyRotarySpectrumWKB( omega, latitude, N0, 0 );
% We want the integral of this to give us the variance back, so we need to
% divide by 2*pi
S = (1/(2*pi))*[flipud(vmean(Snn,2)); vmean(Spp(2:end,:),2)];
figure, plot(omega,S), ylog
hold on, plot(omega,S_gm/5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute a few benchmarks
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = wavemodel.z;
uvVariance = squeeze(mean(mean(u.*u + v.*v,1),2));
zetaVariance = squeeze(mean(mean(zeta.*zeta,1),2));

zeta2 = squeeze(mean(mean(zeta.*zeta,1),2));
u2 = squeeze(mean(mean(u.*u,1),2)+mean(vmean(v.*v,1),2));
w2 = squeeze(mean(mean(w.*w,1),2));

figure
subplot(1,3,1)
plot([44 44], [z(1) z(end)], 'k' ,'LineWidth', 2), hold on
plot(1e4*uvVariance,z,'LineWidth', 2)
xlabel('cm^2/s^2'), ylabel('depth (m)')

subplot(1,3,2)
plot([53 53], [z(1) z(end)], 'k' ,'LineWidth', 2), hold on
plot(zetaVariance,z,'LineWidth', 2)
xlabel('m^2'), ylabel('depth (m)')

subplot(1,3,3)
plot([30 30], [z(1) z(end)], 'k' ,'LineWidth', 2), hold on
plot(0.5*1e4*(u2 + w2 + N0*N0.*zeta2) ,z,'LineWidth', 2)
title('total energy'), xlabel('cm^2 s^{-2}'), ylabel('depth (m)')

HKE = 0.5*(u.*u + v.*v);
HKE_int = trapz(z,HKE,3);

VKE = 0.5*(w.*w);
VKE_int = trapz(z,VKE,3);

PE = 0.5*(N0^2)*zeta.*zeta;
PE_int = trapz(z,PE,3);

totalGM = mean(mean(HKE_int + VKE_int + PE_int))*1032; % scaled by the density of water
fprintf('The total energy in the water column is %f J/m^2, compared to 3800 J/m^2 expected for GM.\n',totalGM);

AvgHKE = mean(mean(HKE(:,:,end)))*1e4;
fprintf('The average 2*HKE is %f cm^2/s^2 at the surface, compared to 44 cm^2/s^2 for WKB scaled GM.\n',AvgHKE);