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
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx = 5e3;
Ly = 5e3;
Lz = 5000;

Nx = 64;
Ny = 64;
Nz = 64;

latitude = 31;
N0 = 5.2e-3; % Choose your stratification

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

z = wavemodel.z;
uvVariance = squeeze(mean(mean(u.*u + v.*v,1),2));
zetaVariance = squeeze(mean(mean(zeta.*zeta,1),2));

figure
subplot(1,2,1)
plot([44 44], [z(1) z(end)], 'k' ,'LineWidth', 2), hold on
plot(1e4*uvVariance,z,'LineWidth', 2)
xlabel('cm^2/s^2'), ylabel('depth (m)')

subplot(1,2,2)
plot([53 53], [z(1) z(end)], 'k' ,'LineWidth', 2), hold on
plot(zetaVariance,z,'LineWidth', 2)
xlabel('m^2'), ylabel('depth (m)')

HKE = u.*u + v.*v;
HKE_int = trapz(z,HKE,3);

VKE = w.*w;
VKE_int = trapz(z,VKE,3);

PE = (N0^2)*zeta.*zeta;
PE_int = trapz(z,PE,3);

totalGM = mean(mean(HKE_int + VKE_int + PE_int))*1032; % scaled by the density of water
fprintf('The total energy in the water column is %f J/m^2, compared to 3800 J/m^2 expected for GM.\n',totalGM);

AvgHKE = mean(mean(HKE(:,:,end)))*1e4;
fprintf('The average 2*HKE is %f cm^2/s^2 at the surface, compared to 44 cm^2/s^2 for WKB scaled GM.\n',AvgHKE);