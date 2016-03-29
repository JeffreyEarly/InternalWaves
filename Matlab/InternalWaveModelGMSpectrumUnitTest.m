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
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx = 80e3;
Ly = 80e3;
Lz = 5000;

Nx = 256;
Ny = 256;
Nz = 64;

latitude = 35;
N0 = 5.2e-3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavemodel = InternalWaveModel([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);
wavemodel.InitializeWithGMSpectrum(1.0);


[u,v]=wavemodel.VelocityFieldAtTime(0);
[w,zeta] = wavemodel.VerticalFieldsAtTime(0);

z = wavemodel.z;
HKE = u.*u + v.*v;
AvgHKE = mean(mean(mean(HKE)))*1e4;
fprintf('The average HKE is %f cm^2/s^2, compared to 44 cm^2/s^2 for WKB scaled GM.\n',AvgHKE);


HKE_int = trapz(z,HKE,3);

VKE = w.*w;
VKE_int = trapz(z,VKE,3);

PE = (N0^2)*zeta.*zeta;
PE_int = trapz(z,PE,3);