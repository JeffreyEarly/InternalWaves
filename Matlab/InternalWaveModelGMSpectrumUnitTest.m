%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
HKE_int = trapz(z,HKE,3);

VKE = w.*w;
VKE_int = trapz(z,VKE,3);

PE = (N0^2)*zeta.*zeta;
PE_int = trapz(z,PE,3);


mean(mean(HKE_int))*1e4
mean(mean(VKE_int))*1e4
mean(mean(PE_int))*1e4

L_gm = 1.3e3; % thermocline exponential scale, meters
invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
E_gm = 6.3e-5; % non-dimensional energy parameter
E = L_gm*L_gm*L_gm*invT_gm*invT_gm*E_gm;