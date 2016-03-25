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

Lx = 40e3;
Ly = 40e3;
Lz = 5000;

Nx = 32;
Ny = 32;
Nz = 64;

latitude = 35;
N0 = 5.2e-3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavemodel = InternalWaveModel([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create a single plane-wave with the model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k0 = 3; % k=0..Nx/2
l0 = 0; % l=0..Ny/2
j0 = 1; % j=1..nModes, where 1 indicates the 1st baroclinic mode
U = 1.0; % m/s
sign = +1;

wavemodel.InitializeWithPlaneWave(k0,l0,j0,U,sign);

t = 500;
[u,v] = wavemodel.VelocityFieldAtTime(t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create a single plane-wave with the known analytical solution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = wavemodel.k;
X = wavemodel.X;
Z = wavemodel.Z;

% Note that this solution is only valid for l=0.
omega = 2*pi/wavemodel.period;
u_unit = Amplitude*cos( k(k0+1)*X + omega*t ).*cos((j0*pi/Lz)*Z);
v_unit = -(f0/omega)*Amplitude*sin( k(k0+1)*X + omega*t ).*cos((j0*pi/Lz)*Z);

max_speed = max(max(max( sqrt(u.*u + v.*v) )));
u_error = max(max(max(abs(u-u_unit)/max_speed)));
v_error = max(max(max(abs(v-v_unit)/max_speed)));

fprintf('The model solution for (u,v) matches the analytical solution to 1 part in (10^%d, 10^%d)\n', round(abs(log10(u_error))), round(abs(log10(v_error))));