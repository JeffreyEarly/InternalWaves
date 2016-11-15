%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% InternalWaveModelPlaneWaveUnitTest
%
% This script uses the InternalWaveModel to create, and validate, a single
% internal wave.
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

k0 = -1; % k=0..Nx/2
l0 = 0; % l=0..Ny/2
j0 = 1; % j=1..nModes, where 1 indicates the 1st baroclinic mode
U = 0.01; % m/s
sign = -1;

wavemodel.InitializeWithPlaneWave(k0,l0,j0,U,sign);

t = 86400/4;
[u,v] = wavemodel.VelocityFieldAtTime(t);
[w,zeta] = wavemodel.VerticalFieldsAtTime(t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create a single plane-wave with the known analytical solution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f0 = wavemodel.f0;
k = wavemodel.k;
h = wavemodel.h;
x = wavemodel.x;
X = wavemodel.X;
Z = wavemodel.Z;

% Note that this solution is only valid for l=0.
if (k0 < 0)
    k0 = Nx + k0;
end
omega = sign*2*pi/wavemodel.period;
m = j0*pi/Lz;
u_unit = U*cos( k(k0+1)*X + omega*t ).*cos(m*Z);
v_unit = -(f0/omega)*U*sin( k(k0+1)*X + omega*t ) .* cos(m*Z);
w_unit = (U*k(k0+1)/m) * sin(k(k0+1)*X + omega*t) .* sin(m*Z);
zeta_unit = -(U*k(k0+1)/m/omega) * cos(k(k0+1)*X + omega*t) .* sin(m*Z);

% Compute the relative error
max_speed = max(max(max( sqrt(u.*u + v.*v) )));
max_u = max( [max(max(max( u ))), 1e-15] );
u_error = max(max(max(abs(u-u_unit)/max_u)));
max_v = max( [max(max(max( v ))), 1e-15] );
v_error = max(max(max(abs(v-v_unit)/max_v)));
max_w = max( [max(max(max( abs(w) ))), 1e-15] );
w_error = max( [max(max(max(abs(w-w_unit)/max_w))), 1e-15] );
max_zeta = max( [max(max(max( zeta ))), 1e-15] );
zeta_error = max( [max(max(max(abs(zeta-zeta_unit)/max_zeta))), 1e-15] );

fprintf('The model solution for (u,v) matches the analytical solution to 1 part in (10^%d, 10^%d) at time t=%d\n', round((log10(u_error))), round((log10(v_error))),t);
fprintf('The model solution for (w,zeta) matches the analytical solution to 1 part in (10^%d, 10^%d) at time t=%d\n', round((log10(w_error))), round((log10(zeta_error))),t);


% figure, plot( squeeze(w(8,8,:)), z)
% hold on, plot( squeeze(w_unit(8,8,:)), z)