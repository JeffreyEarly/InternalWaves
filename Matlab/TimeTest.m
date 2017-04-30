%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DiffusivityExperiment
%
% This script uses the InternalWaveModel to generate the time-evolution of
% a Garrett-Munk spectrum of internal waves and save the output to a NetCDF
% file. Advects particles as well.
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com
%
% April 12th, 2017      Version 1.0

N = 64;
aspectRatio = 8;

L = 15e3;
Lx = aspectRatio*L;
Ly = L;
Lz = 5000;

Nx = aspectRatio*N;
Ny = N;
Nz = (N/2)+1; % Must include end point to advect at the surface, so use 2^N + 1

latitude = 31;
N0 = 5.2e-3; % Choose your stratification
GMReferenceLevel = 1.0;

kappa = 5e-6;
outputInterval = 15*60;
maxTime = 86400/4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shouldUseGMSpectrum = 1;

wavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

if shouldUseGMSpectrum == 1
    wavemodel.FillOutWaveSpectrum();
    wavemodel.InitializeWithGMSpectrum(GMReferenceLevel);
    wavemodel.ShowDiagnostics();
    period = 2*pi/wavemodel.N0;
    [u,v] = wavemodel.VelocityFieldAtTime(0.0);
    U = max(max(max( sqrt(u.*u + v.*v) )));
else
    j0 = 1; % j=1..nModes, where 1 indicates the 1st baroclinic mode
    U = 0.1; % m/s
    sign = 1;
    phi = 0;
    k0 = 2;
    l0 = 0;
    alpha = atan2(l0,k0);
    k = 2*pi*sqrt(k0^2 + l0^2)/Lx;
    
    period = wavemodel.InitializeWithPlaneWave(k0,l0,j0,U,sign);
    maxTime = period;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create floats/drifters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx = wavemodel.x(2)-wavemodel.x(1);
dy = wavemodel.y(2)-wavemodel.y(1);
N = 10;
x_float = (1:N)*dx;
y_float = (1:N)*dy;
z_float = (0:2)*(-Lz/4);

[x_float,y_float,z_float] = ndgrid(x_float,y_float,z_float);
x_float = reshape(x_float,[],1);
y_float = reshape(y_float,[],1);
z_float = reshape(z_float,[],1);
nFloats = numel(x_float);

% Now let's place the floats along an isopycnal.
[w,zeta] = wavemodel.VerticalFieldsAtTime(0);
isopycnalDeviation = interpn(wavemodel.X,wavemodel.Y,wavemodel.Z,zeta,x_float,y_float,z_float,'spline');
z_isopycnal = z_float + isopycnalDeviation;

ymin = [-Inf -Inf -Lz -Inf -Inf -Lz -Inf -Inf];
ymax = [Inf Inf 0 Inf Inf 0 Inf Inf];
kappa_vector = [0 0 0 kappa kappa kappa 0 0];
p0 = cat(2, x_float, y_float, z_isopycnal, x_float, y_float, z_isopycnal, x_float, y_float);

f = @(t,y) FluxForFloatDiffusiveDrifter(t,y,z_float,wavemodel, 'spline');

tic
for i=1:1
   a=f(0,p0) ;
end
toc
