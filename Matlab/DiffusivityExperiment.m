Lx = 15e3;
Ly = 15e3;
Lz = 5000;

Nx = 64;
Ny = 64;
Nz = 65; % Must include end point to advect at the surface, so use 2^N + 1

latitude = 31;
N0 = 5.2e-3; % Choose your stratification

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shouldUseGMSpectrum = 0;

wavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

if shouldUseGMSpectrum == 1
    wavemodel.InitializeWithGMSpectrum(1.0);
    maxTime = 60*60; %2*pi/wavemodel.f0;
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

dx = wavemodel.x(2)-wavemodel.x(1);
dy = wavemodel.y(2)-wavemodel.y(1);
N = 10;
x_float = (1:N)*dx;
y_float = (1:N)*dy;
z_float = (0:2)*(-Lz/4);

[x,y,z] = ndgrid(x_float,y_float,z_float);
p0 = cat(2, reshape(x,[],1), reshape(y,[],1), reshape(z,[],1));

f = @(t,y) wavemodel.VelocityAtTimePositionVector(t,y);

% Let's do fixed step size integrator.
cfl = 0.25;
advectiveDT = cfl*(wavemodel.x(2)-wavemodel.x(1))/U;
oscillatoryDT = period/8;
if advectiveDT < oscillatoryDT
    fprintf('Using the advective dt: %.2f\n',advectiveDT);
    deltaT = advectiveDT;
else
    fprintf('Using the oscillatory dt: %.2f\n',oscillatoryDT);
    deltaT = oscillatoryDT;
end

t = (0:deltaT:maxTime)'; %(0:60*ceil(deltaT/60):maxTime)';
if t(end) < period
    t(end+1) = period;
end

kappa = 1e-6;
p = ode4kVector(f, kappa*[1 1 1], [-Inf -Inf -Lz], [Inf Inf 0], t, p0);
x = squeeze(p(:,1,:));
y = squeeze(p(:,2,:));
z = squeeze(p(:,3,:));

figure, plot(x',y')