%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GMVelocitySpatialCorrelationPrediction
%

N = 64;
aspectRatio = 1;

L = 500e3;
Lx = aspectRatio*L;
Ly = L;
Lz = 5000;

Nx = aspectRatio*N;
Ny = N;
Nz = N+1; % Must include end point to advect at the surface, so use 2^N + 1

latitude = 31;
N0 = 5.2e-3; % Choose your stratification
GMReferenceLevel = 1.0;

wavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

wavemodel.FillOutWaveSpectrum();
wavemodel.InitializeWithGMSpectrum(GMReferenceLevel,0);

% Kh = sqrt(wavemodel.K.^2 + wavemodel.L.^2);
% k = wavemodel.k;
% dk = k(2)-k(1);
% k_cutoff = 4*dk;
% indices = Kh < k_cutoff & Kh ~= 0;

% Compute total HKE, at each wavenumber
U2 = (wavemodel.u_plus).*conj(wavemodel.u_plus) + (wavemodel.u_minus).*conj(wavemodel.u_minus);
% U2(indices) = 0;
U2 = sum(U2,3);

% Compute the horizontal wavenumber
Kh = sqrt(wavemodel.K.^2 + wavemodel.L.^2);
Kh = Kh(:,:,1);
L = 2*pi./Kh;

x = wavemodel.x;
x(end+1) = Inf;
r = zeros(size(x));
for i=2:length(x)
    r(i-1) = sum(sum(U2( L > x(i-1) & L <= x(i))));
end
correlation = flip(cumsum(flip(r)));

s = 1/1000;
% figure, plot(s*x,r)
figure, plot(s*x,correlation/correlation(1))
xlabel('km')