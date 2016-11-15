file = '/Volumes/Data/InternalWaveSimulations/InternalWavesGMSpectrumWeakFlow_64_64_65.nc';
file = '/Volumes/Data/InternalWavesLatmix_256_256_50_GM_0.013.nc';
file = '/Volumes/jearly/Desktop/InternalWavesLatmix_256_256_50_GM_0.042.nc';
file = '/Volumes/Data/InternalWaveSimulations/InternalWavesGMSpectrumExponentialStratification.nc';
file = '/Volumes/home/jearly/InternalWavesLatmix_256_256_50_GM_0.062.nc';
file = '/Volumes/Data/InternalWavesLatmix_256_64_80_GM_0.125.nc';
file = '/Volumes/home/jearly/InternalWavesLatmixStrained_256_64_80_GM_0.016.nc';
file = '/Volumes/Scratch/InternalWavesGMSpectrumExponentialStratification.nc';
%file = '/Volumes/Data/InternalWavesLatmix_256_256_50_GM_0.062.nc';

x = ncread(file, 'x');
y = ncread(file, 'y');
z = ncread(file, 'z');
t = ncread(file, 'time');

rho_bar = double(ncread(file, 'rho_bar'));
N2 = double(ncread(file, 'N2'));
f0 = ncreadatt(file, '/', 'f0');
latitude = ncreadatt(file, '/', 'latitude');

%N2 = double(ncread(file, 'N2'));
%[~, ~, ~, N2] = InternalWaveModesFromDensityProfile_Spectral( rho_bar, z, z, 0.0, latitude, 'total_energy', 'rigid_lid' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check the frequency spectrum at a given point
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

depth = -50;
[depth_index] = find(z <= depth, 1, 'first');

stride = 4;
t_index = length(t)-1;
u3d = double(squeeze(ncread(file, 'u', [1 1 depth_index 1], [length(y)/stride length(x)/stride 1 t_index], [stride stride 1 1])));
v3d = double(squeeze(ncread(file, 'v', [1 1 depth_index 1], [length(y)/stride length(x)/stride 1 t_index], [stride stride 1 1])));

dt = t(2)-t(1)

[M, N, K] = size(u3d);

% Compute a few 'mooring' time series
cv_mooring = zeros([K 1]);
subsample = 1;
iMooring = 0;
for i=1:subsample:M
	for j=1:subsample:N
		iMooring = iMooring+1;
		cv_mooring(:,iMooring) = squeeze(u3d(i,j,:) + sqrt(-1)*v3d(i,j,:));
	end
end

taper_bandwidth = 2;
psi=[];
%  [psi,lambda]=sleptap(size(cv_mooring,1),taper_bandwidth);
[omega_p, Spp, Snn, Spn] = mspec(dt,cv_mooring,psi);

omega = [ -flipud(omega_p(2:end)); omega_p];
% We want the integral of this to give us the variance back, so we need to
% divide by 2*pi
S = (1/(2*pi))*[flipud(vmean(Snn,2)); vmean(Spp(2:end,:),2)];
%S = (1/(2*pi))*[flipud(Snn(:,1)); Spp(2:end,2)];

N0=5.23E-3;
%[S_gm] = 2.5*GarrettMunkHorizontalKineticEnergyRotarySpectrumWKB( omega, latitude, sqrt(N2(depth_index)) );
[S_gm2] = 0.5*GarrettMunkHorizontalKineticEnergyRotarySpectrum( omega, latitude, z, rho_bar, 3, depth, 0 );

% Factor of two to get variance, instead of energy
[S_gm] = GarrettMunkHorizontalKineticEnergyRotarySpectrumWKB( omega, latitude, N0, 0 );
%S_gm(find(isnan(S_gm))) = 0;
%S_gm = BlurSpectrum( omega, S_gm);
% S_gm2 = BlurSpectrum( omega, S_gm2);

figure
plot( omega, N0/sqrt(N2(depth_index))*S, 'blue', 'LineWidth', 2), ylog
hold on
plot( omega, S_gm, 'black', 'LineWidth', 2), ylog
plot( omega, S_gm2, 'green', 'LineWidth', 2), ylog
vlines( -f0 )
legend('Observed', 'GM81 WKB', 'GM81 JJE')

df = omega(2)-omega(1);
sum(S_gm)*df
sum(S)*df

restrictedIndices = find( abs(omega) > 2*f0 & abs(omega) < 0.75*max(omega) );

sum(S_gm(restrictedIndices))*df
sum(S(restrictedIndices))*df

return;

t=t(1:t_index);
figure
plot(t,real(cv_mooring(:,iMooring)))
hold on
plot(t,imag(cv_mooring(:,iMooring)))

sn = vmean(sn,2);
negativeF = find(fn<0);
positiveF = find(fn>0);
if (length(negativeF) > length(positiveF))
	negativeF(1) = [];
end
s1sided = cat(1, sn(find(fn==0)), flipud(sn(negativeF))+sn(positiveF));
f1sided = fn(find(fn>=0));
figure, plot(fn, sn), ylog
hold on, plot(f1sided, s1sided, 'k')
plot(-f1sided, s1sided, 'k')