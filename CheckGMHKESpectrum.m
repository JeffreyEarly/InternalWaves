file = '/Volumes/Data/InternalWaveSimulations/InternalWavesGMSpectrumWeakFlow_64_64_65.nc';

x = ncread(file, 'x');
y = ncread(file, 'y');
z = ncread(file, 'z');
t = ncread(file, 'time');

rho_bar = double(ncread(file, 'rho_bar'));

f0 = ncreadatt(file, '/', 'f0');
latitude = ncreadatt(file, '/', 'latitude');

%N2 = double(ncread(file, 'N2'));
[~, ~, ~, N2] = InternalWaveModesFromDensityProfile_Spectral( rho_bar, z, z, 0.0, latitude, 'total_energy', 'rigid_lid' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check the frequency spectrum at a given point
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

depth = -50;
[depth_index] = find(z >= depth, 1, 'first');

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

taper_bandwidth = 1;
psi=[];
%[psi,lambda]=sleptap(size(cv_mooring,1),taper_bandwidth);
[omega_p, Spp, Snn, Spn] = mspec(dt,cv_mooring,psi);

omega = [ -flipud(omega_p(2:end)); omega_p];
S = [flipud(vmean(Snn,2)); vmean(Spp(2:end,:),2)];

%[S_gm] = 2.5*GarrettMunkHorizontalKineticEnergyRotarySpectrumWKB( omega, latitude, sqrt(N2(depth_index)) );
%[S_gm] = (2/3.14)*GarrettMunkHorizontalKineticEnergyRotarySpectrum( omega, latitude, z, rho_bar, depth );
[S_gm] = 1.0*GarrettMunkHorizontalKineticEnergyRotarySpectrumWKB( omega, latitude, sqrt(N2(depth_index)) );
S_gm = BlurSpectrum( omega, S_gm);

figure
plot( omega, S, 'blue', 'LineWidth', 2), ylog
hold on
plot( omega, S_gm, 'black', 'LineWidth', 2), ylog
vlines( -f0 )
legend('Observed', 'GM81')

return;

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