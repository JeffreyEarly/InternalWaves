file = '/Users/jearly/Desktop/InternalWavesConstantN_256_256_128_lat31.nc';

x = ncread(file, 'x');
y = ncread(file, 'y');
z = ncread(file, 'z');
t = ncread(file, 'time');
rho_bar = double(ncread(file, 'rho_bar'));
N2 = double(ncread(file, 'N2'));

latitude = 31;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check the frequency spectrum at a given point
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

depth = -150;
[depth_index] = find(z >= depth, 1, 'first');

u3d = double(squeeze(ncread(file, 'u', [1 1 depth_index 1], [length(y) length(x) 1 length(t)], [1 1 1 1])));
v3d = double(squeeze(ncread(file, 'v', [1 1 depth_index 1], [length(y) length(x) 1 length(t)], [1 1 1 1])));

f0 = ncreadatt(file, '/', 'f0');
dt = t(2)-t(1)

[M, N, K] = size(u3d);

% Compute a few 'mooring' time series
cv_mooring = zeros([length(t) 1]);
subsample = 4;
iMooring = 0;
for i=1:subsample:M
	for j=1:subsample:N
		iMooring = iMooring+1;
		cv_mooring(:,iMooring) = squeeze(u3d(i,j,:) + sqrt(-1)*v3d(i,j,:));
	end
end

taper_bandwidth = 3;
psi=[];
%[psi,lambda]=sleptap(size(cv_mooring,1),taper_bandwidth);
[omega_p, Spp, Snn, Spn] = mspec(dt,cv_mooring,psi);

omega = [ -flipud(omega_p(2:end)); omega_p];
S = [flipud(vmean(Snn,2)); vmean(Spp(2:end,:),2)];

[S_gm] = GarrettMunkHorizontalKineticEnergyRotarySpectrum( omega, latitude, sqrt(N2(depth_index)) );


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