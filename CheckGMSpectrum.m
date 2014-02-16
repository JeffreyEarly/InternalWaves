L_GM = 1.3E3;       % [m]
invT_GM = 5.2E-3;   % [1/s]
E_GM = 6.3E-5;      % [unitless]

E_GM_total = L_GM*L_GM*L_GM*invT_GM*invT_GM*E_GM;

file = '/Users/jearly/Desktop/InternalWavesGM.nc';

x = ncread(file, 'x');
y = ncread(file, 'y');
z = ncread(file, 'z');
t = ncread(file, 'time');
rho_bar = double(ncread(file, 'rho_bar'));
dz=z(2)-z(1);
N2 = (-9.81/vmean(rho_bar,1))*vdiff(dz,rho_bar,1);

iTime=1;
zeta3d = double(squeeze(ncread(file, 'zeta', [1 1 1 iTime], [length(y) length(x) length(z) 1], [1 1 1 1])));
u3d = double(squeeze(ncread(file, 'u', [1 1 1 iTime], [length(y) length(x) length(z) 1], [1 1 1 1])));
v3d = double(squeeze(ncread(file, 'v', [1 1 1 iTime], [length(y) length(x) length(z) 1], [1 1 1 1])));


E_p = 0.5*dz*trapz(N2.*squeeze(vmean(vmean(zeta3d.*zeta3d,1),2)));

E_k = 0.5*dz*trapz(squeeze(vmean(vmean(u3d.*u3d+v3d.*v3d,1),2)));

potential_kinetic_ratio = E_p/E_k
GM_relative = (E_p+E_k)/E_GM_total


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check the frequency spectrum at a given point
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u3d = double(squeeze(ncread(file, 'u', [1 1 64 1], [length(y) length(x) 1 length(t)], [1 1 1 1])));
v3d = double(squeeze(ncread(file, 'v', [1 1 64 1], [length(y) length(x) 1 length(t)], [1 1 1 1])));

f0 = ncreadatt(file, '/', 'f0');
dt = t(2)-t(1)

[M, N, K] = size(u3d);

% Compute a few 'mooring' time series
cv_mooring = zeros([length(t) 1]);
subsample = 8;
iMooring = 0;
for i=1:subsample:M
	for j=1:subsample:N
		iMooring = iMooring+1;
		cv_mooring(:,iMooring) = squeeze(u3d(i,j,:) + sqrt(-1)*v3d(i,j,:));
	end
end

[fn, sn] = powspec(dt,cv_mooring);
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


fn = f1sided;
sn = s1sided;
% grab indices between the Coriolis and half the nyquist
fitIndices = find( fn > 1.5*f0/(2*pi) & fn < 0.5*1/(2*dt));
[P,S] = polyfit(log10(fn(fitIndices)), log10(sn(fitIndices)),1);
fit = 10^(P(2))*fn.^(P(1));
hold on
plot( fn(fitIndices), fit(fitIndices), 'r')
slope = P(1)

[psi,lambda]=sleptap(size(cv_mooring,1),3);
[fn,spp,snn,spn]=mspec(dt,cv_mooring,psi);
fn=fn/(2*pi);
figure, plot(fn, vmean(spp,2)), ylog
hold on, plot(fn, vmean(snn,2))
plot(fn,vmean(spp+snn,2))

sn = vmean(spp+snn,2);
% grab indices between the Coriolis and half the nyquist
fitIndices = find( fn > 1.5*f0/(2*pi) & fn < 0.5*1/(2*dt));
[P,S] = polyfit(log10(fn(fitIndices)), log10(sn(fitIndices)),1);
fit = 10^(P(2))*fn.^(P(1));
hold on
plot( fn(fitIndices), fit(fitIndices), 'r')
slope = P(1)


xindices = 1:2:length(x);
yindices = 1:2:length(y);
figure
iTime = 2;
quiver(x(xindices),y(yindices),u3d(xindices,yindices,iTime),v3d(xindices,yindices,iTime),0.8)
xlim([min(x) max(x)])
ylim([min(y) max(y)])