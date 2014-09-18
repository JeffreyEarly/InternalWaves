file = '/Users/jearly/Desktop/InternalWavesConstantN_UnitTest_128_128_65.nc';

outputfile = 'InternalWavesConstantN_UnitTest_high_wavenumber.mat';

x = ncread(file, 'x');
y = ncread(file, 'y');
z = ncread(file, 'z');
t = ncread(file, 'time');
rho_bar = double(ncread(file, 'rho_bar'));
N2 = double(ncread(file, 'N2'));
f0 = ncreadatt(file, '/', 'f0');

iTime=1;

zeta3d = double(squeeze(ncread(file, 'zeta', [1 1 1 iTime], [length(y) length(x) length(z) 1], [1 1 1 1])));
rho3d = double(squeeze(ncread(file, 'rho', [1 1 1 iTime], [length(y) length(x) length(z) 1], [1 1 1 1])));
u3d = double(squeeze(ncread(file, 'u', [1 1 1 iTime], [length(y) length(x) length(z) 1], [1 1 1 1])));
v3d = double(squeeze(ncread(file, 'v', [1 1 1 iTime], [length(y) length(x) length(z) 1], [1 1 1 1])));
w3d = double(squeeze(ncread(file, 'w', [1 1 1 iTime], [length(y) length(x) length(z) 1], [1 1 1 1])));

rho_prime = permute(rho3d, [2 1 3]);
u = permute(u3d, [2 1 3]);
v = permute(v3d, [2 1 3]);
w = permute(w3d, [2 1 3]);

rho_prime = rho_prime - repmat(permute(rho_bar,[3 2 1]), [length(x) length(y) 1]);

save( outputfile, 'f0', 'u', 'v', 'w', 'rho_prime', 'rho_bar' );