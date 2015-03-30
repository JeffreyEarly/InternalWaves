file = '/Users/jearly/Desktop/InternalWaves.nc';
file = '/Users/jearly/Desktop/InternalWavesLatmix_128_128_50_GM_0.013.nc'

xpos=double(ncread(file, 'x-position'));
ypos=double(ncread(file, 'y-position'));
zpos=double(ncread(file, 'z-position'));

i=11;
j=13;
k=1;

figure
plot3( squeeze(xpos(i,j,k,:)), squeeze(ypos(i,j,k,:)), squeeze(zpos(i,j,k,:)))