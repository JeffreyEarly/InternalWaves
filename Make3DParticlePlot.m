file = '/Users/jearly/Desktop/InternalWaves.nc';
file = '/Users/jearly/Desktop/InternalWavesLatmix_64_64_25.nc'

xpos=double(ncread(file, 'x-position'));
ypos=double(ncread(file, 'y-position'));
zpos=double(ncread(file, 'z-position'));

i=5;
j=5;
k=2;

figure
plot3( squeeze(xpos(i,j,k,:)), squeeze(ypos(i,j,k,:)), squeeze(zpos(i,j,k,:)))