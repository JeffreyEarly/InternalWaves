file = '/Volumes/OceanTransfer/InternalWaveModel_2017-04-13T120332_512x64x33.nc';
x = ncread(file,'x-position')';
y = ncread(file,'y-position')';
z = ncread(file,'z-position')';
% x = ncread(file,'x-position-diffusive')';
% y = ncread(file,'y-position-diffusive')';
% z = ncread(file,'z-position-diffusive')';
% x = ncread(file,'x-position-drifter')';
% y = ncread(file,'y-position-drifter')';
% z = ncread(file,'z-position-drifter')';
t = ncread(file,'t');

indices = (1:100) + 0*100;


x = x(:,indices);
y = y(:,indices);
z = z(:,indices);

[x_com, y_com, q, r] = CenterOfMass( x, y );

[minD, maxD, theta] = SecondMomentMatrix( x, y, 'eigen' );

D2 = (minD+maxD)/2;
[p,S,mu]=polyfit(t,D2,1);
kappa_fit = 0.5*p(1)/mu(2);
fprintf('diffusive linear fit: kappa = %f\n', kappa_fit)

figure
plot(q,r)

figure
plot(t/86400, D2)