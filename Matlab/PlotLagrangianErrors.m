file = '/Users/jearly/Desktop/LagrangianErrorExperiment_2017-05-02T125750_256x32x33.nc';
% file = '/Users/jearly/Desktop/LagrangianErrorExperiment_2017-05-02T112657_128x16x17.nc';
file = '/Users/jearly/Desktop/LagrangianErrorExperiment_2017-05-02T133536_128x16x33.nc';
% file = '/Users/jearly/Desktop/LagrangianErrorExperiment_2017-05-02T124907_64x8x9.nc';

t = ncread(file, 't');

x = ncread(file, 'x-position-exact');
y = ncread(file, 'y-position-exact');
z = ncread(file, 'z-position-exact');

xLinear = ncread(file, 'x-position-linear');
yLinear = ncread(file, 'y-position-linear');
zLinear = ncread(file, 'z-position-linear');

xSpline = ncread(file, 'x-position-spline');
ySpline = ncread(file, 'y-position-spline');
zSpline = ncread(file, 'z-position-spline');

% xSplineDiffusion = mean((xSpline(:,end)-x(:,end)).^2)/t(end)
% ySplineDiffusion = mean((ySpline(:,end)-y(:,end)).^2)/t(end)
% zSplineDiffusion = mean((zSpline(:,end)-z(:,end)).^2)/t(end)

% drifterIndices = 101;
% figure
% plot( x(drifterIndices,:), y(drifterIndices,:) ), hold on
% plot( xSpline(drifterIndices,:), ySpline(drifterIndices,:) )


xSplineDiffusion = mean((xSpline-x).^2,1)'./t;
ySplineDiffusion = mean((ySpline-y).^2,1)'./t;
zSplineDiffusion = mean((zSpline(101:300,:)-z(101:300,:)).^2,1)'./t;

kappa_h = (xSplineDiffusion + xSplineDiffusion)/4;
kappa_z = zSplineDiffusion/2;

figure
plot(t,[kappa_h kappa_z]), ylog
legend('horizontal diffusivity', 'vertical diffusivity')

kappa_z0000 = mean((zSpline(1:100,:)-z(1:100,:)).^2,1)'./t;
kappa_z1250 = mean((zSpline(101:200,:)-z(101:200,:)).^2,1)'./t;
kappa_z2500 = mean((zSpline(201:300,:)-z(201:300,:)).^2,1)'./t;
figure
plot(t, [kappa_z1250 kappa_z2500]), ylog
legend('1250 m', '2500 m')
