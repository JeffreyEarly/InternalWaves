file = '/Volumes/MoreStorage/DiffusivityExperiment_2020-01-23T163207_128x128x129.nc';
file = '/Volumes/MoreStorage/DiffusivityExperiment_2020-01-23T182720_64x64x65.nc';
file = '/Volumes/MoreStorage/DiffusivityExperiment_2020-01-23T184203_64x64x129.nc';
file = '/Volumes/MoreStorage/DiffusivityExperiment_2020-01-23T185956_64x64x129.nc';
file = '/Volumes/MoreStorage/DiffusivityExperiment_2020-01-23T192741_64x64x129.nc';
file = '/Volumes/MoreStorage/DiffusivityExperiment_2020-01-23T195242_32x32x129.nc';
file = '/Volumes/MoreStorage/DiffusivityExperiment_2020-01-23T201809_64x64x129.nc';
file = '/Volumes/MoreStorage/DiffusivityExperiment_2020-01-23T211833_32x32x129.nc';
file = '/Volumes/MoreStorage/DiffusivityExperiment_2020-01-23T224856_128x128x129.nc';

x = ncread(file,'x-position');
y = ncread(file,'y-position');
z = ncread(file,'z-position');
t = ncread(file,'t');

validIndices = 1:(length(t)-1);
x = x(:,validIndices);
y = y(:,validIndices);
z = z(:,validIndices);
t = t(validIndices);

x2 = ncread(file,'x-position-spline');
y2 = ncread(file,'y-position-spline');
z2 = ncread(file,'z-position-spline');

x2 = x2(:,validIndices);
y2 = y2(:,validIndices);
z2 = z2(:,validIndices);

dx = x2-x;
dy = y2-y;
figure, plot(t,mean((dx.^2 + dy.^2).',2))

x = x2;
y = y2;

D2_particles = zeros(length(t),1);
r2_particles = zeros(length(t),1);
nLevels = 5;
floatsPerLevel = size(x,1)/nLevels;
for zLevel = 1:nLevels
    nLevels = size(x,1)/floatsPerLevel;
    zLevelIndices = (zLevel-1)*floatsPerLevel + (1:floatsPerLevel);
    
    x_float = x(zLevelIndices,:).';
    y_float = y(zLevelIndices,:).';
    
    % The bin edge of 30km is chosen so that the rms separation of the
    % particles matches the dye.
    [D2_level,r2] = PairwiseRelativeDispersion( t, x_float, y_float, [0 30e3 Inf] );
    
    D2_particles = D2_particles + D2_level(:,1);
end

r2_particles = r2(1);

% one factor of 2 to average m_xx and m_yy, another factor to convert
% from relative diffusivity
D2_particles = (D2_particles/nLevels)/4;

[D2_coeff,D2_err] = LinearLeastSquaresFit(t,D2_particles);
kappa_particles = D2_coeff(2)/2;
kappa_err_particles = D2_err(2)/2;

figure
plot(t/86400,D2_particles*1e-6,'LineWidth',2), hold on
plot(t/86400, (D2_coeff(2)*t + D2_coeff(1))*1e-6)