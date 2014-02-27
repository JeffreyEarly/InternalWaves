file = '/Users/jearly/Desktop/InternalWavesLatmix2011_128_128_64_lat31_unit_test.nc';
FramesFolder ='/Users/jearly/Desktop/InternalWavesLatmix2011_128_128_64_lat31_unit_testFrames';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Make the frames folder
%
if exist(FramesFolder) == 0
	mkdir(FramesFolder);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Read in the problem dimensions
%
xDomain = ncread(file, 'x');
yDomain = ncread(file, 'y');
zDomain = ncread(file, 'z');
[X,Y,Z]=meshgrid(xDomain,yDomain,zDomain);

x = ncread(file, 'x-float');
y = ncread(file, 'y-float');
z = ncread(file, 'z-float');
t = ncread(file, 'time');

% read in the mean density profile
rho_bar = double(ncread(file, 'rho_bar'));

% read in the dynamical variables
u = double(ncread(file, 'u'));
v = double(ncread(file, 'v'));
w = double(ncread(file, 'w'));
rho = double(ncread(file, 'rho'));

deltaX = xDomain(2)-xDomain(1);
minX = min(xDomain);
maxX = max(xDomain+deltaX);

deltaY = yDomain(2)-yDomain(1);
minY = min(yDomain);
maxY = max(yDomain+deltaY);

minZ = min(zDomain);
maxZ = max(zDomain);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	The stride indicates how many floats we will skip
%
stride = 1;
zStride = 1;
floatSize = 5;

% Read in the initial position of the floats.
% We will use this information to maintain a constant color on each float.
xposInitial = double(ncread(file, 'x-position', [ceil(stride/2) ceil(stride/2) 1 1], [length(x)/stride length(y)/stride length(z)/zStride 1], [stride stride zStride 1]));
xposInitial = reshape(xposInitial, length(x)*length(y)*length(z)/(stride*stride*zStride), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Setup the figure
%
figure('Position', [50 50 1920 1080])
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

%for iTime=1:length(t)
for iTime=1:1

	% read in the position of the floats for the given time	
	xpos = double(ncread(file, 'x-position', [ceil(stride/2) ceil(stride/2) ceil(zStride/2) iTime], [length(x)/stride length(y)/stride length(z)/zStride 1], [stride stride zStride 1]));
	ypos = double(ncread(file, 'y-position', [ceil(stride/2) ceil(stride/2) ceil(zStride/2) iTime], [length(x)/stride length(y)/stride length(z)/zStride 1], [stride stride zStride 1]));
	zpos = double(ncread(file, 'z-position', [ceil(stride/2) ceil(stride/2) ceil(zStride/2) iTime], [length(x)/stride length(y)/stride length(z)/zStride 1], [stride stride zStride 1]));
	
	% make everything a column vector
	xpos = reshape(xpos, length(x)*length(y)*length(z)/(stride*stride*zStride), 1);
	ypos = reshape(ypos, length(x)*length(y)*length(z)/(stride*stride*zStride), 1);
	zpos = reshape(zpos, length(x)*length(y)*length(z)/(stride*stride*zStride), 1);
	
	xpos = mod( xpos-minX, maxX-minX ) + minX;
	ypos = mod( ypos-minY, maxY-minY ) + minY;
	
	% default color map is only 128 shades---we need more!
	colormap(jet(1024))	
	
	% now plot the floats, colored by initial position
	% Scatter works, but is substantially slower than using mesh.
	% scatter(xpos, ypos, floatSize*floatSize, xposInitial, 'filled')	
	mesh([xpos';xpos'],[ypos';ypos'],[xposInitial';xposInitial'],'mesh','column','marker','.','MarkerSize',floatSize*floatSize), view(2)
	grid off
	
	% make the axes look better
	set( gca, 'TickDir', 'out');
	set( gca, 'Linewidth', 1.0);
	axis equal tight
	
	% get rid of the xticks because we're going to put a colorbar below with the same info.
	set( gca, 'xtick', [])
	
	xlim([minX maxX])
	ylim([minY maxY])
	
	% label everything
	title( sprintf('Floats advected by a Poincare Wave field, colored by initial position, hour %.1f', t(iTime)/3600), 'fontsize', 28, 'FontName', 'Helvetica' );
	ylabel( 'distance', 'FontSize', 24.0, 'FontName', 'Helvetica');
	
	% add a color bar
	cb = colorbar( 'location', 'SouthOutside' );
	set(get(cb,'xlabel'),'String', 'distance', 'FontSize', 24.0, 'FontName', 'Helvetica');
	set( gca, 'clim', [minX maxX] );
	
	% write everything out	
% 	output = sprintf('%s/%03d', FramesFolder,iTime-1);
% 	print('-depsc2', output)
end