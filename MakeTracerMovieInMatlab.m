file = '/Users/jearly/Desktop/InternalWaves.nc';
FramesFolder ='/Users/jearly/Desktop/InternalWavesUnitTest';

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
% u = double(ncread(file, 'u'));
% v = double(ncread(file, 'v'));
% w = double(ncread(file, 'w'));
% rho = double(ncread(file, 'rho'));

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
floatSize = 25;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Setup the figure
%
figure('Position', [50 50 1920 1080])
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

% Read in the initial position of the floats.
% We will use this information to maintain a constant color on each float.
colormap(jet(128))
coloraxis = linspace(0,1,length(colormap));

xposInitial = double(ncread(file, 'x-position', [ceil(stride/2) ceil(stride/2) ceil(zStride/2) 1], [length(y)/stride length(x)/stride length(z)/zStride 1], [stride stride zStride 1]));
yposInitial = double(ncread(file, 'y-position', [ceil(stride/2) ceil(stride/2) ceil(zStride/2) 1], [length(y)/stride length(x)/stride length(z)/zStride 1], [stride stride zStride 1]));
zposInitial = double(ncread(file, 'z-position', [ceil(stride/2) ceil(stride/2) ceil(zStride/2) 1], [length(y)/stride length(x)/stride length(z)/zStride 1], [stride stride zStride 1]));

xposInitial = reshape(xposInitial, length(x)*length(y)*length(z)/(stride*stride*zStride), 1);
yposInitial = reshape(yposInitial, length(x)*length(y)*length(z)/(stride*stride*zStride), 1);
zposInitial = reshape(zposInitial, length(x)*length(y)*length(z)/(stride*stride*zStride), 1);
	
particle_color = interp1(coloraxis,colormap, (xposInitial-minX)./(maxX-minX) );

particle_size = floatSize*floatSize*ones(size(xposInitial));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% specify the density surface that we want to display, and find their color in the colormap
%
colormap(flipud(jet(128)))
coloraxis = linspace(0,1,length(colormap));
deltaRho = (max(rho_bar)-min(rho_bar))/3;
density_surface = [min(rho_bar) min(rho_bar)+deltaRho min(rho_bar)+2*deltaRho];
density_color = interp1(coloraxis,colormap, (density_surface-min(rho_bar))./(max(rho_bar)-min(rho_bar)) );

%for iTime=1:1%length(t)
for iTime=1:1
	
%	rho3d = double(ncread(file, 'rho', [1 1 1 iTime], [length(yDomain) length(xDomain) length(zDomain) 1], [1 1 1 1]));
% 	rho3d(end+1,:,:) = rho3d(1,:,:);
% 	rho3d(:,end+1,:) = rho3d(:,1,:);

	zeta3d = double(ncread(file, 'zeta', [1 1 1 iTime], [length(yDomain) length(xDomain) length(zDomain) 1], [1 1 1 1]));
	rho3d = zeros(size(zeta3d));
	for m=1:size(rho3d,1)
		for n=1:size(rho3d,2)
			coordinate = squeeze(zeta3d(m,n,:))+zDomain;
			negIndex = find(diff(coordinate)<=0);
			if (length(negIndex) ~= 0)
				for j=(min(negIndex)+1):(max(negIndex)+1)
					if (coordinate(j)<=coordinate(j-1))
						coordinate(j)=coordinate(j-1)+0.001;
					end
				end
			end
			rho3d(m,n,:) = interp1( coordinate, rho_bar, zDomain );
		end
	end
	
	hsurfaces = slice(X,Y,Z,rho3d,[min(xDomain)],[max(yDomain)],[minZ]);
	set(hsurfaces,'FaceColor','interp','EdgeColor','none')
	caxis([min(rho_bar) max(rho_bar)])
		
	view(30,10);

	for i=1:length(density_surface)
		p = patch(isosurface(X,Y,Z,rho3d,density_surface(i)));
		isonormals(X,Y,Z,rho3d,p)
		alpha(p,0.5)
		set(p,'FaceColor',density_color(i,:),'EdgeColor','none');
	end
	
	hold on
	
	% read in the position of the floats for the given time	
	xpos = double(ncread(file, 'x-position', [ceil(stride/2) ceil(stride/2) ceil(zStride/2) iTime], [length(y)/stride length(x)/stride length(z)/zStride 1], [stride stride zStride 1]));
	ypos = double(ncread(file, 'y-position', [ceil(stride/2) ceil(stride/2) ceil(zStride/2) iTime], [length(y)/stride length(x)/stride length(z)/zStride 1], [stride stride zStride 1]));
	zpos = double(ncread(file, 'z-position', [ceil(stride/2) ceil(stride/2) ceil(zStride/2) iTime], [length(y)/stride length(x)/stride length(z)/zStride 1], [stride stride zStride 1]));
	
	% make everything a column vector
	xpos = reshape(xpos, length(x)*length(y)*length(z)/(stride*stride*zStride), 1);
	ypos = reshape(ypos, length(x)*length(y)*length(z)/(stride*stride*zStride), 1);
	zpos = reshape(zpos, length(x)*length(y)*length(z)/(stride*stride*zStride), 1);
	
	xpos = mod( xpos-minX, maxX-minX ) + minX;
	ypos = mod( ypos-minY, maxY-minY ) + minY;
	
	% now plot the floats
	h = scatter3(xpos, ypos, zpos, particle_size, particle_color, 'filled');
	
	xlim([minX maxX])
	ylim([minY maxY])
	%zlim([-60 max(zDomain)])
	zlim([min(zDomain) max(zDomain)])
	
	lighting gouraud
	camlight(30,20)
	camlight
	
	
	hold off
	
	% write everything out	
	output = sprintf('%s/%03d', FramesFolder,iTime-1);
	print('-depsc2', output)
end

error = xpos-xposInitial;
rms_error_x = sqrt(mean(error.*error))

errorZ = zpos-zposInitial;
rms_error_z = sqrt(mean(errorZ.*errorZ))