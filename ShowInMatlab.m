file = '/Users/jearly/Desktop/InternalWaves.nc';

% read in the dimensional variables
x = ncread(file, 'x');
y = ncread(file, 'y');
z = ncread(file, 'z');
t = ncread(file, 'time');

% read in the mean density profile
rho_bar = double(ncread(file, 'rho_bar'));

% read in the dynamical variables
u = double(ncread(file, 'u'));
v = double(ncread(file, 'v'));
w = double(ncread(file, 'w'));
rho = double(ncread(file, 'rho'));


% [X,Y,Z]=meshgrid(x,y,z);
% figure
% colormap(flipud(jet(128)))
% hsurfaces = slice(X,Y,Z,rho3d,[min(x)],[max(y)],[min(z)])
% set(hsurfaces,'FaceColor','interp','EdgeColor','none')
% caxis([min(rho_bar) max(rho_bar)])
% view(30,10);
% 
% deltaRho = (max(rho_bar)-min(rho_bar))/3;
% rho1 = min(rho_bar)+deltaRho;
% rho2 = min(rho_bar)+2*deltaRho;
% coloraxis = linspace(0,1,length(colormap));
% rho1Color = interp1(coloraxis,colormap, 1/3);
% rho2Color = interp1(coloraxis,colormap, 2/3);
% 
% p = patch(isosurface(X,Y,Z,rho3d,rho1));
% %isonormals(X,Y,Z,rho3d,p)
% alpha(p,0.5)
% set(p,'FaceColor',rho1Color,'EdgeColor','none');
% 
% p = patch(isosurface(X,Y,Z,rho3d,rho2));
% %isonormals(X,Y,Z,rho3d,p)
% alpha(p,0.5)
% set(p,'FaceColor',rho2Color,'EdgeColor','none');
% 
% lighting gouraud
% camlight(30,20)
% camlight
% return

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(4) scrsz(4)/2])
for iTime=1:4
	u3d = u(:,:,:,iTime);
	v3d = v(:,:,:,iTime);
	w3d = w(:,:,:,iTime);
	rho3d = rho(:,:,:,iTime);
	
	subplot(3,4,0*4+iTime)
	
	zIndex = 50;
	u2d = squeeze(u3d(:,:,zIndex));
	v2d = squeeze(v3d(:,:,zIndex));
	quiver(x,y,u2d',v2d',0.8)
	xlim([min(x) max(x)])
	ylim([min(y) max(y)])
	
	subplot(3,4,1*4+iTime)
	pcolor(x,z,squeeze(w3d(:,1,:))'),shading flat

	subplot(3,4,2*4+iTime)
	pcolor(x,z,squeeze(rho3d(:,1,:))'),shading flat
end	


[X,Y,Z]=meshgrid(x,y,z);
p = patch(isosurface(X,Y,Z,rho3d,1025.5));
isonormals(X,Y,Z,rho3d,p)
set(p,'FaceColor','red','EdgeColor','none');
% %daspect([1,1,1])
% view(3); axis tight
% camlight 
% lighting gouraud
% xlim([min(x) max(x)])
% ylim([min(y) max(y)])
% zlim([min(z) max(z)])