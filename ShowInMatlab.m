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


% 
% subplot(2,2,1)
% pcolor(x, y, ssh(:,:,1)), axis equal tight, shading interp
% title('SSH of gaussian eddy, day 0')
% subplot(2,2,3)
% pcolor(x, y, ssh(:,:,end)), axis equal tight, shading interp
% title(sprintf('SSH of gaussian eddy, day %d', round(t(end))))
% 
% 
% subplot(2,2,2)
% pcolor(x, y, tracer(:,:,1)), axis equal tight, shading interp
% title('Tracer advected by a gaussian eddy, day 0')
% subplot(2,2,4)
% pcolor(x, y, tracer(:,:,end)), axis equal tight, shading interp
% title(sprintf('Tracer advected by a gaussian eddy, day %d', round(t(end))))