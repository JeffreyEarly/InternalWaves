file = '/Users/jearly/Desktop/InternalWaves.nc';
x = ncread(file, 'x');
y = ncread(file, 'y');
z = ncread(file, 'z');
t = ncread(file, 'time');
u = double(ncread(file, 'u'));
v = double(ncread(file, 'v'));
w = double(ncread(file, 'w'));

figure, pcolor(x,z,squeeze(w(:,1,:))'),shading flat
% 
% figure, pcolor(x,y,squeeze(u(:,1,:))'),shading flat

figure

zIndex = 50;
u2d = squeeze(u(:,:,zIndex));
v2d = squeeze(v(:,:,zIndex));
quiver(x,y,u2d',v2d',0.8)
xlim([min(x) max(x)])
ylim([min(y) max(y)])

% scrsz = get(0,'ScreenSize');
% figure('Position',[1 scrsz(4)/2 scrsz(4) scrsz(4)/2])
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