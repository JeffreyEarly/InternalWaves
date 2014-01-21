file = '/Users/jearly/Desktop/InternalWaves.nc';
FramesFolder ='/Users/jearly/Desktop/InternalWavesVelocityfieldFrames';

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
x = ncread(file, 'x');
y = ncread(file, 'y');
t = ncread(file, 'time');
u = double(ncread(file, 'u'));
v = double(ncread(file, 'v'));
xpos = ncread(file, 'x-position');
ypos = ncread(file, 'y-position');
zpos = ncread(file, 'z-position');

minX = min(x);
maxX = max(x)+x(2)-x(1);
minY = min(y);
maxY = max(y)+y(2)-y(1);

xpos = mod( xpos-minX, maxX-minX ) + minX;
ypos = mod( ypos-minY, maxY-minY ) + minY;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Setup the figure
%
figure('Position', [50 50 1920 1080])
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

i=floor(size(xpos,1)/2);
j=floor(size(ypos,2)/2);
k=size(ypos,3);

for iTime=1:length(t)

	quiver(x,y,u(:,:,end,iTime)',v(:,:,end,iTime)',0.8)
	xlim([min(x) max(x)])
	ylim([min(y) max(y)])
	
	hold on
	scatter(squeeze(xpos(i,j,k,iTime)), squeeze(ypos(i,j,k,iTime)),10^2,'MarkerEdgeColor','b','MarkerFaceColor','c','LineWidth',1.5)
	hold off
	
	% write everything out	
	output = sprintf('%s/Day_%03d', FramesFolder,iTime-1);
	print('-depsc2', output)
end