latitude = 31;
f0 = 2 * 7.2921E-5 * sin( latitude*pi/180 );
N0 = 5.2e-3;
g = 9.81;

xAxisMax = 15*f0;
% xAxisMax = N0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Model dimensions and grid
%

N = 512;
Lmin = 50e3;
aspectRatio = 4;
maxGapTime = 48*3600; % How densely to fill the gaps

Lx = aspectRatio*Lmin;
Ly = Lmin;
Lz = 5000;

Nx = aspectRatio*N;
Ny = N;
Nz = N/2+1;

dx = Lx/Nx;
dy = Ly/Ny;
dz = Lz/Nz;

fprintf('Resolution is %.2fm x %.2fm. There are 2^%d points in the horizontal.\n',dx,dz,round(log(Nx*Ny)/log(2)));

x = dx*(0:Nx-1)'; % periodic basis
y = dy*(0:Ny-1)'; % periodic basis
z = dz*(0:Nz-1)'; % cosine basis (not your usual dct basis, however)

% Spectral domain, in radians
dk = 1/Lx;          % fourier frequency
k = 2*pi*([0:ceil(Nx/2)-1 -floor(Nx/2):-1]*dk)';

dl = 1/Ly;          % fourier frequency
l = 2*pi*([0:ceil(Ny/2)-1 -floor(Ny/2):-1]*dl)';

nModes = Nz;
j = (1:nModes)';

[K,L,J] = ndgrid(k,l,j);
[X,Y,Z] = ndgrid(x,y,z);

M = J*pi/Lz;        % Vertical wavenumber
K2 = K.*K + L.*L;   % Square of the horizontal wavenumber
Kh = sqrt(K2);

C2 = (N0*N0-f0*f0)./(M.*M+K2);
C = sqrt( C2 );                         % Mode speed
h = C2/g;                               % Mode depth
Omega = sqrt(C.*C.*K2 + f0*f0);         % Mode frequency
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Energy spectrum
%
j_star = 15;
L_gm = 1.3e3; % thermocline exponential scale, meters
invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
E_gm = 6.3e-5; % non-dimensional energy parameter
E = L_gm*L_gm*L_gm*invT_gm*invT_gm*E_gm;
N0 = invT_gm;

H = (j_star+(1:3000)).^(-5/2);
H_norm = 1/sum(H);
B_norm = 1/atan(sqrt(N0*N0/(f0*f0)-1));

GM2D_int = @(omega0,omega1,j) E*H_norm*B_norm*((j+j_star).^(-5/2))*(atan(f0/sqrt(omega0*omega0-f0*f0)) - atan(f0/sqrt(omega1*omega1-f0*f0)));
GM2D_uv_int = @(omega0,omega1,j) E*H_norm*B_norm*((j+j_star).^(-5/2))*( f0*sqrt(omega1*omega1-f0*f0)/(2*omega1*omega1) - (3/2)*atan(f0/sqrt(omega1*omega1-f0*f0)) - f0*sqrt(omega0*omega0-f0*f0)/(2*omega0*omega0) + (3/2)*atan(f0/sqrt(omega0*omega0-f0*f0)));
GM2D_w_int = @(omega0,omega1,j) E*H_norm*B_norm*((j+j_star).^(-5/2))*( f0*sqrt(omega1*omega1-f0*f0) + f0*f0*atan(f0/sqrt(omega1*omega1-f0*f0)) - f0*sqrt(omega0*omega0-f0*f0) - f0*f0*atan(f0/sqrt(omega0*omega0-f0*f0)));
GM2D_zeta_int = @(omega0,omega1,j) E*H_norm*B_norm*((j+j_star).^(-5/2))*( ((omega1*omega1-f0*f0)^(3/2))/(2*f0*omega1*omega1) - (1/2)*atan(f0/sqrt(omega1*omega1-f0*f0)) - sqrt(omega1*omega1-f0*f0)/(2*f0) - ((omega0*omega0-f0*f0)^(3/2))/(2*f0*omega0*omega0) + (1/2)*atan(f0/sqrt(omega0*omega0-f0*f0)) + sqrt(omega0*omega0-f0*f0)/(2*f0) );

omegaAxis = linspace(f0,xAxisMax,1000)';
modeAxis = (1:15);

TE = zeros(length(omegaAxis),length(modeAxis));
HKE = zeros(length(omegaAxis),length(modeAxis));
VKE = zeros(length(omegaAxis),length(modeAxis));
IE = zeros(length(omegaAxis),length(modeAxis));
for j = modeAxis
   for iOmega = 1:(length(omegaAxis)-1)
       TE(iOmega,j) = GM2D_int(omegaAxis(iOmega),omegaAxis(iOmega+1),j);
       HKE(iOmega,j) = GM2D_uv_int(omegaAxis(iOmega),omegaAxis(iOmega+1),j);
       VKE(iOmega,j) = GM2D_w_int(omegaAxis(iOmega),omegaAxis(iOmega+1),j);
       IE(iOmega,j) = GM2D_zeta_int(omegaAxis(iOmega),omegaAxis(iOmega+1),j);
   end
end

% Now we create extra points to fill in the gaps
dOmegaInitial = 0.05*f0;
maxdOmega = 2*pi/maxGapTime;
Ln = -1/log(1-dOmegaInitial/maxdOmega);
dOmegas = (1-exp(-(1:100)'/Ln));
gapOmegas = f0 + cumsum(maxdOmega*dOmegas);
omegaExt = [];
jExt = [];
for iMode = 1:nModes
    omegas = sort(reshape(abs(Omega(:,:,iMode)),[],1));
    
    % fill in the lower triangle
    indices = find(gapOmegas < omegas(2));
    jExt = cat(1,jExt,iMode*ones(length(indices),1));
    omegaExt = cat(1,omegaExt,gapOmegas(indices));
    
    diffOmega = diff(omegas);
    gapIndices = find(diffOmega>maxdOmega);
    for i=2:length(gapIndices)
        n = ceil(diffOmega(gapIndices(i))/maxdOmega);
        newOmegas = linspace(omegas(gapIndices(i)),omegas(gapIndices(i)+1),n+1)';
        jExt = cat(1,jExt,iMode*ones(n-1,1));
        omegaExt = cat(1,omegaExt,newOmegas(2:end-1));
    end
end
fprintf('Added %d external waves.\n', length(omegaExt));

% TE = cat(2,TE(:,1),TE);
% HKE = cat(2,HKE(:,1),HKE);
% VKE = cat(2,VKE(:,1),VKE);
% IE = cat(2,IE(:,1),IE);
% modeAxis = cat(2,zeros(1,1),modeAxis);

% This shift is applied to the points along the omega axis for visual aid.
omega_epsilon = 0.0;

ticks = linspace(f0,xAxisMax,5);

labels = cell(length(ticks),1);
labels{1} = 'f_0';
for i=2:(length(ticks)-1)
   labels{i} = sprintf('%df_0',round(ticks(i)/f0));
end
if xAxisMax == N0
    labels{length(ticks)} = 'N_0';
else
    labels{length(ticks)} = sprintf('%df_0',round(xAxisMax/f0));
end

vticks = (1.5:1:max(modeAxis))';
vlabels = cell(length(vticks),1);
for i=1:length(vticks)
   vlabels{i} = sprintf('%d',floor(vticks(i)));
end

scale = @(a) log10(a);

figHandle = figure('Position',[100 100 700 600]);

sp1 = subplot(2,2,1);
pcolor(omegaAxis,modeAxis,scale(TE/max(max(TE)))'), hold on
caxis([-2 0])
shading flat
for j = modeAxis
    omega = sort(reshape(Omega(:,:,j),1,[]));
    scatter(omega+omega_epsilon,j*ones(size(omega))+0.5,16*ones(size(omega)),'filled', 'MarkerFaceColor', 0*[1 1 1])
end
scatter(omegaExt+omega_epsilon,jExt+0.5,16*ones(size(omegaExt)),'filled', 'MarkerFaceColor', 1*[1 1 1])
title('total')
ylabel('vertical mode')
% xlabel('frequency')
xticks([])
% xticklabels(labels)
yticks(vticks)
yticklabels(vlabels)

sp2 = subplot(2,2,2);
pcolor(omegaAxis,modeAxis,scale(HKE/max(max(HKE)))'), hold on
caxis([-2 0])
shading flat
for j = modeAxis
    omega = sort(reshape(Omega(:,:,j),1,[]));
    scatter(omega+omega_epsilon,j*ones(size(omega))+0.5,16*ones(size(omega)),'filled', 'MarkerFaceColor', 0*[1 1 1])
end
scatter(omegaExt+omega_epsilon,jExt+0.5,16*ones(size(omegaExt)),'filled', 'MarkerFaceColor', 1*[1 1 1])
title('horizontal')
% ylabel('vertical mode')
yticks([])
% xlabel('frequency')
xticks([])
% xticklabels(labels)

sp3 = subplot(2,2,3);
pcolor(omegaAxis,modeAxis,scale(VKE/max(max(VKE)))'), hold on
caxis([-2 0])
shading flat
for j = modeAxis
    omega = sort(reshape(Omega(:,:,j),1,[]));
    scatter(omega+omega_epsilon,j*ones(size(omega))+0.5,16*ones(size(omega)),'filled', 'MarkerFaceColor', 0*[1 1 1])
end
scatter(omegaExt+omega_epsilon,jExt+0.5,16*ones(size(omegaExt)),'filled', 'MarkerFaceColor', 1*[1 1 1])
title('vertical')
ylabel('vertical mode')
xlabel('frequency')
xticks(ticks)
xticklabels(labels)
yticks(vticks)
yticklabels(vlabels)

sp4 = subplot(2,2,4);
pcolor(omegaAxis,modeAxis,scale(IE/max(max(IE)))'), hold on
caxis([-2 0])
shading flat
for j = modeAxis
    omega = sort(reshape(Omega(:,:,j),1,[]));
    scatter(omega+omega_epsilon,j*ones(size(omega))+0.5,16*ones(size(omega)),'filled', 'MarkerFaceColor', 0*[1 1 1])
end
scatter(omegaExt+omega_epsilon,jExt+0.5,16*ones(size(omegaExt)),'filled', 'MarkerFaceColor', 1*[1 1 1])
title('isopycnal')
% ylabel('vertical mode')
yticks([])
xlabel('frequency')
xticks(ticks)
xticklabels(labels)
% c = colorbar;

dy = 0.045;
sp1.Position(2) = sp1.Position(2) - dy;
sp1.Position(4) = sp1.Position(4) + dy;
sp3.Position(4) = sp3.Position(4) + dy;
sp2.Position(2) = sp2.Position(2) - dy;
sp2.Position(4) = sp2.Position(4) + dy;
sp4.Position(4) = sp4.Position(4) + dy;

dx = 0.045;
sp1.Position(3) = sp1.Position(3) + dx;
sp3.Position(3) = sp3.Position(3) + dx;
sp2.Position(1) = sp2.Position(1) - dx;
sp2.Position(3) = sp2.Position(3) + dx;
sp4.Position(1) = sp4.Position(1) - dx;
sp4.Position(3) = sp4.Position(3) + dx;

figHandle.NextPlot = 'add';
a = axes; 

%// Set the title and get the handle to it
ht = title(sprintf('%dkm x %dkm x %dkm (%dx%dx%d)',round(Lx/1e3),round(Ly/1e3),round(Lz/1e3),Nx,Ny,Nz),'FontSize', 24);
ht.Position = [0.5 1.04 0.5];

%// Turn the visibility of the axes off
a.Visible = 'off';

%// Turn the visibility of the title on
ht.Visible = 'on';

% if xAxisMax == N0
%     print('ResolutionAndVariance.png','-r200','-dpng')
% else
%     print('ResolutionAndVarianceZoomed.png','-r200','-dpng')
% end
