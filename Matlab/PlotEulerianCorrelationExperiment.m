
% k=8*dk
file = '/Volumes/OceanTransfer/EulerianCorrelationExperiment_2017-05-17T111803_128x128x128.nc';

% k=6*dk
% file = '/Volumes/OceanTransfer/EulerianCorrelationExperiment_2017-05-17T112038_128x128x128.nc';

% k=4*dk
% file = '/Volumes/OceanTransfer/EulerianCorrelationExperiment_2017-05-17T112306_128x128x128.nc';

% k=3*dk
% file = '/Volumes/OceanTransfer/EulerianCorrelationExperiment_2017-05-17T112536_128x128x128.nc';

% k=2*dk
% file = '/Volumes/OceanTransfer/EulerianCorrelationExperiment_2017-05-17T112849_128x128x128.nc';

% k=1*dk
% file = '/Volumes/OceanTransfer/EulerianCorrelationExperiment_2017-05-17T113129_128x128x128.nc';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read in the problem dimensions and parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = ncread(file, 'x');
y = ncread(file, 'y');
z = ncread(file, 'z');
t = ncread(file, 't');

maxWavelength = ncreadatt(file,'/','max-wavelength-in-spectrum');
maxPeriod = ncreadatt(file,'/','max-period-in-spectrum');

titlestring = sprintf('Largest wave, L=%.2f km, T=%.1f min',maxWavelength/1e3, maxPeriod/60);

depth_index = 1;


f = @(u0,un) mean((u0 - mean(u0)).*(un - mean(un)))/(std(u0,1)*std(un,1));

stride = 1;
r = zeros(size(x));
AC = zeros(size(t));
DOF = zeros(size(t));
yLoop = 1:stride:length(y);
for m=yLoop
    u0 = double(squeeze(ncread(file, 'u', [m 1 depth_index 1], [1 1 1 length(t)], [1 1 1 1])));
    v0 = double(squeeze(ncread(file, 'v', [m 1 depth_index 1], [1 1 1 length(t)], [1 1 1 1])));
    for n=1:length(x)
        un = double(squeeze(ncread(file, 'u', [m n depth_index 1], [1 1 1 length(t)], [1 1 1 1])));
        vn = double(squeeze(ncread(file, 'v', [m n depth_index 1], [1 1 1 length(t)], [1 1 1 1])));
        
        r(n) = r(n) + 0.5*(f(u0,un) + f(v0,vn));
    end
    
    % Compute the autocorrelation of the time series
    [ACu, DOFu] = Autocorrelation(u0, length(t)-1);
    AC = AC + ACu;
    DOF = DOF + DOFu;
end
r = r/length(yLoop);
AC = AC/length(yLoop);

SE_indep = t(2:end);
SE =  sqrt((1 + 2*cumsum(AC.^2))./DOF);
SE(1) = sqrt(1/DOF(1)); % first point is lag 1
SE(end) = []; % there is no end point

s = 1e-3;

figure('Position', [50 50 1000 500] )
subplot(1,2,1)
plot(s*x,r)
vlines(s*[maxWavelength max(x)-maxWavelength],'g--')
xlabel('distance (km)')
ylabel('spatial correlation')
title(titlestring)

s = 1/60;
maxT = max(t);

subplot(1,2,2)
plot(s*t,AC, 'LineWidth',1,'Color',0.0*[1.0 1.0 1.0])
hold on
plot(s*SE_indep, [3*SE,-3*SE], 'LineWidth', 1.5, 'Color',0.4*[1.0 1.0 1.0] )
vlines(s*maxPeriod,'g--')
xlabel('time lag (minutes)')
ylabel('temporal correlation')
xlim([0 s*maxT])
ylim([-0.2 1.0])

print('-depsc2', sprintf('L=%d_km.eps',round(maxWavelength/1e3)))

% This experiment does suggest that you can only use half the domain width
% before things start correlating again.