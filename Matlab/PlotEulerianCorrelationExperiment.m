
% This case has k = nyquist
file = '/Volumes/OceanTransfer/EulerianCorrelationExperiment_2017-05-13T135843_128x128x128.nc';

% This case has k = nyquist/4
file = '/Volumes/OceanTransfer/EulerianCorrelationExperiment_2017-05-13T141304_128x128x128.nc';

% % This case has a full wave field, with domain size L=15km
% file = '/Volumes/OceanTransfer/EulerianCorrelationExperiment_2017-05-13T143411_128x128x128.nc';
% 
% % This case has a full wave field, with domain size L=150km
% file = '/Volumes/OceanTransfer/EulerianCorrelationExperiment_2017-05-13T143903_128x128x128.nc';
% 
% % This case has k_cutof = 4*dk with domain size L=150km
file = '/Volumes/OceanTransfer/EulerianCorrelationExperiment_2017-05-16T104551_128x128x128.nc';

% This case has k_cutof = 8*dk, \pm dk (instead of dk/2) with domain size L=150km
% file = '/Volumes/OceanTransfer/EulerianCorrelationExperiment_2017-05-16T135758_128x128x128.nc';
% 
% % This case has k_cutof = 12*dk, \pm dk (instead of dk/2) with domain size L=150km
% file = '/Volumes/OceanTransfer/EulerianCorrelationExperiment_2017-05-16T140248_128x128x128.nc';
% 
% % This case has k_cutof = 4*dk with domain size L=150km, but this time \pm
% % dk instead of dk/2.
% file = '/Volumes/OceanTransfer/EulerianCorrelationExperiment_2017-05-16T142836_128x128x128.nc';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read in the problem dimensions and parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = ncread(file, 'x');
y = ncread(file, 'y');
z = ncread(file, 'z');
t = ncread(file, 't');

% u1 = double(squeeze(ncread(file, 'u', [1 1 depth_index 1], [length(x) length(y) 1 1], [1 1 1 1])));
% u10 = double(squeeze(ncread(file, 'u', [1 1 depth_index 10], [length(x) length(y) 1 1], [1 1 1 1])));

depth_index = 1;


f = @(u0,un) mean((u0 - mean(u0)).*(un - mean(un)))/(std(u0,1)*std(un,1));

stride = 4;
r = zeros(size(x));
for m=1:stride:length(y)
    u0 = double(squeeze(ncread(file, 'u', [m 1 depth_index 1], [1 1 1 length(t)], [1 1 1 1])));
    v0 = double(squeeze(ncread(file, 'v', [m 1 depth_index 1], [1 1 1 length(t)], [1 1 1 1])));
    for n=1:length(x)
        un = double(squeeze(ncread(file, 'u', [m n depth_index 1], [1 1 1 length(t)], [1 1 1 1])));
        vn = double(squeeze(ncread(file, 'v', [m n depth_index 1], [1 1 1 length(t)], [1 1 1 1])));
        
        %     corrLength=length(u0)+length(un)-1;
        %     N = ceil(corrLength/2);
        %     unbiased = [1:N N-(1:(N-1))]';
        %     c=fftshift(ifft(fft(u0,corrLength).*conj(fft(un,corrLength))))/(mean(u0.*un))./unbiased; % std(u0,1)*std(un,1)
        
        %     r(n) = c(N);
        r(n) = r(n) + 0.5*(f(u0,un) + f(v0,vn));
    end
end
r = r/length(1:stride:length(y));

figure, plot(x,r)

% This experiment does suggest that you can only use half the domain width
% before things start correlating again.