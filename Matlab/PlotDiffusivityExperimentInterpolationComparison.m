files{1} = '/Volumes/MoreStorage/DiffusivityExperiment_2020-01-24T045517_128x128x129.nc';
files{2} = '/Volumes/MoreStorage/DiffusivityExperiment_2020-01-24T110233_64x64x129.nc';
files{3} = '/Volumes/MoreStorage/DiffusivityExperiment_2020-01-24T123859_32x32x129.nc';
files{4} = '/Volumes/MoreStorage/DiffusivityExperiment_2020-01-24T133617_64x64x129.nc';
files{5} = '/Volumes/MoreStorage/DiffusivityExperiment_2020-01-25T193808_256x256x129.nc';

nFiles = 5;

file = files{1};
x = ncread(file,'x-position').';
y = ncread(file,'y-position').';
z = ncread(file,'z-position').';
t = ncread(file,'t');

nLevels = 5;
floatsPerLevel = size(x,2)/nLevels;

zLevel = 3;
zLevelIndices = (zLevel-1)*floatsPerLevel + (1:floatsPerLevel);
edges = CreateBinEdgesForInitialSeparation(t,x(:,zLevelIndices),y(:,zLevelIndices));

figure('Name','LateralDiffusivity with resolution changes')
for iFile=1:nFiles
    x = ncread(files{iFile},'x-position').';
    y = ncread(files{iFile},'y-position').';
    t = ncread(files{iFile},'t');
    
    if iFile == 5
       validIndices = 1:length(t)-1;
       t=t(validIndices);
       x=x(validIndices,:);
       y=y(validIndices,:);
    end
    
    floatsPerLevel = size(x,2)/nLevels;
    zLevel = 3;
    zLevelIndices = (zLevel-1)*floatsPerLevel + (1:floatsPerLevel);
    
    x_float = x(:,zLevelIndices);
    y_float = y(:,zLevelIndices);
    [r2, kappa_r, kappa_err ] = PairwiseRelativeDiffusivityFromSlope(t, x(:,zLevelIndices), y(:,zLevelIndices), edges );
    errorbar(sqrt(r2)/1000,kappa_r, 2*kappa_err), hold on
end
xlabel('distance (m)')
ylabel('kappa (m^2/s)')
legend('128x128x129','64x64x129','32x32x129','64x64x129 @ L/2','256x256x129','Location','northwest')
print('-depsc', 'LateralDiffusivity-vs-Resolution.eps')

return

t_particles = cell(nFiles,1);
D2_particles = cell(nFiles,1);
r2_particles = zeros(nFiles,1);
for iFile=5:nFiles
   x = ncread(files{iFile},'x-position').';
   y = ncread(files{iFile},'y-position').';
   t = ncread(files{iFile},'t');
   
   t_particles{iFile} = t;
   D2_particles{iFile} = zeros(length(t),length(edges));
   
   floatsPerLevel = size(x,2)/nLevels;
   for zLevel = 1:nLevels
       nLevels = size(x,2)/floatsPerLevel;
       zLevelIndices = (zLevel-1)*floatsPerLevel + (1:floatsPerLevel);
       
       x_float = x(:,zLevelIndices);
       y_float = y(:,zLevelIndices);
       
       % The bin edge of 30km is chosen so that the rms separation of the
       % particles matches the dye.
       [D2_level,r2] = PairwiseRelativeDispersion( t, x_float, y_float, edges );
       
       D2_particles{iFile} = D2_particles{iFile} + D2_level;
   end
   
   r2_particles(iFile) = r2(1);
   
   % one factor of 2 to average m_xx and m_yy, another factor to convert
   % from relative diffusivity
   D2_particles{iFile} = (D2_particles{iFile}/nLevels)/4;
   

end

figure('Name','LateralDispersion-vs-Energy')
for iFile=1:nFiles
   subplot(nFiles,1,iFile)
   plot(t_particles{nFiles-iFile+1}/86400,D2_particles{nFiles-iFile+1}*1e-6,'LineWidth',2), hold on
   plot(t_tracer{nFiles-iFile+1}/86400,D2_tracer{nFiles-iFile+1}*1e-6,'LineWidth',2)
   if iFile==1
      title(sprintf('Lateral Dispersion at (%d km)^2',round(sqrt(mean(r2_particles))*1e-3)))
      legend('particles','tracer','Location','northwest') 
   end
   ylabel('km^2')
end
xlabel('days')
packfig(nFiles,1)

return

[r2, kappa_r, kappa_err ] = PairwiseRelativeDiffusivityFromSlope(t, x(:,zLevelIndices), y(:,zLevelIndices), [0 30e3 Inf] );
[D2,r2,r0] = PairwiseRelativeDispersion( t, x(:,zLevelIndices), y(:,zLevelIndices), [0 30e3 Inf] );

figure
plot(t/86400,D2*1e-6)

figure
errorbar(sqrt(r2)/1000,kappa_r, 2*kappa_err)
ylim([0 max(kappa_r+2*kappa_err)])
xlabel('distance (m)')
ylabel('kappa (m^2/s)')
