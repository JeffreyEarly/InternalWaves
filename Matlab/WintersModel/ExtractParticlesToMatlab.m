model_dir = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_NL_35e-11_36000s_restart/';
model_dir = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_LIN_unforced_3600000s_restart/';

eulerian_file = [model_dir 'input/SaveIC_EarlyIWmodel.nc'];
lagrangian_dir = [model_dir 'output/lagrangian/'];
floatsPerLevel = 100;

% First we need to read Lx and Ly from the Eulerian files, in order to
% unwrap the Lagrangian trajectories, so that they don't jump across the
% domain.
x_e = ncread(eulerian_file,'x');
Lx = (x_e(end)-x_e(1)) + (x_e(2)-x_e(1)); % They're periodic, so you need an extra dx
y_e = ncread(eulerian_file,'y');
Ly = (y_e(end)-y_e(1)) + (y_e(2)-y_e(1));

UniqueParticleFiles = dir([lagrangian_dir 'particles_*.nc']);

for iFile = 1:length(UniqueParticleFiles)
    file = [lagrangian_dir UniqueParticleFiles(iFile).name];
    if (iFile == 1)
        t = ncread(file,'t_secs');
        x = ncread(file,'x')';
        y = ncread(file,'y')';
        z = ncread(file,'z')';
    else
        x = cat(2,x,ncread(file,'x')');
        y = cat(2,y,ncread(file,'y')');
        z = cat(2,z,ncread(file,'z')');
    end
end

outputfile = [lagrangian_dir 'particles.mat'];
outputfile = '/Users/jearly/Documents/ProjectRepositories/InternalWaves/Matlab/WintersModel/particles_LIN.mat';
save(outputfile,'x','y','z','t','floatsPerLevel', 'model_dir');
