%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SaveModelOutputToNetCDF
%
% This script uses the InternalWaveModel to generate the time-evolution of
% a Garrett-Munk spectrum of internal waves and save the output to a NetCDF
% file.
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com
%
% December 6th, 2016      Version 1.0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx = 30e3;
Ly = 15e3;
Lz = 5000;

Nx = 512;
Ny = 256;
Nz = 64;

% Lx = 15e3;
% Ly = 15e3;
% Lz = 5000;
% 
% Nx = 32;
% Ny = 32;
% Nz = 16;

latitude = 31;
N0 = 5.2e-3;
GMReferenceLevel = 1.0;

timeStep = 15*60; % in seconds
maxTime = 4*86400;

outputfolder = '/Volumes/OceanTransfer';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavemodel = InternalWaveModel([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);
wavemodel.InitializeWithGMSpectrum(GMReferenceLevel);

t = (0:timeStep:maxTime)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create a NetCDF file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filepath = sprintf('%s/InternalWave5.nc', outputfolder);

% Apple uses 1e9 bytes as 1 GB (not the usual multiples of 2 definition)
totalFields = 4;
bytePerDouble = 8;
totalSize = totalFields*bytePerDouble*length(t)*(wavemodel.Nx)*(wavemodel.Ny)*(wavemodel.Nz)/1e9;
fprintf('Writing output file to %s\nExpected file size is %.2f GB.\n',filepath,totalSize);

ncid = netcdf.create(filepath, 'CLOBBER');

% Define the dimensions
xDimID = netcdf.defDim(ncid, 'x', wavemodel.Nx);
yDimID = netcdf.defDim(ncid, 'y', wavemodel.Ny);
zDimID = netcdf.defDim(ncid, 'z', wavemodel.Nz);
tDimID = netcdf.defDim(ncid, 't', netcdf.getConstant('NC_UNLIMITED'));

% Define the coordinate variables
xVarID = netcdf.defVar(ncid, 'x', 'double', xDimID);
yVarID = netcdf.defVar(ncid, 'y', 'double', yDimID);
zVarID = netcdf.defVar(ncid, 'z', 'double', zDimID);
tVarID = netcdf.defVar(ncid, 't', 'double', tDimID);
netcdf.putAtt(ncid,xVarID, 'units', 'm');
netcdf.putAtt(ncid,yVarID, 'units', 'm');
netcdf.putAtt(ncid,zVarID, 'units', 'm');
netcdf.putAtt(ncid,tVarID, 'units', 's');

% Define the dynamical variables
uVarID = netcdf.defVar(ncid, 'u', 'double', [xDimID,yDimID,zDimID,tDimID]);
vVarID = netcdf.defVar(ncid, 'v', 'double', [xDimID,yDimID,zDimID,tDimID]);
wVarID = netcdf.defVar(ncid, 'w', 'double', [xDimID,yDimID,zDimID,tDimID]);
zetaVarID = netcdf.defVar(ncid, 'zeta', 'double', [xDimID,yDimID,zDimID,tDimID]);
netcdf.putAtt(ncid,uVarID, 'units', 'm/s');
netcdf.putAtt(ncid,vVarID, 'units', 'm/s');
netcdf.putAtt(ncid,wVarID, 'units', 'm/s');
netcdf.putAtt(ncid,zetaVarID, 'units', 'm');

% Write some metadata
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'latitude', latitude);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'N0', N0);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'GMReferenceLevel', GMReferenceLevel);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'Model', 'Created from InternalWaveModel.m written by Jeffrey J. Early.');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'ModelVersion', wavemodel.version);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'CreationDate', datestr(datetime('now')));

% End definition mode
netcdf.endDef(ncid);

% Add the data for the coordinate variables
netcdf.putVar(ncid, xVarID, wavemodel.x);
netcdf.putVar(ncid, yVarID, wavemodel.y);
netcdf.putVar(ncid, zVarID, wavemodel.z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Run the model, and write the output to NetCDF
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iTime=1:length(t)
    if mod(iTime,10) == 0
        fprintf('writing values time step %d of %d to file.\n', iTime, length(t));
    end
    [u,v]=wavemodel.VelocityFieldAtTime(t(iTime));
    [w,zeta] = wavemodel.VerticalFieldsAtTime(t(iTime));
    netcdf.putVar(ncid, uVarID, [0 0 0 iTime-1], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], u);
    netcdf.putVar(ncid, vVarID, [0 0 0 iTime-1], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], v);
    netcdf.putVar(ncid, wVarID, [0 0 0 iTime-1], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], w);
    netcdf.putVar(ncid, zetaVarID, [0 0 0 iTime-1], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], zeta);
    netcdf.putVar(ncid, tVarID, iTime-1, 1, t(iTime));
end

netcdf.close(ncid);	
