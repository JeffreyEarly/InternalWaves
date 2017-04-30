%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DiffusivityExperiment
%
% This script uses the InternalWaveModel to generate the time-evolution of
% a Garrett-Munk spectrum of internal waves and save the output to a NetCDF
% file. Advects particles as well.
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com
%
% April 12th, 2017      Version 1.0

N = 64;
aspectRatio = 8;

L = 15e3;
Lx = aspectRatio*L;
Ly = L;
Lz = 5000;

Nx = aspectRatio*N;
Ny = N;
Nz = (N/2)+1; % Must include end point to advect at the surface, so use 2^N + 1

latitude = 31;
N0 = 5.2e-3; % Choose your stratification
GMReferenceLevel = 1.0;

kappa = 5e-6;
outputInterval = 15*60;
maxTime = 86400/4;

outputfolder = '/Volumes/OceanTransfer';

precision = 'single';

if strcmp(precision,'single')
    ncPrecision = 'NC_FLOAT';
    setprecision = @(x) single(x);
    bytePerFloat = 4;
else
    ncPrecision = 'NC_DOUBLE';
    setprecision = @(x) double(x);
    bytePerFloat = 8;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shouldUseGMSpectrum = 1;

wavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

if shouldUseGMSpectrum == 1
    wavemodel.FillOutWaveSpectrum();
    wavemodel.InitializeWithGMSpectrum(GMReferenceLevel);
    wavemodel.ShowDiagnostics();
    period = 2*pi/wavemodel.N0;
    [u,v] = wavemodel.VelocityFieldAtTime(0.0);
    U = max(max(max( sqrt(u.*u + v.*v) )));
else
    j0 = 1; % j=1..nModes, where 1 indicates the 1st baroclinic mode
    U = 0.1; % m/s
    sign = 1;
    phi = 0;
    k0 = 2;
    l0 = 0;
    alpha = atan2(l0,k0);
    k = 2*pi*sqrt(k0^2 + l0^2)/Lx;
    
    period = wavemodel.InitializeWithPlaneWave(k0,l0,j0,U,sign);
    maxTime = period;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create floats/drifters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx = wavemodel.x(2)-wavemodel.x(1);
dy = wavemodel.y(2)-wavemodel.y(1);
N = 10;
x_float = (1:N)*dx;
y_float = (1:N)*dy;
z_float = (0:2)*(-Lz/4);

[x_float,y_float,z_float] = ndgrid(x_float,y_float,z_float);
x_float = reshape(x_float,[],1);
y_float = reshape(y_float,[],1);
z_float = reshape(z_float,[],1);
nFloats = numel(x_float);

% Now let's place the floats along an isopycnal.
[w,zeta] = wavemodel.VerticalFieldsAtTime(0);
isopycnalDeviation = interpn(wavemodel.X,wavemodel.Y,wavemodel.Z,zeta,x_float,y_float,z_float,'spline');
z_isopycnal = z_float + isopycnalDeviation;

ymin = [-Inf -Inf -Lz -Inf -Inf -Lz -Inf -Inf];
ymax = [Inf Inf 0 Inf Inf 0 Inf Inf];
kappa_vector = [0 0 0 kappa kappa kappa 0 0];
p0 = cat(2, x_float, y_float, z_isopycnal, x_float, y_float, z_isopycnal, x_float, y_float);

f = @(t,y) FluxForFloatDiffusiveDrifter(t,y,z_float,wavemodel, 'exact');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Determine the proper time interval
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfl = 0.25;
advectiveDT = cfl*(wavemodel.x(2)-wavemodel.x(1))/U;
oscillatoryDT = period/8;
if advectiveDT < oscillatoryDT
    fprintf('Using the advective dt: %.2f\n',advectiveDT);
    deltaT = advectiveDT;
else
    fprintf('Using the oscillatory dt: %.2f\n',oscillatoryDT);
    deltaT = oscillatoryDT;
end

deltaT = outputInterval/ceil(outputInterval/deltaT);
fprintf('Rounding to match the output interval dt: %.2f\n',deltaT);

t = (0:deltaT:maxTime)';
if t(end) < period
    t(end+1) = period;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create a NetCDF file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filepath = sprintf('%s/InternalWaveModel_%s_%dx%dx%d.nc', outputfolder,datestr(datetime('now'),'yyyy-mm-ddTHHMMSS'),Nx,Ny,Nz);

% Apple uses 1e9 bytes as 1 GB (not the usual multiples of 2 definition)
totalFields = 4;
totalSize = totalFields*bytePerFloat*length(t)*(wavemodel.Nx)*(wavemodel.Ny)*(wavemodel.Nz)/1e9;
fprintf('Writing output file to %s\nExpected file size is %.2f GB.\n',filepath,totalSize);

ncid = netcdf.create(filepath, 'CLOBBER');

% Define the dimensions
xDimID = netcdf.defDim(ncid, 'x', wavemodel.Nx);
yDimID = netcdf.defDim(ncid, 'y', wavemodel.Ny);
zDimID = netcdf.defDim(ncid, 'z', wavemodel.Nz);
tDimID = netcdf.defDim(ncid, 't', netcdf.getConstant('NC_UNLIMITED'));

% Define the coordinate variables
xVarID = netcdf.defVar(ncid, 'x', ncPrecision, xDimID);
yVarID = netcdf.defVar(ncid, 'y', ncPrecision, yDimID);
zVarID = netcdf.defVar(ncid, 'z', ncPrecision, zDimID);
tVarID = netcdf.defVar(ncid, 't', ncPrecision, tDimID);
netcdf.putAtt(ncid,xVarID, 'units', 'm');
netcdf.putAtt(ncid,yVarID, 'units', 'm');
netcdf.putAtt(ncid,zVarID, 'units', 'm');
netcdf.putAtt(ncid,tVarID, 'units', 's');

% Define the dynamical variables
uVarID = netcdf.defVar(ncid, 'u', ncPrecision, [xDimID,yDimID,zDimID,tDimID]);
vVarID = netcdf.defVar(ncid, 'v', ncPrecision, [xDimID,yDimID,zDimID,tDimID]);
wVarID = netcdf.defVar(ncid, 'w', ncPrecision, [xDimID,yDimID,zDimID,tDimID]);
zetaVarID = netcdf.defVar(ncid, 'zeta', ncPrecision, [xDimID,yDimID,zDimID,tDimID]);
netcdf.putAtt(ncid,uVarID, 'units', 'm/s');
netcdf.putAtt(ncid,vVarID, 'units', 'm/s');
netcdf.putAtt(ncid,wVarID, 'units', 'm/s');
netcdf.putAtt(ncid,zetaVarID, 'units', 'm');

% Define the *float* dimensions
floatDimID = netcdf.defDim(ncid, 'float_id', nFloats);
xFloatID = netcdf.defVar(ncid, 'x-position', ncPrecision, [floatDimID,tDimID]);
yFloatID = netcdf.defVar(ncid, 'y-position', ncPrecision, [floatDimID,tDimID]);
zFloatID = netcdf.defVar(ncid, 'z-position', ncPrecision, [floatDimID,tDimID]);
netcdf.putAtt(ncid,xFloatID, 'units', 'm');
netcdf.putAtt(ncid,yFloatID, 'units', 'm');
netcdf.putAtt(ncid,zFloatID, 'units', 'm');

% Define the *float* dimensions
xDiffusiveFloatID = netcdf.defVar(ncid, 'x-position-diffusive', ncPrecision, [floatDimID,tDimID]);
yDiffusiveFloatID = netcdf.defVar(ncid, 'y-position-diffusive', ncPrecision, [floatDimID,tDimID]);
zDiffusiveFloatID = netcdf.defVar(ncid, 'z-position-diffusive', ncPrecision, [floatDimID,tDimID]);
netcdf.putAtt(ncid,xDiffusiveFloatID, 'units', 'm');
netcdf.putAtt(ncid,yDiffusiveFloatID, 'units', 'm');
netcdf.putAtt(ncid,zDiffusiveFloatID, 'units', 'm');

% Define the *float* dimensions
xDrifterID = netcdf.defVar(ncid, 'x-position-drifter', ncPrecision, [floatDimID,tDimID]);
yDrifterID = netcdf.defVar(ncid, 'y-position-drifter', ncPrecision, [floatDimID,tDimID]);
zDrifterID = netcdf.defVar(ncid, 'z-position-drifter', ncPrecision, [floatDimID,tDimID]);
netcdf.putAtt(ncid,xDrifterID, 'units', 'm');
netcdf.putAtt(ncid,yDrifterID, 'units', 'm');
netcdf.putAtt(ncid,zDrifterID, 'units', 'm');

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
netcdf.putVar(ncid, setprecision(xVarID), wavemodel.x);
netcdf.putVar(ncid, setprecision(yVarID), wavemodel.y);
netcdf.putVar(ncid, setprecision(zVarID), wavemodel.z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Run the model, and write the output to NetCDF
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startTime = datetime('now');
fprintf('Starting numerical simulation on %s\n', datestr(startTime));
integrator = IntegratorWithDiffusivity( f, p0, deltaT, kappa_vector, ymin, ymax);
for iTime=1:length(t)
    if iTime == 2 || mod(iTime,10) == 0
        timePerStep = (datetime('now')-startTime)/(iTime-1);
        timeRemaining = (length(t)-iTime+1)*timePerStep;   
        fprintf('\twriting values time step %d of %d to file. Estimated finish time %s (%s from now)\n', iTime, length(t), datestr(datetime('now')+timeRemaining), datestr(timeRemaining, 'HH:MM:SS')) ;
    end
    [u,v]=wavemodel.VelocityFieldAtTime(t(iTime));
    [w,zeta] = wavemodel.VerticalFieldsAtTime(t(iTime));
    
    netcdf.putVar(ncid, setprecision(uVarID), [0 0 0 iTime-1], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], u);
    netcdf.putVar(ncid, setprecision(vVarID), [0 0 0 iTime-1], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], v);
    netcdf.putVar(ncid, setprecision(wVarID), [0 0 0 iTime-1], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], w);
    netcdf.putVar(ncid, setprecision(zetaVarID), [0 0 0 iTime-1], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], zeta);
    netcdf.putVar(ncid, setprecision(tVarID), iTime-1, 1, t(iTime));
    
    p = integrator.StepForwardToTime(t(iTime));
    netcdf.putVar(ncid, setprecision(xFloatID), [0 iTime-1], [nFloats 1], p(:,1));
    netcdf.putVar(ncid, setprecision(yFloatID), [0 iTime-1], [nFloats 1], p(:,2));
    netcdf.putVar(ncid, setprecision(zFloatID), [0 iTime-1], [nFloats 1], p(:,3));
    
    netcdf.putVar(ncid, setprecision(xDiffusiveFloatID), [0 iTime-1], [nFloats 1], p(:,4));
    netcdf.putVar(ncid, setprecision(yDiffusiveFloatID), [0 iTime-1], [nFloats 1], p(:,5));
    netcdf.putVar(ncid, setprecision(zDiffusiveFloatID), [0 iTime-1], [nFloats 1], p(:,6));
    
    netcdf.putVar(ncid, setprecision(xDrifterID), [0 iTime-1], [nFloats 1], p(:,7));
    netcdf.putVar(ncid, setprecision(yDrifterID), [0 iTime-1], [nFloats 1], p(:,8));
    netcdf.putVar(ncid, setprecision(zDrifterID), [0 iTime-1], [nFloats 1], z_float);
    

end
fprintf('Ending numerical simulation on %s\n', datestr(datetime('now')));

netcdf.close(ncid);	
