%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LagrangianErrorExperiment
%
% This script uses the InternalWaveModel to generate the time-evolution of
% a Garrett-Munk spectrum of internal waves and save the output to a NetCDF
% file. It advects particles using the exact (spectral) interpolation,
% linear interpolation and cubic-spline interpolation. The goal is to
% assess the errors of these techniques.
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com
%
% May 2nd, 2017      Version 1.0

N = 64;
aspectRatio = 1;

L = 50e3;
Lx = aspectRatio*L;
Ly = L;
Lz = 5000;

Nx = aspectRatio*N;
Ny = N;
Nz = N+1; % Must include end point to advect at the surface, so use 2^N + 1

latitude = 31;
N0 = 5.2e-3; % Choose your stratification
GMReferenceLevel = 1.0 * Lz/1300;

outputInterval = 60;
maxTime = 86400/8;

outputfolder = '/Volumes/OceanTransfer';
% outputfolder = '/Users/jearly/Desktop';

precision = 'double';

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

shouldUseGMSpectrum = 0;

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
    U = 0.025; % m/s
    sign = 1;
    phi = 0;
    k0 = 2;
    l0 = 0;
    alpha = atan2(l0,k0);
    k = 2*pi*sqrt(k0^2 + l0^2)/Lx;
    
    period = wavemodel.InitializeWithPlaneWave(k0,l0,j0,U,sign);
    wavemodel.AddGriddedWavesWithWavemodes(floor(N/3),l0,floor(3*N/4),0,U/4,sign);
    
%     [omega, alpha, k, l, mode, phi, A] = wavemodel.WaveCoefficientsFromGriddedWaves();
%     wavemodel.RemoveAllGriddedWaves();
%     wavemodel.SetExternalWavesWithFrequencies(omega, alpha, mode, phi, A,Normalization.kConstant);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create floats/drifters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx = wavemodel.x(2)-wavemodel.x(1);
dy = wavemodel.y(2)-wavemodel.y(1);
N = 1;
nLevels = 1;
x_float = (1:N)*dx;
y_float = (1:N)*dy;
z_float = -Lz/2.2165;

[x_float,y_float,z_float] = ndgrid(x_float,y_float,z_float);
x_float = reshape(x_float,[],1);
y_float = reshape(y_float,[],1);
z_float = reshape(z_float,[],1);
nFloats = numel(x_float);

% Now let's place the floats along an isopycnal.
isopycnalDeviation = wavemodel.IsopycnalDisplacementAtTimePosition(0,x_float,y_float,z_float, 'exact');
z_isopycnal = z_float + isopycnalDeviation;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Determine the proper time interval
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


t_in = (0:outputInterval:maxTime)';
p0 = [x_float, y_float, z_float];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create a NetCDF file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filepath = sprintf('%s/LagrangianErrorExperiment_%s_%dx%dx%d.nc', outputfolder,datestr(datetime('now'),'yyyy-mm-ddTHHMMSS'),Nx,Ny,Nz);

% Apple uses 1e9 bytes as 1 GB (not the usual multiples of 2 definition)
totalFields = 3;
totalSize = totalFields*bytePerFloat*length(t_in)*3*nFloats/1e6;
fprintf('Writing output file to %s\nExpected file size is %.2f MB.\n',filepath,totalSize);

cmode = netcdf.getConstant('CLOBBER');
cmode = bitor(cmode,netcdf.getConstant('SHARE'));
ncid = netcdf.create(filepath, cmode);

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

% Define the *float* dimensions
floatDimID = netcdf.defDim(ncid, 'float_id', nFloats);
xFloatID = netcdf.defVar(ncid, 'x-position-exact', ncPrecision, [floatDimID,tDimID]);
yFloatID = netcdf.defVar(ncid, 'y-position-exact', ncPrecision, [floatDimID,tDimID]);
zFloatID = netcdf.defVar(ncid, 'z-position-exact', ncPrecision, [floatDimID,tDimID]);
densityFloatID = netcdf.defVar(ncid, 'density-exact', ncPrecision, [floatDimID,tDimID]);
netcdf.putAtt(ncid,xFloatID, 'units', 'm');
netcdf.putAtt(ncid,yFloatID, 'units', 'm');
netcdf.putAtt(ncid,zFloatID, 'units', 'm');

% Define the *float* dimensions
xLinearFloatID = netcdf.defVar(ncid, 'x-position-linear', ncPrecision, [floatDimID,tDimID]);
yLinearFloatID = netcdf.defVar(ncid, 'y-position-linear', ncPrecision, [floatDimID,tDimID]);
zLinearFloatID = netcdf.defVar(ncid, 'z-position-linear', ncPrecision, [floatDimID,tDimID]);
densityLinearFloatID = netcdf.defVar(ncid, 'density-linear', ncPrecision, [floatDimID,tDimID]);
netcdf.putAtt(ncid,xLinearFloatID, 'units', 'm');
netcdf.putAtt(ncid,yLinearFloatID, 'units', 'm');
netcdf.putAtt(ncid,zLinearFloatID, 'units', 'm');

% Define the *float* dimensions
xSplineFloatID = netcdf.defVar(ncid, 'x-position-spline', ncPrecision, [floatDimID,tDimID]);
ySplineFloatID = netcdf.defVar(ncid, 'y-position-spline', ncPrecision, [floatDimID,tDimID]);
zSplineFloatID = netcdf.defVar(ncid, 'z-position-spline', ncPrecision, [floatDimID,tDimID]);
densitySplineFloatID = netcdf.defVar(ncid, 'density-spline', ncPrecision, [floatDimID,tDimID]);
netcdf.putAtt(ncid,xSplineFloatID, 'units', 'm');
netcdf.putAtt(ncid,ySplineFloatID, 'units', 'm');
netcdf.putAtt(ncid,zSplineFloatID, 'units', 'm');

% Write some metadata
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'latitude', latitude);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'N0', N0);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'GMReferenceLevel', GMReferenceLevel);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'Model', 'Created from InternalWaveModel.m written by Jeffrey J. Early.');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'ModelVersion', wavemodel.version);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'CreationDate', datestr(datetime('now')));

netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'nFloatLevels', nLevels);

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

netcdf.putVar(ncid, tVarID, 0, length(t_in), t_in);

fprintf('Spectral interpolation...\n');
f = @(t,y) wavemodel.VelocityAtTimePositionVector(t,y,'exact');
[t,p] = ode45(f,t_in, p0,odeset('RelTol',1e-11,'AbsTol',1e-8));
x = p(:,1);
y = p(:,2);
z = p(:,3);
rho = zeros(size(x));
for iTime=1:length(t)
    rho(iTime) = wavemodel.DensityAtTimePosition(t(iTime),x(iTime),y(iTime),z(iTime), 'exact') - wavemodel.rho0;
end
netcdf.putVar(ncid, xFloatID, setprecision(x));
netcdf.putVar(ncid, yFloatID, setprecision(y));
netcdf.putVar(ncid, zFloatID, setprecision(z));
netcdf.putVar(ncid, densityFloatID, setprecision(rho));


fprintf('Linear interpolation...\n');
f = @(t,y) wavemodel.VelocityAtTimePositionVector(t,y,'linear');
[t,p] = ode45(f,t_in, p0,odeset('RelTol',1e-11,'AbsTol',1e-8));
x = p(:,1);
y = p(:,2);
z = p(:,3);
rho = zeros(size(x));
for iTime=1:length(t)
    rho(iTime) = wavemodel.DensityAtTimePosition(t(iTime),x(iTime),y(iTime),z(iTime), 'exact') - wavemodel.rho0;
end
netcdf.putVar(ncid, xLinearFloatID, setprecision(x));
netcdf.putVar(ncid, yLinearFloatID, setprecision(y));
netcdf.putVar(ncid, zLinearFloatID, setprecision(z));
netcdf.putVar(ncid, densityLinearFloatID, setprecision(rho));


fprintf('Spline interpolation...\n');
f = @(t,y) wavemodel.VelocityAtTimePositionVector(t,y,'spline');
[t,p] = ode45(f,t_in, p0,odeset('RelTol',1e-11,'AbsTol',1e-8));
x = p(:,1);
y = p(:,2);
z = p(:,3);
rho = zeros(size(x));
for iTime=1:length(t)
    rho(iTime) = wavemodel.DensityAtTimePosition(t(iTime),x(iTime),y(iTime),z(iTime), 'exact') - wavemodel.rho0;
end
netcdf.putVar(ncid, xSplineFloatID, setprecision(x));
netcdf.putVar(ncid, ySplineFloatID, setprecision(y));
netcdf.putVar(ncid, zSplineFloatID, setprecision(z));
netcdf.putVar(ncid, densitySplineFloatID, setprecision(rho));

fprintf('Ending numerical simulation on %s\n', datestr(datetime('now')));

netcdf.close(ncid);	
