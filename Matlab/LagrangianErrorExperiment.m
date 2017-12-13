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

outputInterval = 15*60;
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

shouldUseGMSpectrum = 1;

wavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

if shouldUseGMSpectrum == 1
    wavemodel.FillOutWaveSpectrum();
    wavemodel.InitializeWithGMSpectrum(GMReferenceLevel);
    
    %     Ap = wavemodel.Amp_plus;
    %     Am = wavemodel.Amp_minus;
    %     Ap(:,:,(N+1):end) = 0;
    %     Am(:,:,(N+1):end) = 0;
    %     wavemodel.GenerateWavePhases(Ap,Am);
    
%     Kh = sqrt(wavemodel.K.^2 + wavemodel.L.^2);
%     k = wavemodel.k;
%     k_nyquist = max(k) - 2*(k(2)-k(1));
%     dk = k(2)-k(1);
%     %     indices = Kh < (k_nyquist - dk/2) | Kh > (k_nyquist + dk/2);
%     k_cutoff = 12*dk;
%     indices = Kh < k_cutoff; % & Kh ~= 0;
%     
%     A_plus = wavemodel.Amp_plus;
%     A_minus = wavemodel.Amp_minus;
%     A_plus(indices) = 0;
%     A_minus(indices) = 0;
%     
%     wavemodel.GenerateWavePhases(A_plus,A_minus);
    
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
    %     maxTime = period;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create floats/drifters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx = wavemodel.x(2)-wavemodel.x(1);
dy = wavemodel.y(2)-wavemodel.y(1);
N = 5;
nLevels = 1;
x_float = (1:N)*dx;
y_float = (1:N)*dy;
% z_float = (0:nLevels-1)*(-Lz/(2*(nLevels-1)));
z_float = -Lz/2;

[x_float,y_float,z_float] = ndgrid(x_float,y_float,z_float);
x_float = reshape(x_float,[],1);
y_float = reshape(y_float,[],1);
z_float = reshape(z_float,[],1);
nFloats = numel(x_float);

% Now let's place the floats along an isopycnal.
isopycnalDeviation = wavemodel.IsopycnalDisplacementAtTimePosition(0,x_float,y_float,z_float, 'exact');
z_isopycnal = z_float + isopycnalDeviation;

% Iteratively place floats on the isopycnal surface. Overkill, probably.
for zLevel = 1:nLevels
    zLevelIndices = (zLevel-1)*N*N + (1:(N*N));
    for i=1:15
        rho = wavemodel.DensityAtTimePosition(0,x_float(zLevelIndices),y_float(zLevelIndices),z_isopycnal(zLevelIndices), 'exact');
        dRho = rho - mean(rho);
        dz = dRho * 9.81/(N0*N0*wavemodel.rho0);
        z_isopycnal(zLevelIndices) = z_isopycnal(zLevelIndices)+dz;
    end
    
    rho = wavemodel.DensityAtTimePosition(0,x_float(zLevelIndices),y_float(zLevelIndices),z_isopycnal(zLevelIndices), 'exact');
    dRho = rho - mean(rho);
    dz = dRho * 9.81/(N0*N0*wavemodel.rho0);
    fprintf('All floats are within %.2g meters of the isopycnal at z=%.1f meters\n',max(abs(dz)),z_float((zLevel-1)*N*N+1))
end

ymin = [-Inf -Inf -Lz -Inf -Inf -Lz -Inf -Inf -Lz];
ymax = [Inf Inf 0 Inf Inf 0 Inf Inf 0];
kappa_vector = zeros(size(ymin));
p0 = cat(2, x_float, y_float, z_isopycnal, x_float, y_float, z_isopycnal, x_float, y_float, z_isopycnal);

f = @(t,y) FluxForLagrangianErrorExperiment(t,y,wavemodel);

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

filepath = sprintf('%s/LagrangianErrorExperiment_%s_%dx%dx%d.nc', outputfolder,datestr(datetime('now'),'yyyy-mm-ddTHHMMSS'),Nx,Ny,Nz);

% Apple uses 1e9 bytes as 1 GB (not the usual multiples of 2 definition)
totalFields = 3;
totalSize = totalFields*bytePerFloat*length(t)*3*nFloats/1e6;
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

integrator = IntegratorWithDiffusivity( f, p0, deltaT, kappa_vector, ymin, ymax);
startTime = datetime('now');
fprintf('Starting numerical simulation on %s\n', datestr(startTime));
for iTime=1:length(t)
    if iTime == 2
       startTime = datetime('now'); 
    end
    if iTime == 3 || mod(iTime,10) == 0
        timePerStep = (datetime('now')-startTime)/(iTime-2);
        timeRemaining = (length(t)-iTime+1)*timePerStep;   
        fprintf('\twriting values time step %d of %d to file. Estimated finish time %s (%s from now)\n', iTime, length(t), datestr(datetime('now')+timeRemaining), datestr(timeRemaining, 'HH:MM:SS')) ;
    end

    netcdf.putVar(ncid, setprecision(tVarID), iTime-1, 1, t(iTime));
    
    p = integrator.StepForwardToTime(t(iTime));
    netcdf.putVar(ncid, setprecision(xFloatID), [0 iTime-1], [nFloats 1], p(:,1));
    netcdf.putVar(ncid, setprecision(yFloatID), [0 iTime-1], [nFloats 1], p(:,2));
    netcdf.putVar(ncid, setprecision(zFloatID), [0 iTime-1], [nFloats 1], p(:,3));
    netcdf.putVar(ncid, setprecision(densityFloatID), [0 iTime-1], [nFloats 1], wavemodel.DensityAtTimePosition(t(iTime),p(:,1),p(:,2),p(:,3), 'exact')-wavemodel.rho0);
    
    netcdf.putVar(ncid, setprecision(xLinearFloatID), [0 iTime-1], [nFloats 1], p(:,4));
    netcdf.putVar(ncid, setprecision(yLinearFloatID), [0 iTime-1], [nFloats 1], p(:,5));
    netcdf.putVar(ncid, setprecision(zLinearFloatID), [0 iTime-1], [nFloats 1], p(:,6));
    netcdf.putVar(ncid, setprecision(densityLinearFloatID), [0 iTime-1], [nFloats 1], wavemodel.DensityAtTimePosition(t(iTime),p(:,4),p(:,5),p(:,6), 'exact')-wavemodel.rho0);
    
    netcdf.putVar(ncid, setprecision(xSplineFloatID), [0 iTime-1], [nFloats 1], p(:,7));
    netcdf.putVar(ncid, setprecision(ySplineFloatID), [0 iTime-1], [nFloats 1], p(:,8));
    netcdf.putVar(ncid, setprecision(zSplineFloatID), [0 iTime-1], [nFloats 1], p(:,9));
    netcdf.putVar(ncid, setprecision(densitySplineFloatID), [0 iTime-1], [nFloats 1], wavemodel.DensityAtTimePosition(t(iTime),p(:,7),p(:,8),p(:,9), 'exact')-wavemodel.rho0);
end
fprintf('Ending numerical simulation on %s\n', datestr(datetime('now')));

netcdf.close(ncid);	
