//
//  main.m
//  InternalWavesNoFloats
//
//  Created by  on 8/22/16.
//  Copyright © 2016 Jeffrey J. Early. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>
#import <GLOceanKit/GLOceanKit.h>

typedef NS_ENUM(NSUInteger, ExperimentType) {
    kSingleModeExperimentType = 0,
    kGMSpectrumExperimentType = 1
};

typedef NS_ENUM(NSUInteger, StratificationType) {
    kConstantStratificationType = 0,
    kExponentialStratificationType = 1
};

int main(int argc, const char * argv[])
{
    @autoreleasepool {
        ExperimentType experiment = kGMSpectrumExperimentType;
        StratificationType stratification = kExponentialStratificationType;
        
        GLFloat latitude = 31;
        GLFloat N0 = 5.23e-3;
        GLFloat depth, width, height;
        NSUInteger Nx, Ny, Nz;
        GLFloat maxWavePeriods = 1;
        NSString *filename;
        GLFloat amplitude;
        NSUInteger modeUnit = 1;
        NSUInteger kUnit = 1;
        NSUInteger lUnit = 1;
        NSInteger omegaSign = 1;
        NSString *strat = stratification == kConstantStratificationType ? @"Constant" : @"Exponential";
        if (experiment == kSingleModeExperimentType) {
            depth = 5000;
            width = 15e3;
            height = 7.5e3;
            Nx = 32;
            Ny = 16;
            Nz = 64;
            amplitude = 0.25; // Max horizontal wave speed, m/s
            filename = [NSString stringWithFormat: @"InternalWaveSingleMode%@Stratification.nc",strat];
        } else {
            depth = 5000;
            width = 8e3;
            height = 8e3;
            Nx = 256;
            Ny = 256;
            Nz = 64;
            maxWavePeriods = 1;
            amplitude = 1.0; // GM reference energy level
            filename = [NSString stringWithFormat: @"InternalWavesGMSpectrum%@Stratification.nc",strat];
        }
        
        GLFloat f0 = 2*(7.2921e-5)*sin(latitude*M_PI/180);
        GLFloat rho0 = 1025;
        GLFloat g = 9.81;
        
        GLFloat sampleTimeInMinutes = 10; // This will be overriden for the unit test.
        NSUInteger numStepsPerCycle = 61; // Only relevant for the single mode test, otherwise the sampleTimeInMinutes will be used.
        GLFloat omega = 0.0; // Don't set this value, it will be set for you based on the modes.
        
        // We load the transformation matrices from file, if possible.
        GLInternalWaveInitialization *wave;
        GLEquation *equation;
        
        NSString *initialConditionsFile = [[NSSearchPathForDirectoriesInDomains(NSDesktopDirectory, NSUserDomainMask, YES) objectAtIndex:0] stringByAppendingPathComponent: [NSString stringWithFormat: @"InternalWavesLatmix_%lu_%lu_%lu_%luKM_%@Stratification.internalwaves", Nx, Ny, Nz,(NSUInteger) (width/1e3), strat]];
        NSFileManager *manager = [[NSFileManager alloc] init];
        if ([manager fileExistsAtPath: initialConditionsFile isDirectory: nil])
        {
            if (!(wave = [NSKeyedUnarchiver unarchiveObjectWithFile: initialConditionsFile])) {
                NSLog(@"Failed to load wave file.");
                return 0;
            }
            equation = wave.equation;
        }
        else
        {
            equation = [[GLEquation alloc] init];
            GLDimension *zDim = [[GLDimension alloc] initDimensionWithGrid: kGLChebyshevEndpointGrid nPoints: Nz domainMin: -depth length: depth];
            zDim.name = @"z";
            GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: Nx domainMin: -width/2 length: width];
            xDim.name = @"x";
            GLDimension *yDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: Ny domainMin: -height/2 length: height];
            yDim.name = @"y";
            
            GLFunction *z = [GLFunction functionOfRealTypeFromDimension:zDim withDimensions:@[zDim] forEquation:equation];
            GLFunction *rho_bar;
            if (stratification == kConstantStratificationType) {
                rho_bar = [[z times: @(-N0*N0*rho0/g)] plus: @(rho0)];
            } else {
                // rho0*N0*N0*1300/(2*g)*(1-exp(2*z/1300))+rho0;
                rho_bar = [[[[[[z times: @(2./1300.)] exponentiate] negate] plus: @(1.)] times: @(rho0*N0*N0*1300./(2*g))] plus: @(rho0)];
            }
            rho_bar.name = @"rho_bar";
            [rho_bar solve];
            
            // The ordering the dimensions is very deliberate here, for two reasons:
            // 1. The z-dimension is in the first position, so that the horizontal fft will act on contiguous chunks of memory, and
            // 2. The last two dimensions are ordered (x,y) to appease pcolor, meshgrid, and all the standard matlab formating.
            if (experiment == kSingleModeExperimentType) {
                wave = [[GLInternalWaveInitialization alloc] initWithDensityProfile: rho_bar fullDimensions:@[zDim, xDim, yDim] latitude:latitude maxMode: modeUnit+2 equation:equation];
            } else {
                wave = [[GLInternalWaveInitialization alloc] initWithDensityProfile: rho_bar fullDimensions:@[zDim, xDim, yDim] latitude:latitude equation:equation];
            }
            
            if (![NSKeyedArchiver archiveRootObject: wave toFile: initialConditionsFile]) {
                NSLog(@"Failed to save restart file.");
            }
        }
        
        if (experiment == kSingleModeExperimentType) {
            omega = [wave createUnitWaveWithSpeed: amplitude verticalMode: modeUnit k: kUnit l: lUnit omegaSign: omegaSign];
        } else {
            [wave createGarrettMunkSpectrumWithEnergy: amplitude];
            [wave showDiagnostics];
        }
        
        // The time dimension must get set after we know what omega is.
        GLFloat maxTime = experiment==kSingleModeExperimentType ? maxWavePeriods*2*M_PI/omega : maxWavePeriods*2*M_PI/f0;
        GLFloat sampleTime = experiment==kSingleModeExperimentType ? maxTime/numStepsPerCycle : sampleTimeInMinutes*60;
        NSUInteger nTimePoints = 1 + round(maxTime/sampleTime);
        NSLog(@"Maximum wave period %d @ %02d:%02d (HH:MM)", (int) floor(maxTime/86400), ((int) floor(maxTime/3600))%24, ((int) floor(maxTime/60))%60);
        
        /************************************************************************************************/
        /*		Create the dynamical variables from the analytical solution								*/
        /************************************************************************************************/
        
        // We should check that we optimizations in places for these. They should be purely imaginary, and the multiplication and exponentiation should take advantage of that.
        GLFunction *iOmega = [[wave.eigenfrequencies swapComplex] makeHermitian];
        GLFunction *minusiOmega = [[[wave.eigenfrequencies swapComplex] negate] makeHermitian];
        
        NSArray * (^timeToUV) (GLScalar *) = ^( GLScalar *t ) {
            GLFunction *time_phase_plus = [[iOmega multiply: t] exponentiate];
            GLFunction *time_phase_minus = [[minusiOmega multiply: t] exponentiate];
            GLFunction *u = [[[wave.Sprime transform: [[wave.u_plus multiply: time_phase_plus] plus: [wave.u_minus multiply: time_phase_minus]]] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis), @(kGLDeltaBasis)]]scaleVariableBy: 1.0 withUnits: @"m/s" dimensionsBy: 1.0 units: @"m"];
            GLFunction *v = [[[wave.Sprime transform: [[wave.v_plus multiply: time_phase_plus] plus: [wave.v_minus multiply: time_phase_minus]]] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis), @(kGLDeltaBasis)]]scaleVariableBy: 1.0 withUnits: @"m/s" dimensionsBy: 1.0 units: @"m"];
            GLFunction *w = [[[wave.S transform: [[wave.w_plus multiply: time_phase_plus] plus: [wave.w_minus multiply: time_phase_minus]]] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis), @(kGLDeltaBasis)]]scaleVariableBy: 1.0 withUnits: @"m/s" dimensionsBy: 1.0 units: @"m"];
            GLFunction *rho = [[[[[wave.S transform: [[wave.rho_plus multiply: time_phase_plus] plus: [wave.rho_minus multiply: time_phase_minus]]] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis), @(kGLDeltaBasis)]] times: wave.N2] plus: wave.rho] scaleVariableBy: 1.0 withUnits: @"kg/m^3" dimensionsBy: 1.0 units: @"m"];
            GLFunction *zeta = [[[wave.S transform: [[wave.zeta_plus multiply: time_phase_plus] plus: [wave.zeta_minus multiply: time_phase_minus]]] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis), @(kGLDeltaBasis)]] scaleVariableBy: 1.0 withUnits: @"m" dimensionsBy: 1.0 units: @"m"];
            
            return @[u,v,w,rho,zeta];
        };
        
        GLScalar *t = [GLScalar scalarWithValue: 0.0*2*M_PI/f0 forEquation: equation];
        NSArray *uv = timeToUV(t);
        GLFunction *u = uv[0];
        GLFunction *v = uv[1];
        GLFunction *w = uv[2];
        GLFunction *rho = uv[3];
        GLFunction *zeta = uv[4];
        GLFunction *speed = [[[[u times: u] plus: [v times: v]] plus: [w times: w]] sqrt];
        GLFloat maxSpeed = [speed maxNow];
        NSLog(@"Initial maximum speed: %f", maxSpeed);
        

        
        /************************************************************************************************/
        /*		Create a NetCDF file and mutable variables in order to record some of the time steps.	*/
        /************************************************************************************************/
        
//        NSString *path = [[NSSearchPathForDirectoriesInDomains(NSDesktopDirectory, NSUserDomainMask, YES) objectAtIndex:0] stringByAppendingPathComponent:filename];
        NSString *path = [NSString stringWithFormat: @"/Volumes/RadiativeTr/%@", filename];
        GLNetCDFFile *netcdfFile = [[GLNetCDFFile alloc] initWithURL: [NSURL URLWithString: path] forEquation: equation overwriteExisting: YES];
        
        [netcdfFile setGlobalAttribute: @(width) forKey: @"L_domain"];
        [netcdfFile setGlobalAttribute: @(latitude) forKey: @"latitude"];
        [netcdfFile setGlobalAttribute: @(f0) forKey: @"f0"];
        [netcdfFile setGlobalAttribute: @(depth) forKey: @"D"];
        [netcdfFile setGlobalAttribute: @(amplitude) forKey: @"amplitude"];
        
        GLFunction *rhoScaled = [wave.rho scaleVariableBy: 1.0 withUnits: @"kg/m^3" dimensionsBy: 1.0 units: @"m"];
        rhoScaled.name = wave.rho.name;
        [netcdfFile addVariable: rhoScaled];
        
        GLFunction *n2Scaled = [wave.N2 scaleVariableBy: 1.0 withUnits: @"radians/s^2" dimensionsBy: 1.0 units: @"m"];
        n2Scaled.name = wave.N2.name;
        [netcdfFile addVariable: n2Scaled];
        
        GLMutableDimension *tDim = [[GLMutableDimension alloc] initWithPoints: @[@(0.0)]];
        tDim.name = @"time"; tDim.units = @"s";
        
        GLMutableVariable *uHistory = [u variableByAddingDimension: tDim];
        uHistory.name = @"u";
        uHistory = [netcdfFile addVariable: uHistory];
        
        GLMutableVariable *vHistory = [v variableByAddingDimension: tDim];
        vHistory.name = @"v";
        vHistory = [netcdfFile addVariable: vHistory];
        
        GLMutableVariable *wHistory = [w variableByAddingDimension: tDim];
        wHistory.name = @"w";
        wHistory = [netcdfFile addVariable: wHistory];
        
        GLMutableVariable *rhoHistory = [rho variableByAddingDimension: tDim];
        rhoHistory.name = @"rho";
        rhoHistory = [netcdfFile addVariable: rhoHistory];
        
        GLMutableVariable *zetaHistory = [zeta variableByAddingDimension: tDim];
        zetaHistory.name = @"zeta";
        zetaHistory = [netcdfFile addVariable: zetaHistory];
        
        for (NSUInteger iTime = 1; iTime < nTimePoints; iTime++ ) {
            @autoreleasepool {
                NSLog(@"Writing time %lu of %lu.",iTime,nTimePoints);
                
                GLFloat time = iTime * sampleTime;
                [tDim addPoint: @(time)];
                NSArray *uv = timeToUV([GLScalar scalarWithValue: time forEquation: equation]);
                
                [uHistory concatenateWithLowerDimensionalVariable: uv[0] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                [vHistory concatenateWithLowerDimensionalVariable: uv[1] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                [wHistory concatenateWithLowerDimensionalVariable: uv[2] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                [rhoHistory concatenateWithLowerDimensionalVariable: uv[3] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                [zetaHistory concatenateWithLowerDimensionalVariable: uv[4] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
            }
        }
        
        [netcdfFile close];
    }
    return 0;
}
