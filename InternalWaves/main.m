//
//  main.m
//  InternalWaves
//
//  Created by Jeffrey J. Early on 12/9/13.
//  Copyright (c) 2013 Jeffrey J. Early. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>
#import <GLOceanKit/GLOceanKit.h>

typedef NS_ENUM(NSUInteger, ExperimentType) {
    kSingleModeExperimentType = 0,
    kGMSpectrumExperimentType = 1
};

int main(int argc, const char * argv[])
{
	@autoreleasepool {
        ExperimentType experiment = kSingleModeExperimentType;
        
        GLFloat latitude = 45;
        GLFloat N2 = 2.5e-3;
        GLFloat depth, width, height;
        NSUInteger Nx, Ny, Nz;
        GLFloat maxWavePeriods = 1;
        NSString *filename;
        if (experiment == kSingleModeExperimentType) {
            depth = 100;
            width = 15e3;
            height = 7.5e3;
            Nx = 32;
            Ny = 16;
            Nz = 24;
            filename = @"InternalWaveSingleMode.nc";
        } else {
            depth = 100;
            width = 15e3;
            height = 15e3;
            Nx = 128;
            Ny = 128;
            Nz = 64;
            maxWavePeriods = 1;
            filename = @"InternalWavesGMSpectrum.nc";
        }

		GLFloat horizontalFloatSpacingInMeters = 1000;
		GLFloat verticalFloatSpacingInMeters = 25;
		
		GLFloat f0 = 2*(7.2921e-5)*sin(latitude*M_PI/180);
		GLFloat rho0 = 1025;
		GLFloat g = 9.81;
        
        GLFloat sampleTimeInMinutes = 10; // This will be overriden for the unit test.
        NSUInteger numStepsPerCycle = 61; // Only relevant for the single mode test, otherwise the sampleTimeInMinutes will be used.
		GLFloat omega = 0.0; // Don't set this value, it will be set for you based on the modes.
        
        /************************************************************************************************/
		/*		Define the problem dimensions															*/
		/************************************************************************************************/


		GLEquation *equation = [[GLEquation alloc] init];
		GLDimension *zDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: Nz domainMin: -depth length: depth];
		zDim.name = @"z";
		GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: Nx domainMin: -width/2 length: width];
		xDim.name = @"x";
		GLDimension *yDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: Ny domainMin: -height/2 length: height];
		yDim.name = @"y";

        /************************************************************************************************/
		/*		Create a density profile and compute the internal wave phases                           */
		/************************************************************************************************/
        GLFunction *z = [GLFunction functionOfRealTypeFromDimension:zDim withDimensions:@[zDim] forEquation:equation];
        GLFunction *rho_bar = [[z times: @(-N2*rho0/g)] plus: @(rho0)];
		rho_bar.name = @"rho_bar";
		
		// The ordering the dimensions is very deliberate here, for two reasons:
		// 1. The z-dimension is in the first position, so that the horizontal fft will act on contiguous chunks of memory, and
		// 2. The last two dimensions are ordered (x,y) to appease pcolor, meshgrid, and all the standard matlab formating.
        GLInternalWaveInitialization *wave;
        
        if (experiment == kSingleModeExperimentType) {
            NSUInteger modeUnit = 1;
            NSUInteger kUnit = 1;
            NSUInteger lUnit = 1;
            NSInteger omegaSign = 1;
            GLFloat U_max = .25;
			wave = [[GLInternalWaveInitialization alloc] initWithDensityProfile: rho_bar fullDimensions:@[zDim, xDim, yDim] latitude:latitude maxMode: modeUnit+2 equation:equation];
            omega = [wave createUnitWaveWithSpeed: U_max verticalMode: modeUnit k: kUnit l: lUnit omegaSign: omegaSign];
        } else {
			wave = [[GLInternalWaveInitialization alloc] initWithDensityProfile: rho_bar fullDimensions:@[zDim, xDim, yDim] latitude:latitude equation:equation];
            [wave createGarrettMunkSpectrumWithEnergy: 0.1];
			
//			NSString *path = [[NSSearchPathForDirectoriesInDomains(NSDesktopDirectory, NSUserDomainMask, YES) objectAtIndex:0] stringByAppendingPathComponent:@"InternalWavesGMSpectrum.internalwaves"];
//            [NSKeyedArchiver archiveRootObject: wave toFile: path];
//            wave = [NSKeyedUnarchiver unarchiveObjectWithFile: path];
		}
        
        // The time dimension must get set after we know what omega is.
        GLFloat maxTime = experiment==kSingleModeExperimentType ? maxWavePeriods*2*M_PI/omega : maxWavePeriods*2*M_PI/f0;
        GLFloat sampleTime = experiment==kSingleModeExperimentType ? maxTime/numStepsPerCycle : sampleTimeInMinutes*60;
        NSLog(@"Maximum wave period %d @ %02d:%02d (HH:MM)", (int) floor(maxTime/86400), ((int) floor(maxTime/3600))%24, ((int) floor(maxTime/60))%60);
        GLDimension *tDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: 1 + round(maxTime/sampleTime)  domainMin:0 length:maxTime];
        tDim.name = @"time"; tDim.units = @"s";
        
		/************************************************************************************************/
		/*		Create the dynamical variables from the analytical solution								*/
		/************************************************************************************************/
		
		// We should check that we optimizations in places for these. They should be purely imaginary, and the multiplication and exponentiation should take advantage of that.
		GLFunction *iOmega = [[wave.eigenfrequencies swapComplex] makeHermitian];
		GLFunction *minusiOmega = [[[wave.eigenfrequencies swapComplex] negate] makeHermitian];

		NSArray * (^timeToUV) (GLScalar *) = ^( GLScalar *t ) {
			GLFunction *time_phase_plus = [[iOmega multiply: t] exponentiate];
			GLFunction *time_phase_minus = [[minusiOmega multiply: t] exponentiate];
			GLFunction *u = [[wave.Sprime transform: [[wave.u_plus multiply: time_phase_plus] plus: [wave.u_minus multiply: time_phase_minus]]] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis), @(kGLDeltaBasis)]];
			GLFunction *v = [[wave.Sprime transform: [[wave.v_plus multiply: time_phase_plus] plus: [wave.v_minus multiply: time_phase_minus]]] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis), @(kGLDeltaBasis)]];
			GLFunction *w = [[wave.S transform: [[wave.w_plus multiply: time_phase_plus] plus: [wave.w_minus multiply: time_phase_minus]]] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis), @(kGLDeltaBasis)]];
            GLFunction *rho = [[[[wave.S transform: [[wave.rho_plus multiply: time_phase_plus] plus: [wave.rho_minus multiply: time_phase_minus]]] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis), @(kGLDeltaBasis)]] times: wave.N2] plus: rho_bar];
            GLFunction *zeta = [[wave.S transform: [[wave.zeta_plus multiply: time_phase_plus] plus: [wave.zeta_minus multiply: time_phase_minus]]] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis), @(kGLDeltaBasis)]];
			
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
		/*		Let's also plop a float at a bunch of grid points.                                      */
		/************************************************************************************************/
        
        GLDimension *xFloatDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: ceil(width/horizontalFloatSpacingInMeters) domainMin: -width/2 length:width];
		xFloatDim.name = @"x-float";
		GLDimension *yFloatDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:ceil(height/horizontalFloatSpacingInMeters) domainMin: -height/2  length:height];
		yFloatDim.name = @"y-float";
		GLDimension *zFloatDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints:ceil(depth/verticalFloatSpacingInMeters) domainMin: -depth  length:depth];
		zFloatDim.name = @"z-float";
        
		// For consistency, we order the float dimensions the same as the dynamical variable dimensions.
        NSArray *floatDimensions = @[zFloatDim, xFloatDim, yFloatDim];
		GLFunction *xFloat = [GLFunction functionOfRealTypeFromDimension: xFloatDim withDimensions: floatDimensions forEquation: equation];
		GLFunction *yFloat = [GLFunction functionOfRealTypeFromDimension: yFloatDim withDimensions: floatDimensions forEquation: equation];
		GLFunction *zFloat = [GLFunction functionOfRealTypeFromDimension: zFloatDim withDimensions: floatDimensions forEquation: equation];
        
		GLFunction *xPosition = [GLFunction functionFromFunction: xFloat];
		GLFunction *yPosition = [GLFunction functionFromFunction: yFloat];
		GLFunction *zPosition = [GLFunction functionFromFunction: zFloat];
        GLFunction *isopycnalDeviation = [zeta interpolateAtPoints:@[zPosition, xPosition, yPosition]];
        zPosition = [zPosition plus: isopycnalDeviation];
		
		CGFloat cfl = 0.25;
        GLFloat cflTimeStep = cfl * xDim.sampleInterval / maxSpeed;
		GLFloat outputTimeStep = experiment==kSingleModeExperimentType ? 2*M_PI/(numStepsPerCycle*omega) : sampleTimeInMinutes*60;
		
		GLFloat timeStep = cflTimeStep > outputTimeStep ? outputTimeStep : outputTimeStep / ceil(outputTimeStep/cflTimeStep);
		
        GLAdaptiveRungeKuttaOperation *integrator = [GLAdaptiveRungeKuttaOperation rungeKutta4AdvanceY: @[xPosition, yPosition, zPosition] stepSize: timeStep fFromTY:^(GLScalar *time, NSArray *yNew) {
			NSArray *uv = timeToUV(time);
			GLFunction *u2 = uv[0];
			GLFunction *v2 = uv[1];
			GLFunction *w2 = uv[2];
			GLSimpleInterpolationOperation *interp = [[GLSimpleInterpolationOperation alloc] initWithFirstOperand: @[u2, v2, w2] secondOperand: @[yNew[2], yNew[0], yNew[1]]];
			return interp.result;
		}];
		//integrator.absoluteTolerance = @[ @(1e0), @(1e0)];
		//integrator.relativeTolerance = @[ @(1e-6), @(1e-6), @(1e-6)];
        
        /************************************************************************************************/
		/*		Create a NetCDF file and mutable variables in order to record some of the time steps.	*/
		/************************************************************************************************/
		
		NSString *path = [[NSSearchPathForDirectoriesInDomains(NSDesktopDirectory, NSUserDomainMask, YES) objectAtIndex:0] stringByAppendingPathComponent:filename];
		GLNetCDFFile *netcdfFile = [[GLNetCDFFile alloc] initWithURL: [NSURL URLWithString: path] forEquation: equation overwriteExisting: YES];
		
		[netcdfFile setGlobalAttribute: @(width) forKey: @"L_domain"];
        [netcdfFile setGlobalAttribute: @(latitude) forKey: @"latitude"];
        [netcdfFile setGlobalAttribute: @(f0) forKey: @"f0"];
        [netcdfFile setGlobalAttribute: @(depth) forKey: @"D"];
		
        GLFunction *rhoScaled = [rho_bar scaleVariableBy: 1.0 withUnits: @"kg/m^3" dimensionsBy: 1.0 units: @"m"];
        rhoScaled.name = rho_bar.name;
		[netcdfFile addVariable: rhoScaled];
		
        integrator.shouldDisplayProgress = YES;
        [integrator integrateAlongDimension: tDim withTimeScale: 1.0 file: netcdfFile output: ^(GLScalar *t, NSArray *y) {
            NSMutableDictionary *scaledVariables = [NSMutableDictionary dictionary];
            NSArray *uv = timeToUV(t);
            
            scaledVariables[@"u"] = [uv[0] scaleVariableBy: 1.0 withUnits: @"m/s" dimensionsBy: 1.0 units: @"m"];
            scaledVariables[@"v"] = [uv[1] scaleVariableBy: 1.0 withUnits: @"m/s" dimensionsBy: 1.0 units: @"m"];
            scaledVariables[@"w"] = [uv[2] scaleVariableBy: 1.0 withUnits: @"m/s" dimensionsBy: 1.0 units: @"m"];
            scaledVariables[@"rho"] = [uv[3] scaleVariableBy: 1.0 withUnits: @"kg/m^3" dimensionsBy: 1.0 units: @"m"];
            scaledVariables[@"zeta"] = [uv[4] scaleVariableBy: 1.0 withUnits: @"m" dimensionsBy: 1.0 units: @"m"];
            scaledVariables[@"x-position"] = [y[0] scaleVariableBy: 1.0 withUnits: @"m" dimensionsBy: 1.0 units: @"float_id"];
            scaledVariables[@"y-position"] = [y[1] scaleVariableBy: 1.0 withUnits: @"m" dimensionsBy: 1.0 units: @"float_id"];
            scaledVariables[@"z-position"] = [y[2] scaleVariableBy: 1.0 withUnits: @"m" dimensionsBy: 1.0 units: @"float_id"];
            
            return scaledVariables;
        }];
        
		[netcdfFile close];
	}
    return 0;
}

