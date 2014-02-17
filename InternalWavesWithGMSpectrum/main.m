//
//  main.m
//  InternalWavesWithGMSpectrum
//
//  Created by Jeffrey Early on 2/11/14.
//  Copyright (c) 2014 Jeffrey J. Early. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>
#import "GLInternalModes.h"
#import "GLInternalWaveInitialization.h"

#define g 9.81

int main(int argc, const char * argv[])
{
    
	@autoreleasepool {
        // @"InternalWavesLatmix2011_16_16_128.internalwaves"
        // @"InternalWavesLatmix2011_128_128_128.internalwaves"
        NSString *restartFile = [[NSSearchPathForDirectoriesInDomains(NSDesktopDirectory, NSUserDomainMask, YES) objectAtIndex:0] stringByAppendingPathComponent:@"InternalWavesLatmix2011_128_128_64_stretched.internalwaves"];
        NSString *outputFile = @"InternalWavesLatmix.nc";
        
        //NSString *restartFile = [[NSSearchPathForDirectoriesInDomains(NSDesktopDirectory, NSUserDomainMask, YES) objectAtIndex:0] stringByAppendingPathComponent:@"InternalWavesUnitTest_128_128_128.internalwaves"];
        
        NSFileManager *manager = [[NSFileManager alloc] init];
        GLInternalWaveInitialization *wave;
        GLDimension *xDim, *yDim, *zDim;
        if ([manager fileExistsAtPath: restartFile isDirectory: nil])
        {
            wave = [NSKeyedUnarchiver unarchiveObjectWithFile: restartFile];
            zDim = wave.fullDimensions[0];
            xDim = wave.fullDimensions[1];
            yDim = wave.fullDimensions[2];
        }
        else
        {
            GLFloat latitude = 45;
            GLFloat N2 = 1e-3; //2.5e-3;
            GLFloat width = 15e3;
            GLFloat height = 15e3;
            GLFloat depth = 300;
            NSUInteger Nx = 128;
            NSUInteger Ny = 128;
            NSUInteger Nz = 64;
            
            /************************************************************************************************/
            /*		Define the problem dimensions															*/
            /************************************************************************************************/
            GLEquation *equation = [[GLEquation alloc] init];
            zDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: Nz domainMin: -depth length: depth];
            zDim.name = @"z";
            xDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: Nx domainMin: -width/2 length: width];
            xDim.name = @"x";
            yDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: Ny domainMin: -height/2 length: height];
            yDim.name = @"y";
            
            /************************************************************************************************/
            /*		Create a density profile and compute the internal wave phases                           */
            /************************************************************************************************/
            GLFunction *rho_bar;
            if (0) {
                GLNetCDFFile *profile = [[GLNetCDFFile alloc] initWithURL:[NSURL URLWithString: @"/Users/jearly/Documents/Models/InternalWaves/Latmix2011Site1Profile.nc"] forEquation:equation];
                GLFunction *rho_profile = profile.variables[0];
                GLFunction *z = [GLFunction functionOfRealTypeFromDimension:zDim withDimensions:@[zDim] forEquation:equation];
                rho_bar = [rho_profile interpolateAtPoints:@[z]];
                rho_bar.name = @"rho_bar";
            }  else if (1) {
                GLNetCDFFile *profile = [[GLNetCDFFile alloc] initWithURL:[NSURL URLWithString: @"/Users/jearly/Documents/Models/InternalWaves/Latmix2011Site1Profile_Stretched_64.nc"] forEquation:equation];
                rho_bar = profile.variables[0];
                rho_bar.name = @"rho_bar";
                zDim = rho_bar.dimensions[0];
            } else {
                GLFloat rho0 = 1025;
                GLFunction *z = [GLFunction functionOfRealTypeFromDimension:zDim withDimensions:@[zDim] forEquation:equation];
                rho_bar = [[z times: @(-N2*rho0/g)] plus: @(rho0)];
                rho_bar.name = @"rho_bar";
            }
            
            // The ordering the dimensions is very deliberate here, for two reasons:
            // 1. The z-dimension is in the first position, so that the horizontal fft will act on contiguous chunks of memory, and
            // 2. The last two dimensions are ordered (x,y) to appease pcolor, meshgrid, and all the standard matlab formating.
            wave = [[GLInternalWaveInitialization alloc] initWithDensityProfile: rho_bar fullDimensions:@[zDim, xDim, yDim] latitude:latitude equation:equation]; //@[zDim, xDim, yDim]
            
            [NSKeyedArchiver archiveRootObject: wave toFile: restartFile];
        }
        
        wave.maximumModes = 2;
        //wave.maxDepth = -300;
        [wave createGarrettMunkSpectrumWithEnergy: 0.01];
        //[wave createUnitWaveWithSpeed: 0.01 verticalMode: 1 k: 1 l: 0 omegaSign: 1];
        zDim = wave.rho.dimensions[0];
        
        GLMutableDimension *tDim = [[GLMutableDimension alloc] initWithPoints: @[@(0.0)]];
        tDim.name = @"time";

		GLFloat maxWavePeriods = 1; // The wave period is the inertial period for the GM spectrum initialization, or omega for the unit test initialization
		GLFloat sampleTimeInMinutes = 10; // This will be overriden for the unit test.
		GLFloat horizontalFloatSpacingInMeters = 1000;
		GLFloat verticalFloatSpacingInMeters = 17;
        
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
            GLFunction *rho = [[[[wave.S transform: [[wave.rho_plus multiply: time_phase_plus] plus: [wave.rho_minus multiply: time_phase_minus]]] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis), @(kGLDeltaBasis)]] times: wave.N2] plus: wave.rho];
            GLFunction *zeta = [[wave.S transform: [[wave.zeta_plus multiply: time_phase_plus] plus: [wave.zeta_minus multiply: time_phase_minus]]] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis), @(kGLDeltaBasis)]];
            
            
			
			return @[u,v,w,rho,zeta];
		};
        
		GLScalar *t = [GLScalar scalarWithValue: 0.0*2*M_PI/wave.f0 forEquation: wave.equation];
		NSArray *uv = timeToUV(t);
		GLFunction *u = uv[0];
		GLFunction *v = uv[1];
		GLFunction *w = uv[2];
        GLFunction *rho = uv[3];
        GLFunction *zeta = uv[4];
        GLFunction *speed = [[[[u times: u] plus: [v times: v]] plus: [w times: w]] sqrt];
		GLFloat maxSpeed = [speed maxNow];
		NSLog(@"Initial maximum speed: %g", maxSpeed);
		
		/************************************************************************************************/
		/*		Let's also plop a float at a bunch of grid points.                                      */
		/************************************************************************************************/
        
        GLDimension *xFloatDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: ceil(xDim.domainLength/horizontalFloatSpacingInMeters) domainMin: -xDim.domainLength/2 length:xDim.domainLength];
		xFloatDim.name = @"x-float";
		GLDimension *yFloatDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:ceil(yDim.domainLength/horizontalFloatSpacingInMeters) domainMin: -yDim.domainLength/2  length:yDim.domainLength];
		yFloatDim.name = @"y-float";
		GLDimension *zFloatDim = [[GLDimension alloc] initWithPoints: @[ @(-38), @(-31.5), @(-25)]];
		zFloatDim.name = @"z-float";
        
		// For consistency, we order the float dimensions the same as the dynamical variable dimensions.
        NSArray *floatDimensions = @[zFloatDim, xFloatDim, yFloatDim];
		GLFunction *xFloat = [GLFunction functionOfRealTypeFromDimension: xFloatDim withDimensions: floatDimensions forEquation: wave.equation];
		GLFunction *yFloat = [GLFunction functionOfRealTypeFromDimension: yFloatDim withDimensions: floatDimensions forEquation: wave.equation];
		GLFunction *zFloat = [GLFunction functionOfRealTypeFromDimension: zFloatDim withDimensions: floatDimensions forEquation: wave.equation];
        
		GLFunction *xPosition = [GLFunction functionFromFunction: xFloat];
		GLFunction *yPosition = [GLFunction functionFromFunction: yFloat];
		GLFunction *zPosition = [GLFunction functionFromFunction: zFloat];
        GLFunction *isopycnalDeviation = [zeta interpolateAtPoints:@[zPosition, xPosition, yPosition]];
        zPosition = [zPosition plus: isopycnalDeviation];
		
		CGFloat cfl = 0.25;
        GLFloat cflTimeStep = cfl * xDim.sampleInterval / maxSpeed;
		GLFloat outputTimeStep = sampleTimeInMinutes*60;
		
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
		
		NSString *path = [[NSSearchPathForDirectoriesInDomains(NSDesktopDirectory, NSUserDomainMask, YES) objectAtIndex:0] stringByAppendingPathComponent:outputFile];
		GLNetCDFFile *netcdfFile = [[GLNetCDFFile alloc] initWithURL: [NSURL URLWithString: path] forEquation: wave.equation overwriteExisting: YES];
		
        [netcdfFile setGlobalAttribute: @(wave.f0) forKey: @"f0"];
		
		[netcdfFile addVariable: wave.rho];
        
//        [(GLDimension *)wave.internalModes.rho.dimensions[0] setName: @"z_full"];
//        [(GLFunction *)wave.internalModes.rho setName: @"rho_bar_full"];
//        [(GLFunction *)wave.internalModes.N2 setName: @"N2_full"];
//        [netcdfFile addVariable: wave.internalModes.rho];
//        [netcdfFile addVariable: wave.internalModes.N2];
        
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
		
		GLMutableVariable *xPositionHistory = [xPosition variableByAddingDimension: tDim];
		xPositionHistory.name = @"x-position";
		xPositionHistory = [netcdfFile addVariable: xPositionHistory];
        
		GLMutableVariable *yPositionHistory = [yPosition variableByAddingDimension: tDim];
		yPositionHistory.name = @"y-position";
		yPositionHistory = [netcdfFile addVariable: yPositionHistory];
		
		GLMutableVariable *zPositionHistory = [zPosition variableByAddingDimension: tDim];
		zPositionHistory.name = @"z-position";
		zPositionHistory = [netcdfFile addVariable: zPositionHistory];
		
		/************************************************************************************************/
		/*		Time step the analytical solution and the integrator forward, then write out the data.	*/
		/************************************************************************************************/
		
		GLFloat maxTime = maxWavePeriods*2*M_PI/wave.f0;
        NSLog(@"Maximum wave period %d @ %02d:%02d (HH:MM)", (int) floor(maxTime/86400), ((int) floor(maxTime/3600))%24, ((int) floor(maxTime/60))%60);
        GLFloat sampleTime = sampleTimeInMinutes*60;
        for (GLFloat time = sampleTime; time <= maxTime-0*sampleTime; time += sampleTime)
            //while( integrator.currentTime < maxTime + integrator.stepSize/2)
        {
            @autoreleasepool {
                NSLog(@"Logging day %d @ %02d:%02d (HH:MM)", (int) floor(time/86400), ((int) floor(time/3600))%24, ((int) floor(time/60))%60);
                //NSArray *yout = [integrator stepForward];
//				NSArray *yout = [integrator stepForwardToTime: time];
//                
//				//GLFloat time = integrator.currentTime;
//                NSLog(@"Logging day %d @ %02d:%02d (HH:MM), last step size: %02d:%02.1f (MM:SS.S).", (int) floor(time/86400), ((int) floor(time/3600))%24, ((int) floor(time/60))%60, (int)floor(integrator.lastStepSize/60), fmod(integrator.lastStepSize,60));
//				NSLog(@"Logging day %d @ %02d:%02d (HH:MM), last step size: %02d:%02.1f (MM:SS.S).", (int) floor(integrator.currentTime/86400), ((int) floor(integrator.currentTime/3600))%24, ((int) floor(integrator.currentTime/60))%60, (int)floor(integrator.lastStepSize/60), fmod(integrator.lastStepSize,60));
                
                NSArray *uv = timeToUV([GLScalar scalarWithValue: time forEquation: wave.equation]);
				
				[tDim addPoint: @(time)];
				[uHistory concatenateWithLowerDimensionalVariable: uv[0] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                [vHistory concatenateWithLowerDimensionalVariable: uv[1] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
				[wHistory concatenateWithLowerDimensionalVariable: uv[2] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
				[rhoHistory concatenateWithLowerDimensionalVariable: uv[3] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                [zetaHistory concatenateWithLowerDimensionalVariable: uv[4] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
				
//                [xPositionHistory concatenateWithLowerDimensionalVariable: yout[0] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
//                [yPositionHistory concatenateWithLowerDimensionalVariable: yout[1] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
//				[zPositionHistory concatenateWithLowerDimensionalVariable: yout[2] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                
                [netcdfFile waitUntilAllOperationsAreFinished];
            }
        }
        
		[netcdfFile close];
	}
    return 0;
}