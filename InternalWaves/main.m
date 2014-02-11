//
//  main.m
//  InternalWaves
//
//  Created by Jeffrey J. Early on 12/9/13.
//  Copyright (c) 2013 Jeffrey J. Early. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>
#import "GLInternalModes.h"
#import "GLInternalWaveInitialization.h"

int main(int argc, const char * argv[])
{

	@autoreleasepool {
		GLFloat latitude = 45;
		GLFloat N2 = 1e-4;
		GLFloat depth = 100;
		GLFloat width = 15e3;
        GLFloat height = 7.5e3;
		NSUInteger Nx = 32;
        NSUInteger Ny = 16;
		NSUInteger Nz = 100;
		GLFloat maxWavePeriods = 1; // The wave period is the inertial period for the GM spectrum initialization, or omega for the unit test initialization
		GLFloat sampleTimeInMinutes = 10; // This will be overriden for the unit test.
		GLFloat horizontalFloatSpacingInMeters = 100;
		GLFloat verticalFloatSpacingInMeters = 25;
		
		GLFloat f0 = 2*(7.2921e-5)*sin(latitude*M_PI/180);
		GLFloat rho0 = 1025;
		GLFloat g = 9.81;
		
        // This is good for unit testing.
        BOOL shouldUnitTest = YES;
        NSUInteger modeUnit = 1;
		NSUInteger kUnit = 1;
		NSUInteger lUnit = 0;
		NSInteger omegaSign = 1;
        GLFloat U_max = .025;
        NSUInteger numStepsPerCycle = 4;
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
        GLMutableDimension *tDim = [[GLMutableDimension alloc] initWithPoints: @[@(0.0)]];
		tDim.name = @"time";
        
        /************************************************************************************************/
		/*		Create a density profile and compute the internal wave phases                           */
		/************************************************************************************************/
        GLFunction *z = [GLFunction functionOfRealTypeFromDimension:zDim withDimensions:@[zDim] forEquation:equation];
        GLFunction *rho_bar = [[z times: @(-N2*rho0/g)] plus: @(rho0)];
		rho_bar.name = @"rho_bar";
		
		// The ordering the dimensions is very deliberate here, for two reasons:
		// 1. The z-dimension is in the first position, so that the horizontal fft will act on contiguous chunks of memory, and
		// 2. The last two dimensions are ordered (x,y) to appease pcolor, meshgrid, and all the standard matlab formating.
        GLInternalWaveInitialization *wave = [[GLInternalWaveInitialization alloc] initWithDensityProfile: rho_bar fullDimensions:@[zDim, xDim, yDim] latitude:latitude equation:equation];
        
        if (shouldUnitTest) {
            wave.maximumModes = modeUnit+2;
            omega = [wave createUnitWaveWithSpeed: U_max verticalMode: modeUnit k: kUnit l: lUnit omegaSign: omegaSign];
        } else {
            wave.maximumModes = 2;
            [wave createGarrettMunkSpectrumWithEnergy: 1.0];
		}
        
		/************************************************************************************************/
		/*		Create the dynamical variables from the analytical solution								*/
		/************************************************************************************************/
		
		// We should check that we optimizations in places for these. They should be purely imaginary, and the multiplication and exponentiation should take advantage of that.
		GLFunction *iOmega = [[wave.eigenfrequencies swapComplex] makeHermitian];
		GLFunction *minusiOmega = [[[wave.eigenfrequencies swapComplex] negate] makeHermitian];

		NSArray * (^timeToUV) (GLScalar *) = ^( GLScalar *t ) {
			GLFunction *time_phase_plus = [[iOmega multiply: t] exponentiate];
			GLFunction *time_phase_minus = [[minusiOmega multiply: t] exponentiate];
			
            GLFunction *tmp = [wave.Sprime transform: [[wave.u_plus multiply: time_phase_plus] plus: [wave.u_minus multiply: time_phase_minus]]];

            
			GLFunction *u = [[wave.Sprime transform: [[wave.u_plus multiply: time_phase_plus] plus: [wave.u_minus multiply: time_phase_minus]]] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis), @(kGLDeltaBasis)]];
			GLFunction *v = [[wave.Sprime transform: [[wave.v_plus multiply: time_phase_plus] plus: [wave.v_minus multiply: time_phase_minus]]] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis), @(kGLDeltaBasis)]];
			GLFunction *w = [[wave.S transform: [[wave.w_plus multiply: time_phase_plus] plus: [wave.w_minus multiply: time_phase_minus]]] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis), @(kGLDeltaBasis)]];
            GLFunction *rho = [[[[wave.S transform: [[wave.rho_plus multiply: time_phase_plus] plus: [wave.rho_minus multiply: time_phase_minus]]] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis), @(kGLDeltaBasis)]] times: wave.N2] plus: rho_bar];
			
			return @[u,v,w,rho];
		};
        		
		GLScalar *t = [GLScalar scalarWithValue: 0.0*2*M_PI/f0 forEquation: equation];
		NSArray *uv = timeToUV(t);
		GLFunction *u = uv[0];
		GLFunction *v = uv[1];
		GLFunction *w = uv[2];
        GLFunction *rho = uv[3];
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
		
		CGFloat cfl = 0.25;
        GLFloat timeStep = cfl * xDim.sampleInterval / maxSpeed;
		timeStep = shouldUnitTest ? 2*M_PI/(numStepsPerCycle*omega*ceil(2*M_PI/(numStepsPerCycle*omega*timeStep))) : timeStep;
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
		
		NSString *path = [[NSSearchPathForDirectoriesInDomains(NSDesktopDirectory, NSUserDomainMask, YES) objectAtIndex:0] stringByAppendingPathComponent:@"InternalWaves.nc"];
		GLNetCDFFile *netcdfFile = [[GLNetCDFFile alloc] initWithURL: [NSURL URLWithString: path] forEquation: equation overwriteExisting: YES];
		
        [netcdfFile setGlobalAttribute: @(U_max) forKey: @"U_max"];
		[netcdfFile setGlobalAttribute: @(width) forKey: @"L_domain"];
        [netcdfFile setGlobalAttribute: @(latitude) forKey: @"latitude"];
        [netcdfFile setGlobalAttribute: @(f0) forKey: @"f0"];
        [netcdfFile setGlobalAttribute: @(depth) forKey: @"D"];
		
		[netcdfFile addVariable: rho_bar];
        
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
		
		GLFloat maxTime = shouldUnitTest ? maxWavePeriods*2*M_PI/omega : maxWavePeriods*2*M_PI/f0;
        GLFloat sampleTime = shouldUnitTest ? maxTime/numStepsPerCycle : sampleTimeInMinutes*60;
        for (GLFloat time = sampleTime; time <= maxTime-sampleTime; time += sampleTime)
		//while( integrator.currentTime < maxTime + integrator.stepSize/2)
        {
            @autoreleasepool {
                //NSArray *yout = [integrator stepForward];
				NSArray *yout = [integrator stepForwardToTime: time];
                
				GLFloat time = integrator.currentTime;
                NSLog(@"Logging day %d @ %02d:%02d (HH:MM), last step size: %02d:%02.1f (MM:SS.S).", (int) floor(time/86400), ((int) floor(time/3600))%24, ((int) floor(time/60))%60, (int)floor(integrator.lastStepSize/60), fmod(integrator.lastStepSize,60));
				NSLog(@"Logging day %d @ %02d:%02d (HH:MM), last step size: %02d:%02.1f (MM:SS.S).", (int) floor(integrator.currentTime/86400), ((int) floor(integrator.currentTime/3600))%24, ((int) floor(integrator.currentTime/60))%60, (int)floor(integrator.lastStepSize/60), fmod(integrator.lastStepSize,60));
				//NSLog(@"Logging day %d @ %02d:%02d (HH:MM)", (int) floor(time/86400), ((int) floor(time/3600))%24, ((int) floor(time/60))%60);

                
                NSArray *uv = timeToUV([GLScalar scalarWithValue: time forEquation: equation]);
				
				[tDim addPoint: @(time)];
				[uHistory concatenateWithLowerDimensionalVariable: uv[0] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                [vHistory concatenateWithLowerDimensionalVariable: uv[1] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
				[wHistory concatenateWithLowerDimensionalVariable: uv[2] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
				[rhoHistory concatenateWithLowerDimensionalVariable: uv[3] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
				
                [xPositionHistory concatenateWithLowerDimensionalVariable: yout[0] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                [yPositionHistory concatenateWithLowerDimensionalVariable: yout[1] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
				[zPositionHistory concatenateWithLowerDimensionalVariable: yout[2] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                
                [netcdfFile waitUntilAllOperationsAreFinished];
            }
        }
        
		[netcdfFile close];
	}
    return 0;
}

