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
//        NSString *restartFile = [[NSSearchPathForDirectoriesInDomains(NSDesktopDirectory, NSUserDomainMask, YES) objectAtIndex:0] stringByAppendingPathComponent:@"InternalWavesLatmix2011_128_128_64_lat31.internalwaves"];
//        NSString *outputFile = @"InternalWavesLatmix2011_128_128_64_lat31_unit_test_no_diffusivity.nc";
        
        NSString *restartFile = [[NSSearchPathForDirectoriesInDomains(NSDesktopDirectory, NSUserDomainMask, YES) objectAtIndex:0] stringByAppendingPathComponent:@"InternalWavesConstantN_64_64_64_lat31.internalwaves"];
        NSString *outputFile = @"InternalWavesConstantN_64_64_64_lat31_unit_test_no_diffusivity.nc";
        
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
            GLFloat latitude = 31;
            GLFloat N2 = 1e-3; //2.5e-3;
            GLFloat width = 15e3;
            GLFloat height = 15e3;
            GLFloat depth = 300;
            NSUInteger Nx = 64;
            NSUInteger Ny = 64;
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
            }  else if (0) {
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
        
        wave.maximumModes = 60;
        //wave.maxDepth = -100;
        //[wave createGarrettMunkSpectrumWithEnergy: 0.125];
        [wave createUnitWaveWithSpeed: 0.01 verticalMode: 1 k: 1 l: 0 omegaSign: 1];
        zDim = wave.rho.dimensions[0];
        
        GLMutableDimension *tDim = [[GLMutableDimension alloc] initWithPoints: @[@(0.0)]];
        tDim.name = @"time";

		GLFloat maxWavePeriods = 2; // The wave period is the inertial period for the GM spectrum initialization, or omega for the unit test initialization
		GLFloat sampleTimeInMinutes = 15; // This will be overriden for the unit test.
		GLFloat horizontalFloatSpacingInMeters = 500;
        
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
        
        GLDimension *xFloatDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: ceil(4000/horizontalFloatSpacingInMeters) domainMin: -2000 length:4000];
		xFloatDim.name = @"x-float";
		GLDimension *yFloatDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: ceil(4000/horizontalFloatSpacingInMeters) domainMin: -2000 length:4000];
		yFloatDim.name = @"y-float";
		//GLDimension *zFloatDim = [[GLDimension alloc] initWithPoints: @[ @(-38), @(-31.5), @(-25)]];
        GLDimension *zFloatDim = [[GLDimension alloc] initWithPoints: @[ @(-32) ]];
		zFloatDim.name = @"z-float";
        
		// For consistency, we order the float dimensions the same as the dynamical variable dimensions.
        NSArray *floatDimensions = @[zFloatDim, xFloatDim, yFloatDim];
		GLFunction *xFloat = [GLFunction functionOfRealTypeFromDimension: xFloatDim withDimensions: floatDimensions forEquation: wave.equation];
		GLFunction *yFloat = [GLFunction functionOfRealTypeFromDimension: yFloatDim withDimensions: floatDimensions forEquation: wave.equation];
		GLFunction *zFloat = [GLFunction functionOfRealTypeFromDimension: zFloatDim withDimensions: floatDimensions forEquation: wave.equation];
        
		GLFunction *xIsopycnal = [GLFunction functionFromFunction: xFloat];
		GLFunction *yIsopycnal = [GLFunction functionFromFunction: yFloat];
		GLFunction *zIsopycnal = [GLFunction functionFromFunction: zFloat];
        GLFunction *isopycnalDeviation = [zeta interpolateAtPoints:@[zIsopycnal, xIsopycnal, yIsopycnal]];
        zIsopycnal = [zIsopycnal plus: isopycnalDeviation];
        
        GLFunction *xIsopycnalDiffusive = [GLFunction functionFromFunction: xIsopycnal];
		GLFunction *yIsopycnalDiffusive = [GLFunction functionFromFunction: yIsopycnal];
		GLFunction *zIsopycnalDiffusive = [GLFunction functionFromFunction: zIsopycnal];
        
        GLFunction *xFixedDepth = [GLFunction functionFromFunction: xFloat];
		GLFunction *yFixedDepth = [GLFunction functionFromFunction: yFloat];
		GLFunction *zFixedDepth = [GLFunction functionFromFunction: zFloat];
        
        GLFunction *xDrifter = [GLFunction functionFromFunction: xFloat];
		GLFunction *yDrifter = [GLFunction functionFromFunction: yFloat];
		GLFunction *zDrifter = [GLFunction functionFromFunction: zFloat];
        
		
		CGFloat cfl = 0.25;
        GLFloat cflTimeStep = cfl * xDim.sampleInterval / maxSpeed;
		GLFloat outputTimeStep = sampleTimeInMinutes*60;
		GLFloat timeStep = cflTimeStep > outputTimeStep ? outputTimeStep : outputTimeStep / ceil(outputTimeStep/cflTimeStep);
		
        GLFloat drogueMin = -33;
        GLFloat drogueMax = -27;
        NSUInteger drogueMinIndex = NSNotFound;
        NSUInteger drogueMaxIndex = NSNotFound;
        for (NSUInteger iPoint=0; iPoint < zDim.nPoints; iPoint++) {
            if ([zDim valueAtIndex: iPoint] < drogueMin ) drogueMinIndex=iPoint;
            if ([zDim valueAtIndex: iPoint] < drogueMax ) drogueMaxIndex=iPoint;
        }
        NSRange drogueRange = NSMakeRange(drogueMaxIndex, drogueMaxIndex-drogueMinIndex+1);
        
        GLFloat kappa = 5e-6; // m^2/s
        GLFloat norm = sqrt(timeStep*2*kappa);
        norm = sqrt(4)*norm/timeStep; // the integrator multiplies by deltaT, so we account for that here.
        // RK4: dt/3 f(0) + dt/6 f(1) + dt/6 *f(4) + dt/3*f(3)
        // Mean of 1/3 and 1/6? 1/4. It's  the geometric mean to get the same norm, hence, sqrt 4.
        
        NSArray *y=@[xIsopycnal, yIsopycnal, zIsopycnal, xIsopycnalDiffusive, yIsopycnalDiffusive, zIsopycnalDiffusive, xFixedDepth, yFixedDepth, xDrifter, yDrifter];
        GLAdaptiveRungeKuttaOperation *integrator = [GLAdaptiveRungeKuttaOperation rungeKutta4AdvanceY: y stepSize: timeStep fFromTY:^(GLScalar *time, NSArray *yNew) {
			NSArray *uv = timeToUV(time);
            GLFunction *xStep = [GLFunction functionWithNormallyDistributedValueWithDimensions: floatDimensions forEquation: wave.equation];
            GLFunction *yStep = [GLFunction functionWithNormallyDistributedValueWithDimensions: floatDimensions forEquation: wave.equation];
            GLFunction *zStep = [GLFunction functionWithNormallyDistributedValueWithDimensions: floatDimensions forEquation: wave.equation];
            xStep = [xStep times: @(norm)];
            yStep = [yStep times: @(norm)];
            zStep = [zStep times: @(norm)];
            
			GLSimpleInterpolationOperation *interpIso = [[GLSimpleInterpolationOperation alloc] initWithFirstOperand: @[uv[0], uv[1], uv[2]] secondOperand: @[yNew[2], yNew[0], yNew[1]]];
            NSMutableArray *f = [interpIso.result mutableCopy];
            
            GLSimpleInterpolationOperation *interpIsoDiff = [[GLSimpleInterpolationOperation alloc] initWithFirstOperand: @[uv[0], uv[1], uv[2]] secondOperand: @[yNew[5], yNew[3], yNew[4]]];
            NSArray *f2 = @[[interpIsoDiff.result[0] plus: xStep], [interpIsoDiff.result[1] plus: yStep], [interpIsoDiff.result[2] plus: zStep]];
            [f addObjectsFromArray: f2];
            
            GLSimpleInterpolationOperation *interpFixedDepth = [[GLSimpleInterpolationOperation alloc] initWithFirstOperand: @[uv[0], uv[1]] secondOperand: @[zFixedDepth, yNew[6], yNew[7]]];
            NSArray *f3 = @[[interpFixedDepth.result[0] plus: xStep], [interpFixedDepth.result[1] plus: yStep]];
            [f addObjectsFromArray: f3];
            
            GLFunction *uMean = [uv[0] mean: 0 range: drogueRange];
            GLFunction *vMean = [uv[1] mean: 0 range: drogueRange];
            GLSimpleInterpolationOperation *interpDrifter = [[GLSimpleInterpolationOperation alloc] initWithFirstOperand: @[uMean, vMean] secondOperand: @[yNew[8], yNew[9]]];
            NSArray *f4 = @[[interpDrifter.result[0] plus: xStep], [interpDrifter.result[1] plus: yStep]];
            [f addObjectsFromArray: f4];
            
			return f;
		}];
        
        /************************************************************************************************/
		/*		Create a NetCDF file and mutable variables in order to record some of the time steps.	*/
		/************************************************************************************************/
		
		NSString *path = [[NSSearchPathForDirectoriesInDomains(NSDesktopDirectory, NSUserDomainMask, YES) objectAtIndex:0] stringByAppendingPathComponent:outputFile];
		GLNetCDFFile *netcdfFile = [[GLNetCDFFile alloc] initWithURL: [NSURL URLWithString: path] forEquation: wave.equation overwriteExisting: YES];
		
        [netcdfFile setGlobalAttribute: @(wave.f0) forKey: @"f0"];
		
		[netcdfFile addVariable: wave.rho];
        [netcdfFile addVariable: wave.N2];
        
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
		
		GLMutableVariable *xPositionHistory = [xIsopycnal variableByAddingDimension: tDim];
		xPositionHistory.name = @"x-position";
		xPositionHistory = [netcdfFile addVariable: xPositionHistory];
        
		GLMutableVariable *yPositionHistory = [yIsopycnal variableByAddingDimension: tDim];
		yPositionHistory.name = @"y-position";
		yPositionHistory = [netcdfFile addVariable: yPositionHistory];
		
		GLMutableVariable *zPositionHistory = [zIsopycnal variableByAddingDimension: tDim];
		zPositionHistory.name = @"z-position";
		zPositionHistory = [netcdfFile addVariable: zPositionHistory];
        
        GLMutableVariable *xIsopycnalDiffusiveHistory = [xIsopycnalDiffusive variableByAddingDimension: tDim];
		xIsopycnalDiffusiveHistory.name = @"x-position-diffusive";
		xIsopycnalDiffusiveHistory = [netcdfFile addVariable: xIsopycnalDiffusiveHistory];
        
		GLMutableVariable *yIsopycnalDiffusiveHistory = [yIsopycnalDiffusive variableByAddingDimension: tDim];
		yIsopycnalDiffusiveHistory.name = @"y-position-diffusive";
		yIsopycnalDiffusiveHistory = [netcdfFile addVariable: yIsopycnalDiffusiveHistory];
		
		GLMutableVariable *zIsopycnalDiffusiveHistory = [zIsopycnalDiffusive variableByAddingDimension: tDim];
		zIsopycnalDiffusiveHistory.name = @"z-position-diffusive";
		zIsopycnalDiffusiveHistory = [netcdfFile addVariable: zIsopycnalDiffusiveHistory];
        
        GLMutableVariable *xFixedDepthHistory = [xFixedDepth variableByAddingDimension: tDim];
		xFixedDepthHistory.name = @"x-position-fixed-depth";
		xFixedDepthHistory = [netcdfFile addVariable: xFixedDepthHistory];
        
		GLMutableVariable *yFixedDepthHistory = [yFixedDepth variableByAddingDimension: tDim];
		yFixedDepthHistory.name = @"y-position-fixed-depth";
		yFixedDepthHistory = [netcdfFile addVariable: yFixedDepthHistory];
		
		GLMutableVariable *zFixedDepthHistory = [zFixedDepth variableByAddingDimension: tDim];
		zFixedDepthHistory.name = @"z-position-fixed-depth";
		zFixedDepthHistory = [netcdfFile addVariable: zFixedDepthHistory];
        
        GLMutableVariable *xDrifterHistory = [xDrifter variableByAddingDimension: tDim];
		xDrifterHistory.name = @"x-position-drifter";
		xDrifterHistory = [netcdfFile addVariable: xDrifterHistory];
        
		GLMutableVariable *yDrifterHistory = [yDrifter variableByAddingDimension: tDim];
		yDrifterHistory.name = @"y-position-drifter";
		yDrifterHistory = [netcdfFile addVariable: yDrifterHistory];
		
		GLMutableVariable *zDrifterHistory = [zDrifter variableByAddingDimension: tDim];
		zDrifterHistory.name = @"z-position-drifter";
		zDrifterHistory = [netcdfFile addVariable: zDrifterHistory];
		
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
				NSArray *yout = [integrator stepForwardToTime: time];
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
				
                [xPositionHistory concatenateWithLowerDimensionalVariable: yout[0] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                [yPositionHistory concatenateWithLowerDimensionalVariable: yout[1] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
				[zPositionHistory concatenateWithLowerDimensionalVariable: yout[2] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                
                [xIsopycnalDiffusiveHistory concatenateWithLowerDimensionalVariable: yout[3] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                [yIsopycnalDiffusiveHistory concatenateWithLowerDimensionalVariable: yout[4] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
				[zIsopycnalDiffusiveHistory concatenateWithLowerDimensionalVariable: yout[5] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                
                [xFixedDepthHistory concatenateWithLowerDimensionalVariable: yout[6] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                [yFixedDepthHistory concatenateWithLowerDimensionalVariable: yout[7] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
				[zFixedDepthHistory concatenateWithLowerDimensionalVariable: zFixedDepth alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
            
                [xDrifterHistory concatenateWithLowerDimensionalVariable: yout[8] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                [yDrifterHistory concatenateWithLowerDimensionalVariable: yout[9] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
				[zDrifterHistory concatenateWithLowerDimensionalVariable: zDrifter alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                
                [netcdfFile waitUntilAllOperationsAreFinished];
            }
        }
        
		[netcdfFile close];
	}
    return 0;
}