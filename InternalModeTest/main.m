//
//  main.m
//  InternalModeTest
//
//  Created by Jeffrey J. Early on 2/7/14.
//  Copyright (c) 2014 Jeffrey J. Early. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>
#import "GLInternalModes.h"

int main(int argc, const char * argv[])
{

	@autoreleasepool {
		GLFloat latitude = 45;
		GLFloat N2_0 = 1e-2;
		GLFloat depth = 100;
		GLFloat width = 1000;
        GLFloat height = 500;
		NSUInteger Nx = 4;
        NSUInteger Ny = 4;
		NSUInteger Nz = 256;
		
		GLFloat rho0 = 1025;
		GLFloat g = 9.81;
		
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
        GLFunction *rho_bar;
        if (0) {
            GLNetCDFFile *profile = [[GLNetCDFFile alloc] initWithURL:[NSURL URLWithString: @"/Users/jearly/Documents/Models/InternalWaves/Latmix2011Site1Profile.nc"] forEquation:equation];
            GLFunction *rho_profile = profile.variables[0];
            GLFunction *z = [GLFunction functionOfRealTypeFromDimension:zDim withDimensions:@[zDim] forEquation:equation];
            rho_bar = [rho_profile interpolateAtPoints:@[z]];
        } else {
            GLFunction *z = [GLFunction functionOfRealTypeFromDimension:zDim withDimensions:@[zDim] forEquation:equation];
            rho_bar = [[z times: @(-N2_0*rho0/g)] plus: @(rho0)];
        }
	    
		GLInternalModes *internalModes = [[GLInternalModes alloc] init];
 		//[internalModes internalGeostrophicModesFromDensityProfile: rho_bar forLatitude: latitude];
        //[internalModes internalWaveModesFromDensityProfile: rho_bar wavenumber: .01 forLatitude: latitude];
        [internalModes internalWaveModesFromDensityProfile: rho_bar withFullDimensions:@[xDim, yDim, zDim] forLatitude: latitude];
        
		[internalModes.eigendepths dumpToConsole];
        [internalModes.eigenfrequencies dumpToConsole];
		[internalModes.Sprime dumpToConsole];
        [internalModes.S dumpToConsole];
	}
    return 0;
}

