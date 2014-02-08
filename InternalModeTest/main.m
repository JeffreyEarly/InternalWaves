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
		GLFloat N2_0 = 1e-4;
		GLFloat depth = 100;
		GLFloat width = 1000;
        GLFloat height = 500;
		NSUInteger Nx = 4;
        NSUInteger Ny = 4;
		NSUInteger Nz = 10;
		
		GLFloat f0 = 2*(7.2921e-5)*sin(latitude*M_PI/180);
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
        GLFunction *z = [GLFunction functionOfRealTypeFromDimension:zDim withDimensions:@[zDim] forEquation:equation];
        GLFunction *rho_bar = [[z times: @(-N2_0*rho0/g)] plus: @(rho0)];
		rho_bar.name = @"rho_bar";
	    
		GLInternalModes *internalModes = [[GLInternalModes alloc] init];
//        internalModes.maximumModes = 10;
 		//NSArray *system = [internalModes internalModesFromDensityProfile: rho_bar];
        //NSArray *system = [internalModes internalModesFromDensityProfile: rho_bar wavenumber: .1 latitude: latitude];
        NSArray *system = [internalModes internalWaveModesFromDensityProfile: rho_bar withFullDimensions:@[xDim, yDim, zDim ] forLatitude:latitude];
		
        GLFunction *eigendepths = system[0];
		GLLinearTransform *S = system[1];
		GLLinearTransform *Sprime = system[2];
		GLFunction *N2 = internalModes.N2;
        
		[eigendepths dumpToConsole];
//        GLFunction *L_r = [[internalModes.eigendepths times: @(g/(f0*f0))] sqrt];
//        [L_r dumpToConsole];
		[S dumpToConsole];
	}
    return 0;
}

