//
//  main.m
//  InternalModeTest
//
//  Created by Jeffrey J. Early on 2/7/14.
//  Copyright (c) 2014 Jeffrey J. Early. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>
#import <GLOceanKit/GLOceanKit.h>
#import <GLOceanKit/GLInternalModesSpectral.h>

int main(int argc, const char * argv[])
{

	@autoreleasepool {
		GLFloat latitude = 33;
		GLFloat N2_0 = 1.69e-4;
		GLFloat depth = 300;
		GLFloat width = 1000;
        GLFloat height = 500;
		NSUInteger Nx = 4;
        NSUInteger Ny = 4;
		NSUInteger Nz = 32;
		
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
        } else if (0) {
            GLNetCDFFile *profile = [[GLNetCDFFile alloc] initWithURL:[NSURL URLWithString: @"/Users/jearly/Documents/Models/InternalWaves/Latmix2011Site1Profile_Stretched_64.nc"] forEquation:equation];
            rho_bar = profile.variables[0];
            zDim = rho_bar.dimensions[0];
        } else {
            GLFunction *z = [GLFunction functionOfRealTypeFromDimension:zDim withDimensions:@[zDim] forEquation:equation];
            rho_bar = [[z times: @(-N2_0*rho0/g)] plus: @(rho0)];
        }
	    
		GLInternalModes *internalModes = [[GLInternalModes alloc] init];
 		//[internalModes internalGeostrophicModesFromDensityProfile: rho_bar forLatitude: latitude];
        //[internalModes internalWaveModesFromDensityProfile: rho_bar wavenumber: 1 forLatitude: latitude];
        //[internalModes internalWaveModesUsingGEPFromDensityProfile: rho_bar wavenumber: 0.008 forLatitude: latitude];
        [internalModes internalWaveModesFromDensityProfile: rho_bar withFullDimensions:@[xDim, yDim, zDim] forLatitude: latitude];
        
		[internalModes.eigendepths dumpToConsole];
//        [internalModes.eigenfrequencies dumpToConsole];
//		[internalModes.Sprime dumpToConsole];
//        [internalModes.S dumpToConsole];
        
        GLLinearTransform *op = [internalModes.diffOp times: internalModes.S];
        //[op dumpToConsole];
        return 0;
        NSUInteger nRows = internalModes.S.matrixDescription.strides[0].nRows;
        NSUInteger rowStride = internalModes.S.matrixDescription.strides[0].rowStride;
        NSUInteger colStride = internalModes.S.matrixDescription.strides[0].columnStride;
        GLFloat dz = zDim.sampleInterval;
        for (NSUInteger iMode=0; iMode<10; iMode++) {
            GLFloat gSum = 0.0;
            GLFloat fSum = 0.0;
            GLFloat gSum2 = 0.0;
            GLFloat omega = internalModes.eigenfrequencies.pointerValue[iMode];
            for (NSUInteger i=0; i<nRows; i++) {
                GLFloat F = internalModes.Sprime.pointerValue[i*rowStride + iMode*colStride];
                GLFloat G = internalModes.S.pointerValue[i*rowStride + iMode*colStride];
                GLFloat G2 = op ? op.pointerValue[i*rowStride + iMode*colStride] : 0.0;
                GLfloat N2 = internalModes.N2.pointerValue[i];
                
                gSum += (1/g)*(N2-internalModes.f0*internalModes.f0)*G*G;
                gSum2 += (1/g)*(N2-internalModes.f0*internalModes.f0)*G2*G2;
                fSum += fabs((N2-internalModes.f0*internalModes.f0)/(N2-omega*omega))*F*F;
                //fSum += F*F;
            }
            gSum *= dz;
            fSum *= dz;
            gSum2 *= dz;
            GLFloat h = internalModes.eigendepths.pointerValue[iMode];
            
            GLFloat Nmin = sqrt(fabs([internalModes.N2 minNow]));
            printf("G: %g, G2_h: %g, F: %g, h: %g, omega: %g, N_min: %g\n", gSum, 1/sqrt(gSum2), fSum, h, omega, Nmin);
        }
        
//        GLFunction *N = [internalModes.N2 sqrt];
//        [N dumpToConsole];
	}
    return 0;
}

