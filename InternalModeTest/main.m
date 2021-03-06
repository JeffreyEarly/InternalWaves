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


int NOTmain(int argc, const char * argv[])
{
    @autoreleasepool {
		GLFloat N2_0 = 1.69e-4;
		GLFloat rho0 = 1025;
		GLFloat g = 9.81;
		GLFloat latitude = 33.0;
		GLFloat H = 300;
		
		GLEquation *equation = [[GLEquation alloc] init];
		GLDimension *zDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: 32 domainMin: -H length: H]; zDim.name = @"z";
		GLFunction *z = [GLFunction functionOfRealTypeFromDimension:zDim withDimensions:@[zDim] forEquation:equation];
		GLFunction *rho_bar = [[z times: @(-N2_0*rho0/g)] plus: @(rho0)];
		
		GLFloat width = 1000;
		GLFloat height = 500;
		NSUInteger Nx = 4;
		NSUInteger Ny = 4;
		GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: Nx domainMin: -width/2 length: width];
		xDim.name = @"x";
		GLDimension *yDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: Ny domainMin: -height/2 length: height];
		yDim.name = @"y";
		
		GLInternalModesSpectral *internalModes = [[GLInternalModesSpectral alloc] init];
		[internalModes internalWaveModesFromDensityProfile: rho_bar withFullDimensions:@[xDim, yDim, zDim] forLatitude: latitude maximumModes: 31];
    }
    
    return 0;
}

int main(int argc, const char * argv[])
{

	@autoreleasepool {
        
		GLFloat latitude = 33;
		GLFloat N2_0 = 1.69e-4;
		GLFloat depth = 5000;
		GLFloat width = 10e3;
        GLFloat height = 10e3;
		NSUInteger Nx = 64;
        NSUInteger Ny = 64;
		NSUInteger Nz = 64;
		
		GLFloat rho0 = 1025;
		GLFloat g = 9.81;
		
        /************************************************************************************************/
		/*		Define the problem dimensions															*/
		/************************************************************************************************/
		GLEquation *equation = [[GLEquation alloc] init];
		GLDimension *zDim = [[GLDimension alloc] initDimensionWithGrid: kGLChebyshevEndpointGrid nPoints: Nz domainMin: -depth length: depth];
		zDim.name = @"z";
		GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: Nx domainMin: -width/2 length: width];
		xDim.name = @"x";
		GLDimension *yDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: Ny domainMin: -height/2 length: height];
		yDim.name = @"y";
        GLMutableDimension *tDim = [[GLMutableDimension alloc] initWithPoints: @[@(0.0)]];
		tDim.name = @"time";
        GLDimension *zOutDim;
        
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
        }  else if (0) {
            GLNetCDFFile *profile = [[GLNetCDFFile alloc] initWithURL:[NSURL URLWithString: @"/Users/jearly/Documents/Models/InternalWaves/Latmix2011Site1Profile_AllPoints.nc"] forEquation:equation];
            GLFunction *rho_full = profile.variables[0];
            GLFunction *N2_full = profile.variables[1];
            GLDimension *zDim_full = rho_full.dimensions[0];
            
            // Interpolate the density onto the reduced grid (our z input grid).
            zDim = [[GLDimension alloc] initDimensionWithGrid: kGLChebyshevEndpointGrid nPoints: Nz domainMin: zDim_full.domainMin length: zDim_full.domainLength];
            GLFunction *z = [GLFunction functionOfRealTypeFromDimension:zDim withDimensions:@[zDim] forEquation:equation];
            rho_bar = [rho_full interpolateAtPoints:@[z]];
            [rho_bar solve];
            
            // Now create a stretched z-output grid
            GLFunction *N_scaled = [[N2_full dividedBy: [N2_full min]] sqrt];
            GLFunction *s = [N_scaled integrate]; // we now have s(z)
            
            NSUInteger Nz_out = 50;
            GLFloat minDepth = -60;
            GLFloat maxDepth = 0;
            GLFloat sAtMinDepth = 0.0;
            GLFloat sAtMaxDepth = 0.0;
            
            [s solve]; // Must solve before trying to use its data to initialize a dimension.
            GLDimension *sDim = [[GLDimension alloc] initWithNPoints: s.nDataPoints values: s.data];
            GLFunction *zOfs = [GLFunction functionOfRealTypeWithDimensions: @[sDim] forEquation: equation];
            for (NSUInteger i=0; i<zOfs.nDataPoints; i++) {
                zOfs.pointerValue[i] = [zDim_full valueAtIndex: i];
                if ([zDim_full valueAtIndex: i] <= minDepth) {
                    sAtMinDepth = s.pointerValue[i];
                }
                if ([zDim_full valueAtIndex: i] >= maxDepth && !sAtMaxDepth) {
                    sAtMaxDepth = s.pointerValue[i];
                }
            }
            
            GLDimension *sDimOut = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: Nz_out domainMin: sAtMinDepth length: sAtMaxDepth-sAtMinDepth];
            GLFunction *sOut = [GLFunction functionOfRealTypeFromDimension: sDimOut withDimensions:@[sDimOut] forEquation:equation];
            GLFunction *zInterp = [zOfs interpolateAtPoints:@[sOut]];
			
			[zInterp solve]; // required!!!
            zOutDim = [[GLDimension alloc] initWithNPoints: zInterp.nDataPoints values: zInterp.data];
		} else if (1) {
			GLFloat N0 = 5.23e-3;
			GLFunction *z = [GLFunction functionOfRealTypeFromDimension:zDim withDimensions:@[zDim] forEquation:equation];
			rho_bar = [[[[[[z times: @(2./1300.)] exponentiate] negate] plus: @(1.)] times: @(rho0*N0*N0*1300./(2*g))] plus: @(rho0)];
		} else {
            GLFunction *z = [GLFunction functionOfRealTypeFromDimension:zDim withDimensions:@[zDim] forEquation:equation];
            rho_bar = [[z times: @(-N2_0*rho0/g)] plus: @(rho0)];
        }
	    
		GLInternalModesSpectral *internalModes = [[GLInternalModesSpectral alloc] init];
 		//[internalModes internalGeostrophicModesFromDensityProfile: rho_bar forLatitude: latitude];
        //[internalModes internalWaveModesFromDensityProfile: rho_bar wavenumber: 1 forLatitude: latitude];
        //[internalModes internalWaveModesUsingGEPFromDensityProfile: rho_bar wavenumber: 0.008 forLatitude: latitude];
        
        //[internalModes internalWaveModesFromDensityProfile: rho_bar wavenumber: 0.0 forLatitude: latitude];
        [internalModes internalWaveModesFromDensityProfile: rho_bar wavenumber: 0.0 forLatitude: latitude maximumModes: 30 zOutDim: zOutDim];
        //[internalModes internalWaveModesFromDensityProfile: rho_bar withFullDimensions:@[xDim, yDim, zDim] forLatitude: latitude maximumModes: 30 zOutDim: zOutDim];
		
		[internalModes.N2 dumpToConsole];
		
        [internalModes.eigendepths dumpToConsole];
        [internalModes.eigenfrequencies dumpToConsole];
        [internalModes.S dumpToConsole];
        [internalModes.Sprime dumpToConsole];
        
		NSLog(@"Done!");
		
		//[internalModes.eigendepths dumpToConsole];
		
		return 0;
		
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

