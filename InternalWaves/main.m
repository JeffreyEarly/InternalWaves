//
//  main.m
//  InternalWaves
//
//  Created by Jeffrey J. Early on 12/9/13.
//  Copyright (c) 2013 Jeffrey J. Early. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>

@interface GLInternalModes : NSObject

@end

@implementation GLInternalModes

- (NSArray *) internalModesFromDensityProfile: (GLFunction *) rho
{
    if (rho.dimensions.count != 1) {
        [NSException raise:@"InvalidDimensions" format:@"Only one dimension allowed, at this point"];
    }
    
    GLEquation *equation = rho.equation;
    GLDimension *zDim = rho.dimensions[0];
    
    GLLinearTransform *diffZ = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 1 leftBC: kGLNeumannBoundaryCondition rightBC:kGLNeumannBoundaryCondition bandwidth:1 fromDimension:zDim forEquation:equation];
    GLFunction *N2 = [diffZ transform: [[rho dividedBy: [rho mean]] times: @(1)]];
	
    GLFunction *invN2 = [N2 scalarDivide: 1.0];
		
    GLLinearTransform *diffZZ = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 2 leftBC: kGLDirichletBoundaryCondition rightBC:kGLDirichletBoundaryCondition bandwidth:1 fromDimension:zDim forEquation:equation];
		
    GLLinearTransform *invN2_trans = [GLLinearTransform linearTransformFromFunction: invN2];
    
    GLLinearTransform *diffOp = [invN2_trans multiply: diffZZ];
	
    NSArray *system = [diffOp eigensystem];
    
    return system;
}

- (NSArray *) internalModesFromDensityProfile: (GLFunction *) rho wavenumber: (GLFloat) k latitude: (GLFloat) latitude
{
    if (rho.dimensions.count != 1) {
        [NSException raise:@"InvalidDimensions" format:@"Only one dimension allowed, at this point"];
    }
	
	GLFloat f0 = 2*(7.2921e-5)*sin(latitude*M_PI/180);
	GLFloat g = 9.81;
    
    GLEquation *equation = rho.equation;
    GLDimension *zDim = rho.dimensions[0];
    
	// First construct N^2
    GLLinearTransform *diffZ = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 1 leftBC: kGLNeumannBoundaryCondition rightBC:kGLNeumannBoundaryCondition bandwidth:1 fromDimension:zDim forEquation:equation];
    GLFunction *N2 = [diffZ transform: [[rho dividedBy: [rho mean]] times: @(-g)]];
	
	// Now construct A = k*k*eye(N) - Diff2;
	GLLinearTransform *k2 = [GLLinearTransform transformOfType: kGLRealDataFormat withFromDimensions: @[zDim] toDimensions: @[zDim] inFormat: @[@(kGLDiagonalMatrixFormat)] forEquation:equation matrix:^( NSUInteger *row, NSUInteger *col ) {
		return (GLFloatComplex) (row[0]==col[0] ? k*k : 0);
	}];
    GLLinearTransform *diffZZ = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 2 leftBC: kGLDirichletBoundaryCondition rightBC:kGLDirichletBoundaryCondition bandwidth:1 fromDimension:zDim forEquation:equation];
	GLLinearTransform *A = [k2 minus: diffZZ];
    
	// Now construct B = k*k*diag(N2) - f0*f0*Diff2;
	GLLinearTransform *B = [[[GLLinearTransform linearTransformFromFunction: N2] times: @(k*k)] minus: [diffZZ times: @(f0*f0)]];

    NSArray *system = [B generalizedEigensystemWith: A];
    
    return system;
}

- (NSArray *) internalWaveModesFromDensityProfile: (GLFunction *) rho  withHorizontalDimensions: (NSArray *) horizDims forLatitude: (GLFloat) latitude
{
	// rho will be a function of z only, so we need to create the full dimensions.
	NSMutableArray *dimensions = [rho.dimensions mutableCopy];
	[dimensions addObjectsFromArray: horizDims];
	
	// create an array with the intended transformation (this is agnostic to dimension ordering).
	NSMutableArray *basis = [NSMutableArray array];
	GLDimension *zDim;
	for (GLDimension *dim in dimensions) {
		if ( [dim.name isEqualToString: @"x"] || [dim.name isEqualToString: @"y"]) {
			[basis addObject: @(kGLExponentialBasis)];
		} else {
			zDim = dim;
			[basis addObject: @(dim.basisFunction)];
		}
	}
	
	NSArray *transformedDimensions = [GLDimension dimensionsForRealFunctionWithDimensions: dimensions transformedToBasis: basis];
	GLDimension *kDim, *lDim;
	for (GLDimension *dim in transformedDimensions) {
		if ( [dim.name isEqualToString: @"k"]) {
			kDim = dim;
		} else if ( [dim.name isEqualToString: @"l"]) {
			lDim = dim;
		}
	}
		
	GLEquation *equation = rho.equation;
	GLFunction *k = [GLFunction functionOfRealTypeFromDimension: kDim withDimensions: transformedDimensions forEquation: equation];
	GLFunction *l = [GLFunction functionOfRealTypeFromDimension: lDim withDimensions: transformedDimensions forEquation: equation];
	GLFunction *K2 = [[k multiply: k] plus: [l multiply: l]];
	
	GLFloat f0 = 2*(7.2921e-5)*sin(latitude*M_PI/180);
    
	// First construct N^2
    GLLinearTransform *diffZ1D = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 1 leftBC: kGLNeumannBoundaryCondition rightBC:kGLNeumannBoundaryCondition bandwidth:1 fromDimension:zDim forEquation:equation];
    GLFunction *N2 = [diffZ1D transform: [[rho dividedBy: [rho mean]] times: @(1)]];

	// Now construct A = k*k*eye(N) - Diff2;
    GLLinearTransform *diffZZ1D = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 2 leftBC: kGLDirichletBoundaryCondition rightBC:kGLDirichletBoundaryCondition bandwidth:1 fromDimension:zDim forEquation:equation];
	GLLinearTransform *diffZZ = [diffZZ1D expandedWithFromDimensions: transformedDimensions toDimensions: transformedDimensions];
	GLLinearTransform *A = [[GLLinearTransform linearTransformFromFunction:K2] minus: diffZZ];

	// Now construct B = k*k*diag(N2) - f0*f0*Diff2;
	GLLinearTransform *B = [[GLLinearTransform linearTransformFromFunction: [K2 multiply: N2]] minus: [diffZZ times: @(f0*f0)]];
	
    NSArray *system = [B generalizedEigensystemWith: A];
    
    return system;
}

@end


int main(int argc, const char * argv[])
{

	@autoreleasepool {
		GLFloat latitude = 45;
		GLFloat N2 = 1e-4;
		GLFloat depth = 100;
		GLFloat width = 1000;
		NSUInteger Nx = 8;
		NSUInteger Nz = 64;
		
		GLFloat f0 = 2*(7.2921e-5)*sin(latitude*M_PI/180);
		GLFloat rho0 = 1025;
		GLFloat g = 9.81;
		
		GLEquation *equation = [[GLEquation alloc] init];
		GLDimension *zDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: Nz domainMin: -depth length: depth];
		zDim.name = @"z";
		GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: Nx domainMin: -width/2 length: width];
		xDim.name = @"x";
		GLDimension *yDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: Nx domainMin: -width/2 length: width];
		yDim.name = @"y";
		
		// We create the z variable with dimensions in reverse order so that the fft will act on contiguous chunks of memory.
        GLFunction *z = [GLFunction functionOfRealTypeFromDimension:zDim withDimensions:@[zDim] forEquation:equation];
        GLFunction *rho = [[z times: @(-N2*rho0/g)] plus: @(rho0)];

        GLInternalModes *internalModes = [[GLInternalModes alloc] init];
        //NSArray *system = [internalModes internalModesFromDensityProfile: rho];
		//NSArray *system = [internalModes internalWaveModesFromDensityProfile: rho withHorizontalDimensions: @[yDim, xDim] forLatitude:latitude];
		NSArray *system = [internalModes internalModesFromDensityProfile: rho wavenumber: 0 latitude: 45.0];
        GLFunction *eigenvalues = system[0];
		GLLinearTransform *S = system[1];
        
		eigenvalues = [[[eigenvalues times: @(1/(f0*f0))] abs] sqrt];
        [eigenvalues dumpToConsole];
		[S dumpToConsole];
	}
    return 0;
}

