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
    GLScalar *rho_0 = [rho mean: 0];
    GLFunction *N2 = [diffZ transform: [rho dividedBy: rho_0]];
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
    GLScalar *rho_0 = [rho mean: 0];
    GLFunction *N2 = [diffZ transform: [rho dividedBy: rho_0]];
	
	// Now construct A = k*k*eye(N) - Diff2;
	GLLinearTransform *k2 = [GLLinearTransform transformOfType: kGLRealDataFormat withFromDimensions: @[zDim] toDimensions: @[zDim] inFormat: @[@(kGLDiagonalMatrixFormat)] forEquation:equation matrix:^( NSUInteger *row, NSUInteger *col ) {
		return (GLFloatComplex) (row[0]==col[0] ? k*k : 0);
	}];
    GLLinearTransform *diffZZ = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 2 leftBC: kGLDirichletBoundaryCondition rightBC:kGLDirichletBoundaryCondition bandwidth:1 fromDimension:zDim forEquation:equation];
	GLLinearTransform *A = [k2 minus: diffZZ];
    
	// Now construct B = k*k*diag(N2) - f0*f0*Diff2;
	GLLinearTransform *B = [[[GLLinearTransform linearTransformFromFunction: N2] times: @(k*k)] minus: [diffZZ times: @(f0*f0)]];

    NSArray *system = [A generalizedEigensystemWith: B];
    
    return system;
}

@end


int main(int argc, const char * argv[])
{

	@autoreleasepool {
	    
		GLEquation *equation = [[GLEquation alloc] init];
		NSUInteger N=32;
		GLDimension *zDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: N domainMin: 0.0 length: N-1];
		zDim.name = @"z";
        GLFunction *z = [GLFunction functionOfRealTypeFromDimension:zDim withDimensions:@[zDim] forEquation:equation];
        
        GLFunction *rho = [z times: @(1)];
        GLInternalModes *internalModes = [[GLInternalModes alloc] init];
        //NSArray *system = [internalModes internalModesFromDensityProfile: rho];
		NSArray *system = [internalModes internalModesFromDensityProfile: rho wavenumber: 1.0 latitude: 45.0];
        GLFunction *eigenvalues = system[0];
		GLLinearTransform *S = system[1];
        
        [eigenvalues dumpToConsole];
		[S dumpToConsole];
	}
    return 0;
}

