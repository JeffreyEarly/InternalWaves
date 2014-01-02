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
    
    GLLinearTransform *diffZ = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 1 leftBC: kGLDirichletBoundaryCondition rightBC:kGLDirichletBoundaryCondition bandwidth:1 fromDimension:zDim forEquation:equation];
    
    GLLinearTransform *diffZZ = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 2 leftBC: kGLDirichletBoundaryCondition rightBC:kGLDirichletBoundaryCondition bandwidth:1 fromDimension:zDim forEquation:equation];
    
    GLScalar *rho_0 = [rho mean: 0];
    GLFunction *N2 = [diffZ transform: [rho dividedBy: rho_0]];
    
    GLFunction *invN2 = [N2 scalarDivide: 1.0];
    GLLinearTransform *invN2_trans = [GLLinearTransform linearTransformFromFunction: invN2];
    
    GLFunction *invN2_z = [diffZ transform: invN2];
    GLLinearTransform *invN2_z_trans = [GLLinearTransform linearTransformFromFunction: invN2_z];
    
    GLLinearTransform *diffOp = [[invN2_z_trans multiply: diffZ] plus: [invN2_trans multiply: diffZZ]];
    NSArray *system = [diffOp eigensystem];
    
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
		
		GLLinearTransform *diffZZ = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 2 leftBC: kGLDirichletBoundaryCondition rightBC:kGLDirichletBoundaryCondition bandwidth:1 fromDimension:zDim forEquation:equation];
		diffZZ = [diffZZ densified];
		[diffZZ dumpToConsole];
		
		NSArray *system = [diffZZ eigensystem];
		GLFunction *eigenvalues = system[0];
		GLLinearTransform *S = system[1];
		
		[eigenvalues dumpToConsole];
		[S dumpToConsole];
        
        GLFunction *rho = [z times: @(1)];
        GLInternalModes *internalModes = [[GLInternalModes alloc] init];
        system = [internalModes internalModesFromDensityProfile: rho];
        eigenvalues = system[0];
		S = system[1];
        
        [eigenvalues dumpToConsole];
		[S dumpToConsole];
	}
    return 0;
}

