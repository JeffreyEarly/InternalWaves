//
//  main.m
//  InternalWaves
//
//  Created by Jeffrey J. Early on 12/9/13.
//  Copyright (c) 2013 Jeffrey J. Early. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>

int main(int argc, const char * argv[])
{

	@autoreleasepool {
	    
		GLEquation *equation = [[GLEquation alloc] init];
		NSUInteger N=32;
		GLDimension *zDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: N domainMin: 0.0 length: N-1];
		zDim.name = @"z";
		
		GLLinearTransform *diffZZ = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 2 leftBC: kGLDirichletBoundaryCondition rightBC:kGLDirichletBoundaryCondition bandwidth:1 fromDimension:zDim forEquation:equation];
		diffZZ = [diffZZ densified];
		[diffZZ dumpToConsole];
		
		NSArray *system = [diffZZ eigensystem];
		GLFunction *eigenvalues = system[0];
		GLLinearTransform *S = system[1];
		
		[eigenvalues dumpToConsole];
		[S dumpToConsole];
	}
    return 0;
}

