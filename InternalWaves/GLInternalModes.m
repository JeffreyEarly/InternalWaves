//
//  GLInternalModes.m
//  InternalWaves
//
//  Created by Jeffrey J. Early on 1/14/14.
//  Copyright (c) 2014 Jeffrey J. Early. All rights reserved.
//

#import "GLInternalModes.h"

@implementation GLInternalModes

- (NSArray *) internalModesFromDensityProfile: (GLFunction *) rho
{
    if (rho.dimensions.count != 1) {
        [NSException raise:@"InvalidDimensions" format:@"Only one dimension allowed, at this point"];
    }
    
	GLFloat g = 9.81;
    GLScalar *rho0 = [rho mean];
	
    GLEquation *equation = rho.equation;
    GLDimension *zDim = rho.dimensions[0];
    
	// First construct N^2
    GLLinearTransform *diffZ = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 1 leftBC: kGLNeumannBoundaryCondition rightBC:kGLNeumannBoundaryCondition bandwidth:1 fromDimension:zDim forEquation:equation];
    self.N2 = [diffZ transform: [[rho dividedBy: rho0] times: @(-g)]];
	
	
    GLFunction *invN2 = [self.N2 scalarDivide: -g];
	
    GLLinearTransform *diffZZ = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 2 leftBC: kGLDirichletBoundaryCondition rightBC:kGLDirichletBoundaryCondition bandwidth:1 fromDimension:zDim forEquation:equation];
    GLLinearTransform *invN2_trans = [GLLinearTransform linearTransformFromFunction: invN2];
    GLLinearTransform *diffOp = [invN2_trans multiply: diffZZ];
	
    NSArray *system = [diffOp eigensystemWithOrder: NSOrderedAscending];
	
	GLFunction *lambda = [system[0] makeRealIfPossible];
    GLLinearTransform *S = [system[1] makeRealIfPossible];
	
	if (self.maximumModes) {
		S = [S reducedFromDimensions: [NSString stringWithFormat: @"0:%lu", self.maximumModes-1] toDimension: @":"];
		lambda = [lambda variableFromIndexRangeString:[NSString stringWithFormat: @"0:%lu", self.maximumModes-1]];
	}
	
	self.eigendepths = [lambda scalarDivide: 1.0];
    self.eigenfrequencies = [self.eigendepths times: @(0)];
    self.S = [S normalizeWithFunction: [self.N2 times: @(1/g)]];
	self.Sprime = [diffZ multiply: S];
	
    return @[self.eigendepths, self.S, self.Sprime];
}

- (NSArray *) internalModesFromDensityProfile: (GLFunction *) rho wavenumber: (GLFloat) k latitude: (GLFloat) latitude
{
    if (rho.dimensions.count != 1) {
        [NSException raise:@"InvalidDimensions" format:@"Only one dimension allowed, at this point"];
    }
    
    GLFloat f0 = 2*(7.2921e-5)*sin(latitude*M_PI/180);
	GLFloat g = 9.81;
    GLScalar *rho0 = [rho mean];
	
    GLEquation *equation = rho.equation;
    GLDimension *zDim = rho.dimensions[0];
    
	// First construct N^2
    GLLinearTransform *diffZ = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 1 leftBC: kGLNeumannBoundaryCondition rightBC:kGLNeumannBoundaryCondition bandwidth:1 fromDimension:zDim forEquation:equation];
    self.N2 = [diffZ transform: [[rho dividedBy: rho0] times: @(-g)]];
	
	
    GLFunction *invN2 = [[self.N2 minus: @(f0*f0)] scalarDivide: -g];
	
    GLLinearTransform *diffZZ = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 2 leftBC: kGLDirichletBoundaryCondition rightBC:kGLDirichletBoundaryCondition bandwidth:1 fromDimension:zDim forEquation:equation];
    GLLinearTransform *invN2_trans = [GLLinearTransform linearTransformFromFunction: invN2];
    GLLinearTransform *k2 = [GLLinearTransform transformOfType: kGLRealDataFormat withFromDimensions: @[zDim] toDimensions: @[zDim] inFormat: @[@(kGLDiagonalMatrixFormat)] forEquation:equation matrix:^( NSUInteger *row, NSUInteger *col ) {
		return (GLFloatComplex) (row[0]==col[0] ? k*k : 0);
	}];
    GLLinearTransform *diffOp = [invN2_trans multiply: [diffZZ minus: k2]];
	
    NSArray *system = [diffOp eigensystemWithOrder: NSOrderedAscending];
	
	GLFunction *lambda = [system[0] makeRealIfPossible];
    GLLinearTransform *S = [system[1] makeRealIfPossible];
	
	if (self.maximumModes) {
		S = [S reducedFromDimensions: [NSString stringWithFormat: @"0:%lu", self.maximumModes-1] toDimension: @":"];
		lambda = [lambda variableFromIndexRangeString:[NSString stringWithFormat: @"0:%lu", self.maximumModes-1]];
	}
	   
	self.eigendepths = [lambda scalarDivide: 1.0];
    self.eigenfrequencies = [[[[self.eigendepths abs] times: @(g*k*k)] plus: @(f0*f0)] sqrt];
    self.S = [S normalizeWithFunction: [[self.N2 minus: @(f0*f0)] times: @(1/g)]];
	self.Sprime = [diffZ multiply: S];
	
    return @[self.eigendepths, self.S, self.Sprime];
}

- (NSArray *) internalWaveModesFromDensityProfile: (GLFunction *) rho withFullDimensions: (NSArray *) dimensions forLatitude: (GLFloat) latitude
{
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
	self.k = k;
	self.l = l;
	
	GLFloat f0 = 2*(7.2921e-5)*sin(latitude*M_PI/180);
    GLFloat g = 9.81;
	GLScalar *rho0 = [rho mean];
	
	// First construct N^2
    GLLinearTransform *diffZ1D = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 1 leftBC: kGLNeumannBoundaryCondition rightBC:kGLNeumannBoundaryCondition bandwidth:1 fromDimension:zDim forEquation:equation];
    self.N2 = [diffZ1D transform: [[rho dividedBy: rho0] times: @(-g)]];
	
    GLFunction *invN2 = [[self.N2 minus: @(f0*f0)] scalarDivide: -g];
    
	// Now construct A = k*k*eye(N) - Diff2;
    GLLinearTransform *diffZZ1D = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 2 leftBC: kGLDirichletBoundaryCondition rightBC:kGLDirichletBoundaryCondition bandwidth:1 fromDimension:zDim forEquation:equation];
	GLLinearTransform *diffZZ = [diffZZ1D expandedWithFromDimensions: transformedDimensions toDimensions: transformedDimensions];

    GLLinearTransform *invN2_trans = [[GLLinearTransform linearTransformFromFunction: invN2] expandedWithFromDimensions: transformedDimensions toDimensions: transformedDimensions];
    GLLinearTransform *diffOp = [invN2_trans multiply: [diffZZ minus: [GLLinearTransform linearTransformFromFunction:K2]]];
	
    NSArray *system = [diffOp eigensystemWithOrder: NSOrderedAscending];
    
	GLFunction *lambda = [system[0] makeRealIfPossible];
	GLLinearTransform *S = [system[1] makeRealIfPossible];
	   
	if (self.maximumModes) {
//		S = [S reducedFromDimensions: [NSString stringWithFormat: @"0:%lu,:,:", self.maximumModes-1] toDimension: @":,:,:"];
//		lambda = [lambda variableFromIndexRangeString:[NSString stringWithFormat: @"0:%lu,:,:", self.maximumModes-1]];
        S = [S reducedFromDimensions: [NSString stringWithFormat: @":,:,0:%lu", self.maximumModes-1] toDimension: @":,:,:"];
		lambda = [lambda variableFromIndexRangeString:[NSString stringWithFormat: @":,:,0:%lu", self.maximumModes-1]];
	}
    
    GLLinearTransform *diffZ = [diffZ1D expandedWithFromDimensions: S.toDimensions toDimensions:S.toDimensions];
	   
    self.eigendepths = [lambda scalarDivide: 1.0];
    
    
    
    //self.eigenfrequencies = [[[[self.eigendepths abs] multiply: [K2 times: @(g)]] plus: @(f0*f0)] sqrt];
    self.S = [S normalizeWithFunction: [[self.N2 minus: @(f0*f0)] times: @(1/g)]];
	self.Sprime = [diffZ multiply: S];
	
    [self.S dumpToConsole];
    
    return @[self.eigendepths, self.S, self.Sprime];
}


- (NSArray *) internalModesGIPFromDensityProfile: (GLFunction *) rho wavenumber: (GLFloat) k latitude: (GLFloat) latitude
{
    if (rho.dimensions.count != 1) {
        [NSException raise:@"InvalidDimensions" format:@"Only one dimension allowed, at this point"];
    }
	
	GLFloat f0 = 2*(7.2921e-5)*sin(latitude*M_PI/180);
	GLFloat g = 9.81;
    GLScalar *rho0 = [rho mean];
	
    GLEquation *equation = rho.equation;
    GLDimension *zDim = rho.dimensions[0];
    
	// First construct N^2
    GLLinearTransform *diffZ = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 1 leftBC: kGLNeumannBoundaryCondition rightBC:kGLNeumannBoundaryCondition bandwidth:1 fromDimension:zDim forEquation:equation];
    self.N2 = [diffZ transform: [[rho dividedBy: rho0] times: @(-g)]];
	
	// Now construct A = k*k*eye(N) - Diff2;
	GLLinearTransform *k2 = [GLLinearTransform transformOfType: kGLRealDataFormat withFromDimensions: @[zDim] toDimensions: @[zDim] inFormat: @[@(kGLDiagonalMatrixFormat)] forEquation:equation matrix:^( NSUInteger *row, NSUInteger *col ) {
		return (GLFloatComplex) (row[0]==col[0] ? k*k : 0);
	}];
    GLLinearTransform *diffZZ = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 2 leftBC: kGLDirichletBoundaryCondition rightBC:kGLDirichletBoundaryCondition bandwidth:1 fromDimension:zDim forEquation:equation];
	GLLinearTransform *A = [k2 minus: diffZZ];
    
	// Now construct B = k*k*diag(N2) - f0*f0*Diff2;
	GLLinearTransform *B = [[[GLLinearTransform linearTransformFromFunction: self.N2] times: @(k*k)] minus: [diffZZ times: @(f0*f0)]];
	
    NSArray *system = [B generalizedEigensystemWith: A];
	
	GLFunction *lambda = [system[0] makeRealIfPossible];
    GLLinearTransform *S = [system[1] makeRealIfPossible];
	
	if (self.maximumModes) {
		S = [S reducedFromDimensions: [NSString stringWithFormat: @"0:%lu", self.maximumModes-1] toDimension: @":"];
		lambda = [lambda variableFromIndexRangeString:[NSString stringWithFormat: @"0:%lu", self.maximumModes-1]];
	}
	
    S = [S normalizeWithFunction: [[self.N2 minus: @(f0*f0)] times: rho0]];
	GLLinearTransform *Sprime = [diffZ multiply: S];
	
	GLFunction *omega = [[lambda abs] sqrt];
	
    return @[omega, S, Sprime];
}

- (NSArray *) internalWaveModesGIPFromDensityProfile: (GLFunction *) rho withFullDimensions: (NSArray *) dimensions forLatitude: (GLFloat) latitude
{
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
	self.k = k;
	self.l = l;
	
	GLFloat f0 = 2*(7.2921e-5)*sin(latitude*M_PI/180);
    GLFloat g = 9.81;
	GLScalar *rho0 = [rho mean];
	
	// First construct N^2
    GLLinearTransform *diffZ1D = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 1 leftBC: kGLNeumannBoundaryCondition rightBC:kGLNeumannBoundaryCondition bandwidth:1 fromDimension:zDim forEquation:equation];
    self.N2 = [diffZ1D transform: [[rho dividedBy: rho0] times: @(-g)]];
	
	// Now construct A = k*k*eye(N) - Diff2;
    GLLinearTransform *diffZZ1D = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 2 leftBC: kGLDirichletBoundaryCondition rightBC:kGLDirichletBoundaryCondition bandwidth:1 fromDimension:zDim forEquation:equation];
	GLLinearTransform *diffZZ = [diffZZ1D expandedWithFromDimensions: transformedDimensions toDimensions: transformedDimensions];
	GLLinearTransform *A = [[GLLinearTransform linearTransformFromFunction:K2] minus: diffZZ];
	
	// Now construct B = k*k*diag(N2) - f0*f0*Diff2;
	GLLinearTransform *B = [[GLLinearTransform linearTransformFromFunction: [K2 multiply: self.N2]] minus: [diffZZ times: @(f0*f0)]];
	
	NSArray *system = [B generalizedEigensystemWith: A];
	
	GLFunction *lambda = [system[0] makeRealIfPossible];
	GLLinearTransform *S = [system[1] makeRealIfPossible];
	
	if (self.maximumModes) {
		S = [S reducedFromDimensions: [NSString stringWithFormat: @"0:%lu,:,:", self.maximumModes-1] toDimension: @":,:,:"];
		lambda = [lambda variableFromIndexRangeString:[NSString stringWithFormat: @"0:%lu,:,:", self.maximumModes-1]];
//        S = [S reducedFromDimensions: [NSString stringWithFormat: @":,:,0:%lu", self.maximumModes-1] toDimension: @":,:,:"];
//		lambda = [lambda variableFromIndexRangeString:[NSString stringWithFormat: @":,:,0:%lu", self.maximumModes-1]];
	}
	
	lambda = [lambda setValue: 0.0 atIndices: @":,0,0"];
	//GLFloat deltaK = lDim.nPoints * kDim.nPoints;
    S = [S normalizeWithFunction: [[self.N2 minus: @(f0*f0)] times: rho0]];
	
	NSUInteger index = 0;
	//	NSUInteger totalVectors = S.matrixDescription.nPoints / S.matrixDescription.strides[index].nPoints;
	//	NSUInteger vectorStride = S.matrixDescription.strides[index].columnStride;
	NSUInteger vectorLength = S.matrixDescription.strides[index].nRows;
	NSUInteger vectorElementStride = S.matrixDescription.strides[index].rowStride;
	//	NSUInteger complexStride = S.matrixDescription.strides[index].complexStride;
	
	for (NSUInteger i=0; i<vectorLength; i++) {
		S.pointerValue[i*vectorElementStride] = 0;
	}
    
    GLLinearTransform *diffZ = [diffZ1D expandedWithFromDimensions: S.toDimensions toDimensions:S.toDimensions];
    GLLinearTransform *Sprime = [diffZ multiply: S];
	
	GLFunction *omega = [[lambda abs] sqrt];
	
    return @[omega, S, Sprime];
}

@end
