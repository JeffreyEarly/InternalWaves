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
@property NSUInteger maximumModes;
@end

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
    GLFunction *N2 = [diffZ transform: [[rho dividedBy: rho0] times: @(1)]];
	
	
    GLFunction *invN2 = [N2 scalarDivide: 1.0];
		
    GLLinearTransform *diffZZ = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 2 leftBC: kGLDirichletBoundaryCondition rightBC:kGLDirichletBoundaryCondition bandwidth:1 fromDimension:zDim forEquation:equation];
    GLLinearTransform *invN2_trans = [GLLinearTransform linearTransformFromFunction: invN2];
    GLLinearTransform *diffOp = [invN2_trans multiply: diffZZ];
	
    NSArray *system = [diffOp eigensystem];
	
	GLFunction *lambda = [system[0] makeRealIfPossible];
    GLLinearTransform *S = [system[1] makeRealIfPossible];
	
	if (self.maximumModes) {
		S = [S reducedFromDimensions: [NSString stringWithFormat: @"0:%lu", self.maximumModes-1] toDimension: @":"];
	}
	
    S = [S normalizeWithFunction: [N2 times: rho0]];
	GLLinearTransform *Sprime = [diffZ multiply: S];
	
    return @[lambda, S, Sprime];
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
    GLFunction *N2 = [diffZ transform: [[rho dividedBy: rho0] times: @(-g)]];
	
	// Now construct A = k*k*eye(N) - Diff2;
	GLLinearTransform *k2 = [GLLinearTransform transformOfType: kGLRealDataFormat withFromDimensions: @[zDim] toDimensions: @[zDim] inFormat: @[@(kGLDiagonalMatrixFormat)] forEquation:equation matrix:^( NSUInteger *row, NSUInteger *col ) {
		return (GLFloatComplex) (row[0]==col[0] ? k*k : 0);
	}];
    GLLinearTransform *diffZZ = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 2 leftBC: kGLDirichletBoundaryCondition rightBC:kGLDirichletBoundaryCondition bandwidth:1 fromDimension:zDim forEquation:equation];
	GLLinearTransform *A = [k2 minus: diffZZ];
    
	// Now construct B = k*k*diag(N2) - f0*f0*Diff2;
	GLLinearTransform *B = [[[GLLinearTransform linearTransformFromFunction: N2] times: @(k*k)] minus: [diffZZ times: @(f0*f0)]];
		
    NSArray *system = [B generalizedEigensystemWith: A];
	
	GLFunction *lambda = [system[0] makeRealIfPossible];
    GLLinearTransform *S = [system[1] makeRealIfPossible];
	
	if (self.maximumModes) {
		S = [S reducedFromDimensions: [NSString stringWithFormat: @"0:%lu", self.maximumModes-1] toDimension: @":"];
	}
	
    S = [S normalizeWithFunction: [[N2 minus: @(f0*f0)] times: rho0]];
	GLLinearTransform *Sprime = [diffZ multiply: S];
	
    return @[lambda, S, Sprime];
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
    GLFloat g = 9.81;
	GLScalar *rho0 = [rho mean];
	
	// First construct N^2
    GLLinearTransform *diffZ1D = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 1 leftBC: kGLNeumannBoundaryCondition rightBC:kGLNeumannBoundaryCondition bandwidth:1 fromDimension:zDim forEquation:equation];
    GLFunction *N2 = [diffZ1D transform: [[rho dividedBy: rho0] times: @(-g)]];

	// Now construct A = k*k*eye(N) - Diff2;
    GLLinearTransform *diffZZ1D = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 2 leftBC: kGLDirichletBoundaryCondition rightBC:kGLDirichletBoundaryCondition bandwidth:1 fromDimension:zDim forEquation:equation];
	GLLinearTransform *diffZZ = [diffZZ1D expandedWithFromDimensions: transformedDimensions toDimensions: transformedDimensions];
	GLLinearTransform *A = [[GLLinearTransform linearTransformFromFunction:K2] minus: diffZZ];

	// Now construct B = k*k*diag(N2) - f0*f0*Diff2;
	GLLinearTransform *B = [[GLLinearTransform linearTransformFromFunction: [K2 multiply: N2]] minus: [diffZZ times: @(f0*f0)]];
	
	NSArray *system = [B generalizedEigensystemWith: A];
	
	GLFunction *lambda = [system[0] makeRealIfPossible];
	GLLinearTransform *S = [system[1] makeRealIfPossible];
	
	if (self.maximumModes) {
		S = [S reducedFromDimensions: [NSString stringWithFormat: @"0:%lu,:,:", self.maximumModes-1] toDimension: @":,:,:"];
	}
	
	lambda = [system[0] setValue: 0.0 atIndices: @":,0,0"];
    S = [S normalizeWithFunction: [[N2 minus: @(f0*f0)] times: rho0]];
	
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
	
    return @[lambda, S, Sprime];
}

@end


@interface GLInternalWaveInitialization : NSObject
@property(retain) NSArray *spectralDimensions;
@property(retain) GLEquation *equation;
@end

@implementation GLInternalWaveInitialization
- (id) initWithSpectralDimensions: (NSArray *) dimensions equation: (GLEquation *) equation
{
    if ((self=[super init])) {
        self.spectralDimensions=dimensions;
        self.equation=equation;
    }
    return self;
}

- (GLFunction *) spectrumWithSpeed: (GLFloat) U_max verticalMode: (NSUInteger) mode k: (NSUInteger) kUnit l: (NSUInteger) lUnit
{
    GLFunction *U_mag = [GLFunction functionOfRealTypeWithDimensions: self.spectralDimensions forEquation: self.equation];
    [U_mag zero];
    U_mag = [U_mag setValue: U_max/(4*sqrt(2.0)) atIndices: [NSString stringWithFormat:@"%lu,%lu,%lu", mode, lUnit, kUnit]];
    
    NSUInteger zDimNPoints = [U_mag.dimensions[0] nPoints];
    NSUInteger kDimNPoints = [U_mag.dimensions[1] nPoints];
    NSUInteger lDimNPoints = [U_mag.dimensions[2] nPoints];
        
    GLSplitComplex C = splitComplexFromData( U_mag.data );
    
    // index=i*ny*nz+j*nz+k
    // index=(i*ny+j)*nz+k
    // my notation, (z*kDimNPoints+k)*lDimNPoints+l
    for (NSUInteger z=0; z<zDimNPoints; z++) {
        // Hermitian conjugates
        for (NSUInteger i=1; i<kDimNPoints/2; i++) {
            C.realp[(z*kDimNPoints+(kDimNPoints-i))*lDimNPoints] = C.realp[(z*kDimNPoints+i)*lDimNPoints];
            C.imagp[(z*kDimNPoints+(kDimNPoints-i))*lDimNPoints] = -C.imagp[(z*kDimNPoints+i)*lDimNPoints];
            
            C.realp[(z*kDimNPoints+(kDimNPoints-i))*lDimNPoints+(lDimNPoints-1)] = C.realp[(z*kDimNPoints+(kDimNPoints-i))*lDimNPoints+(lDimNPoints-1)];
            C.imagp[(z*kDimNPoints+(kDimNPoints-i))*lDimNPoints+(lDimNPoints-1)] = -C.imagp[(z*kDimNPoints+(kDimNPoints-i))*lDimNPoints+(lDimNPoints-1)];
        }
        
        // For the four self-conjugate components, that means that there can be no imaginary component
        C.imagp[z*kDimNPoints*lDimNPoints+0] = 0;
        C.imagp[z*kDimNPoints*lDimNPoints+(lDimNPoints-1)] = 0;
        C.imagp[(z*kDimNPoints+(kDimNPoints/2))*lDimNPoints+0] = 0;
        C.imagp[(z*kDimNPoints+(kDimNPoints/2))*lDimNPoints+(lDimNPoints-1)] = 0;
        
        // But that their real components should be doubled, to make all else equal.
        C.realp[z*kDimNPoints*lDimNPoints+0] *= 2;
        C.realp[z*kDimNPoints*lDimNPoints+(lDimNPoints-1)] *= 2;
        C.realp[(z*kDimNPoints+(kDimNPoints/2))*lDimNPoints+0] *= 2;
        C.realp[(z*kDimNPoints+(kDimNPoints/2))*lDimNPoints+(lDimNPoints-1)] *= 2;
    }
    
    return U_mag;
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
		NSUInteger Nz = 100;
		
		GLFloat f0 = 2*(7.2921e-5)*sin(latitude*M_PI/180);
		GLFloat rho0 = 1025;
		GLFloat g = 9.81;
		
        // This is good for unit testing.
        BOOL shouldUnitTest = NO;
        NSUInteger modeUnit = 0;
		NSUInteger kUnit = 0;
		NSUInteger lUnit = 0;
		NSInteger omegaSign = 0;
        GLFloat U_max = .10;
        
        /************************************************************************************************/
		/*		Define the problem dimensions															*/
		/************************************************************************************************/
		GLEquation *equation = [[GLEquation alloc] init];
		GLDimension *zDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: Nz domainMin: -depth length: depth];
		zDim.name = @"z";
		GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: Nx domainMin: -width/2 length: width];
		xDim.name = @"x";
		GLDimension *yDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: Nx domainMin: -width/2 length: width];
		yDim.name = @"y";
        GLMutableDimension *tDim = [[GLMutableDimension alloc] initWithPoints: @[@(0.0)]];
		tDim.name = @"time";
        
        /************************************************************************************************/
		/*		Create a density profile and compute the internal wave modes                            */
		/************************************************************************************************/
        GLFunction *z = [GLFunction functionOfRealTypeFromDimension:zDim withDimensions:@[zDim] forEquation:equation];
        GLFunction *rho = [[z times: @(-N2*rho0/g)] plus: @(rho0)];

        GLInternalModes *internalModes = [[GLInternalModes alloc] init];
		internalModes.maximumModes = 10;
		NSArray *system = [internalModes internalWaveModesFromDensityProfile: rho withHorizontalDimensions: @[yDim, xDim] forLatitude:latitude];
        //NSArray *system = [internalModes internalModesFromDensityProfile: rho];
		//NSArray *system = [internalModes internalModesFromDensityProfile: rho wavenumber: 0.01 latitude: 45.0];
        GLFunction *eigenvalues = system[0];
		GLLinearTransform *S = system[1];
		GLLinearTransform *Sprime = system[2];
        
        NSArray *spectralDimensions = eigenvalues.dimensions;
        
		eigenvalues = [[[eigenvalues times: @(1/(f0*f0))] abs] sqrt];
        [eigenvalues dumpToConsole];
		
		//[S dumpToConsole];
		[Sprime dumpToConsole];
        
        /************************************************************************************************/
		/*		Now set the magnitude of the wave at each component.									*/
		/************************************************************************************************/
		
		// Set the magnitude of the components
//        GLInternalWaveInitialization *initialization = [[GLInternalWaveInitialization alloc] initWithSpectralDimensions: spectralDimensions equation:equation];
//		GLFunction *U_mag = [initialization spectrumWithSpeed: U_max verticalMode: modeUnit k: kUnit l: lUnit];
        
        // We create the z variable with dimensions in reverse order so that the fft will act on contiguous chunks of memory.
	}
    return 0;
}

