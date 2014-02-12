//
//  GLInternalWaveInitialization.m
//  InternalWaves
//
//  Created by Jeffrey J. Early on 1/14/14.
//  Copyright (c) 2014 Jeffrey J. Early. All rights reserved.
//

#import "GLInternalWaveInitialization.h"
#import "GLInternalModes.h"

#define g 9.81

@interface GLInternalWaveInitialization ()

@property(strong) GLDimension *kDim;
@property(strong) GLDimension *lDim;
@property(strong) GLDimension *modeDim;
@property(strong) NSMutableArray *horizontalDimensions;

- (void) generateWavePhasesFromPositive: (GLFunction *) G_plus negative: (GLFunction *) G_minus;

@end

static NSString *GLInternalWaveMaximumModesKey = @"GLInternalWaveMaximumModesKey";
static NSString *GLInternalWaveEquationKey = @"GLInternalWaveEquationKey";
static NSString *GLInternalWaveFullDimensionsKey = @"GLInternalWaveFullDimensionsKey";
static NSString *GLInternalWaveSpectralDimensionsKey = @"GLInternalWaveSpectralDimensionsKey";
static NSString *GLInternalWaveLatitudeKey = @"GLInternalWaveLatitudeKey";
static NSString *GLInternalWaveF0Key = @"GLInternalWaveF0Key";
static NSString *GLInternalWaveRhoKey = @"GLInternalWaveRhoKey";
static NSString *GLInternalWaveN2Key = @"GLInternalWaveN2Key";
static NSString *GLInternalWaveEigenfrequenciesKey = @"GLInternalWaveEigenfrequenciesKey";
static NSString *GLInternalWaveSKey = @"GLInternalWaveSKey";
static NSString *GLInternalWaveSprimeKey = @"GLInternalWaveSprimeKey";
static NSString *GLInternalWaveZetaPlusKey = @"GLInternalWaveZetaPlusKey";
static NSString *GLInternalWaveZetaMinusKey = @"GLInternalWaveZetaMinusKey";
static NSString *GLInternalWaveRhoPlusKey = @"GLInternalWaveRhoPlusKey";
static NSString *GLInternalWaveRhoMinusKey = @"GLInternalWaveRhoMinusKey";
static NSString *GLInternalWaveUPlusKey = @"GLInternalWaveUPlusKey";
static NSString *GLInternalWaveUMinusKey = @"GLInternalWaveUMinusKey";
static NSString *GLInternalWaveVPlusKey = @"GLInternalWaveVPlusKey";
static NSString *GLInternalWaveVMinusKey = @"GLInternalWaveVMinusKey";
static NSString *GLInternalWaveWPlusKey = @"GLInternalWaveWPlusKey";
static NSString *GLInternalWaveWMinusKey = @"GLInternalWaveWMinusKey";

@implementation GLInternalWaveInitialization

- (void)encodeWithCoder:(NSCoder *)coder
{
    [coder encodeObject: @(self.maximumModes) forKey:GLInternalWaveMaximumModesKey];
    [coder encodeObject: self.equation forKey:GLInternalWaveEquationKey];
    [coder encodeObject: self.fullDimensions forKey:GLInternalWaveFullDimensionsKey];
    [coder encodeObject: self.spectralDimensions forKey:GLInternalWaveSpectralDimensionsKey];
    [coder encodeObject: @(self.latitude) forKey:GLInternalWaveLatitudeKey];
    [coder encodeObject: @(self.f0) forKey:GLInternalWaveF0Key];
    [coder encodeObject: self.rho forKey:GLInternalWaveRhoKey];
    [coder encodeObject: self.N2 forKey:GLInternalWaveN2Key];
    [coder encodeObject: self.eigenfrequencies forKey:GLInternalWaveEigenfrequenciesKey];
    [coder encodeObject: self.S forKey:GLInternalWaveSKey];
    [coder encodeObject: self.Sprime forKey:GLInternalWaveSprimeKey];
    [coder encodeObject: self.zeta_plus forKey:GLInternalWaveZetaPlusKey];
    [coder encodeObject: self.zeta_minus forKey:GLInternalWaveZetaMinusKey];
    [coder encodeObject: self.rho_plus forKey:GLInternalWaveRhoPlusKey];
    [coder encodeObject: self.rho_minus forKey:GLInternalWaveRhoMinusKey];
    [coder encodeObject: self.u_plus forKey:GLInternalWaveUPlusKey];
    [coder encodeObject: self.u_minus forKey:GLInternalWaveUMinusKey];
    [coder encodeObject: self.v_plus forKey:GLInternalWaveVPlusKey];
    [coder encodeObject: self.v_minus forKey:GLInternalWaveVMinusKey];
    [coder encodeObject: self.w_plus forKey:GLInternalWaveWPlusKey];
    [coder encodeObject: self.w_minus forKey:GLInternalWaveWMinusKey];
}

- (id)initWithCoder:(NSCoder *)decoder
{
    if ((self=[super init])) {
        _maximumModes = [[decoder decodeObjectForKey: GLInternalWaveMaximumModesKey] unsignedIntegerValue];
        _equation = [decoder decodeObjectForKey: GLInternalWaveEquationKey];
        _fullDimensions = [decoder decodeObjectForKey: GLInternalWaveFullDimensionsKey];
        _spectralDimensions = [decoder decodeObjectForKey: GLInternalWaveSpectralDimensionsKey];
        _latitude = [[decoder decodeObjectForKey: GLInternalWaveLatitudeKey] doubleValue];
        _f0 = [[decoder decodeObjectForKey: GLInternalWaveF0Key] doubleValue];
        _rho = [decoder decodeObjectForKey: GLInternalWaveRhoKey];
        _N2 = [decoder decodeObjectForKey: GLInternalWaveN2Key];
        _eigenfrequencies = [decoder decodeObjectForKey: GLInternalWaveEigenfrequenciesKey];
        _S = [decoder decodeObjectForKey: GLInternalWaveSKey];
        _Sprime = [decoder decodeObjectForKey: GLInternalWaveSprimeKey];
        _zeta_plus = [decoder decodeObjectForKey: GLInternalWaveZetaPlusKey];
        _zeta_minus = [decoder decodeObjectForKey: GLInternalWaveZetaMinusKey];
        _rho_plus = [decoder decodeObjectForKey: GLInternalWaveRhoPlusKey];
        _rho_minus = [decoder decodeObjectForKey: GLInternalWaveRhoMinusKey];
        _u_plus = [decoder decodeObjectForKey: GLInternalWaveUPlusKey];
        _u_minus = [decoder decodeObjectForKey: GLInternalWaveUMinusKey];
        _v_plus = [decoder decodeObjectForKey: GLInternalWaveVPlusKey];
        _v_minus = [decoder decodeObjectForKey: GLInternalWaveVMinusKey];
        _w_plus = [decoder decodeObjectForKey: GLInternalWaveWPlusKey];
        _w_minus = [decoder decodeObjectForKey: GLInternalWaveWMinusKey];
    }
    return self;
}

- (GLInternalWaveInitialization *) initWithDensityProfile: (GLFunction *) rho fullDimensions: (NSArray *) dimensions latitude: (GLFloat) latitude equation: (GLEquation *) equation
{
	NSUInteger numVerticalDims = 0;
	NSUInteger numHorizontalDims = 0;
	for (GLDimension *dim in dimensions) {
		if ([dim.name isEqualToString: @"x"] || [dim.name isEqualToString: @"y"]) {
			numHorizontalDims++;
		} else if ([dim.name isEqualToString: @"z"]) {
			numVerticalDims++;
		} else {
			[NSException raise: @"BadDimensions" format:@"Dimensions must be name x, y, or z"];
		}
	}
	if (numVerticalDims != 1 || numHorizontalDims == 0 || numHorizontalDims > 2) {
		[NSException raise: @"BadDimensions" format:@"There must be one vertical dimension given, and either 1 or 2 horizontal dimensions."];
	}
	
    if ((self=[super init])) {
        self.fullDimensions=dimensions;
        self.equation=equation;
		self.horizontalDimensions = [NSMutableArray array];
		self.rho = rho;
		self.latitude = latitude;
		self.f0 = 2*(7.2921e-5)*sin(latitude*M_PI/180);
    }
    return self;
}

- (void) computeInternalModes
{
	NSLog(@"Computing internal modes...");
	
	GLInternalModes *internalModes = [[GLInternalModes alloc] init];
	internalModes.maximumModes = self.maximumModes;
	[internalModes internalWaveModesFromDensityProfile: self.rho withFullDimensions: self.fullDimensions forLatitude: self.latitude];
	
	self.eigenfrequencies = internalModes.eigenfrequencies;
    self.eigendepths = internalModes.eigendepths;
    self.rossbyRadius = internalModes.rossbyRadius;
	self.S = internalModes.S;
	self.Sprime = internalModes.Sprime;
	self.N2 = internalModes.N2;
	self.spectralDimensions = self.eigenfrequencies.dimensions;
	
	for (GLDimension *dim in self.spectralDimensions) {
		if ( [dim.name isEqualToString: @"k"]) {
			self.kDim = dim;
			[self.horizontalDimensions addObject: dim];
		} else if ( [dim.name isEqualToString: @"l"]) {
			self.lDim = dim;
			[self.horizontalDimensions addObject: dim];
		} else {
			self.modeDim = dim;
		}
	}
}

- (GLFloat) createUnitWaveWithSpeed: (GLFloat) U_max verticalMode: (NSUInteger) mode k: (NSUInteger) kUnit l: (NSUInteger) lUnit omegaSign: (GLFloat) sign
{
	[self computeInternalModes];
	
    if (mode == 0) {
        [NSException raise: @"UnsupportedMode" format:@"Barotropic mode is not supported. Set the mode to 1 or greater"];
    }
    // The 1st baroclinic mode is in position 0.
    mode = mode-1;
    

    
    GLFunction *U_mag = [GLFunction functionOfRealTypeWithDimensions: self.spectralDimensions forEquation: self.equation];
    [U_mag zero];
	
	// My notation is horrible here. kDim and lDim are actually reversed!
    NSUInteger zDimNPoints = [U_mag.dimensions[0] nPoints];
    NSUInteger kDimNPoints = [U_mag.dimensions[1] nPoints];
    NSUInteger lDimNPoints = [U_mag.dimensions[2] nPoints];

    // We will be setting G, the energy density. So we have to convert U_max into that value.
//    GLFloat k = [self.kDim valueAtIndex: kUnit];
//    GLFloat l = [self.lDim valueAtIndex: lUnit];
	GLFloat omega = self.eigenfrequencies.pointerValue[(mode*kDimNPoints+kUnit)*lDimNPoints+lUnit];
//
//    GLFloat G = U_max*(k*k+l*l)/(k*omega);

    U_mag = [U_mag setValue: 1.0 atIndices: [NSString stringWithFormat:@"%lu,%lu,%lu", mode, kUnit, lUnit]];
	
    GLFloat *C = U_mag.pointerValue;
    
    // index=i*ny*nz+j*nz+k
    // index=(i*ny+j)*nz+k
    // my notation, (z*kDimNPoints+k)*lDimNPoints+l
    for (NSUInteger z=0; z<zDimNPoints; z++) {
        // Hermitian conjugates
        for (NSUInteger i=1; i<kDimNPoints/2; i++) {
            C[(z*kDimNPoints+(kDimNPoints-i))*lDimNPoints] = C[(z*kDimNPoints+i)*lDimNPoints];
            C[(z*kDimNPoints+(kDimNPoints-i))*lDimNPoints+(lDimNPoints-1)] = C[(z*kDimNPoints+(kDimNPoints-i))*lDimNPoints+(lDimNPoints-1)];
        }
        
        // For the four self-conjugate components, that means that there can be no imaginary component
        // But that their real components should be doubled, to make all else equal.
        C[z*kDimNPoints*lDimNPoints+0] *= 2;
        C[z*kDimNPoints*lDimNPoints+(lDimNPoints-1)] *= 2;
        C[(z*kDimNPoints+(kDimNPoints/2))*lDimNPoints+0] *= 2;
        C[(z*kDimNPoints+(kDimNPoints/2))*lDimNPoints+(lDimNPoints-1)] *= 2;
    }
    
	GLScalar *i = [GLScalar scalarWithValue: I forEquation: self.equation];
    GLFunction *G_plus = [U_mag duplicate];
    GLFunction *G_minus = [[U_mag duplicate] times: i]; // Keep the two sides out of phase initially.
    
    if (sign<0) {
        [G_plus zero];
    } else if (sign > 0) {
        [G_minus zero];
    }
    
    [self generateWavePhasesFromPositive: G_plus negative: G_minus];
    
    // Rather than figure out how to properly normalize, we just cheat.
    GLFloat U0 = [[[self.Sprime transform: [self.u_plus plus: self.u_minus]] abs] maxNow];
    GLFloat V0 = [[[self.Sprime transform: [self.v_plus plus: self.v_minus]] abs] maxNow];
    GLFloat Speed = sqrt(U0*U0+V0*V0);
    
    G_plus = [G_plus times: @(U_max/(2*Speed))];
    G_minus = [G_minus times: @(U_max/(2*Speed))];
    
    [self generateWavePhasesFromPositive: G_plus negative: G_minus];
	
	return omega;
}

- (void) createGarrettMunkSpectrumWithEnergy: (GLFloat) energyLevel
{
	[self computeInternalModes];
	
	// We follow Winters and D'Asaro (1997) for the creation of the Garrett-Munk spectrum
	GLFloat j_star = 3;
	GLFloat E = (2.2e-6)*energyLevel; // [ m^{2}/s ]
	
	// The mode dimension, j, starts at zero, but we want it to start at 1... so we add 1!
	GLFunction *j = [[GLFunction functionOfRealTypeFromDimension: self.modeDim withDimensions: self.spectralDimensions forEquation: self.equation] plus: @(1)];
	// H(j) = 2*j_star/(pi*(j^2 + j_star^2)) [unitless]
	GLFunction *H = [[[j plus: @(j_star)] pow: -5/2] scalarDivide: 2*pow(j_star,3/2)/M_PI];
	
	GLFunction *k = [[GLFunction functionOfRealTypeFromDimension: self.kDim withDimensions: self.spectralDimensions forEquation:self.equation] scalarMultiply: 2*M_PI];
	GLFunction *l = [[GLFunction functionOfRealTypeFromDimension: self.lDim withDimensions: self.spectralDimensions forEquation:self.equation] scalarMultiply: 2*M_PI];
    GLFunction *K2 = [[k multiply: k] plus: [l multiply: l]];
	
	// B(alpha,j) = (2/pi) * 1/R_j * 1/(k^2 + R_j^-2) [m]
	GLFunction *invR_j = [self.rossbyRadius scalarDivide: 1.0];
	GLFunction *B = [[invR_j dividedBy: [K2 plus: [invR_j multiply: invR_j]]] multiply: @(2./M_PI)];
	
	// G^2(alpha, j) = E*H(j)*B(alpha,j)	[J/m]
	GLFunction *G2 = [[H multiply: B] multiply: @(E)];
	
	// G(j,l,k) = G(alpha,j)/sqrt(2*pi*alpha) where alpha = sqrt(k^2 + l^2); [sqrt(kg) m/s]
	// We cut the value in half, because we will want the expectation of the the positive and negative sides to add up to this value.
	GLFunction *G = [[[G2 dividedBy: [[K2 sqrt] times: @(2*M_PI)]] sqrt] times: @(0.5)];
	G = [G setValue: 0 atIndices: @":,0,0"]; // We just divided by zero, so we need to get rid of that component.
	
	// <G_+> = 1/2 <G> = <G_->
	GLFunction *G_plus = [GLFunction functionWithNormallyDistributedValueWithDimensions: self.spectralDimensions forEquation: self.equation];
	GLFunction *G_minus = [GLFunction functionWithNormallyDistributedValueWithDimensions: self.spectralDimensions forEquation: self.equation];
	G_plus = [G_plus multiply: G];
	G_minus = [G_minus multiply: G];
	
    [self generateWavePhasesFromPositive: G_plus negative: G_minus];
}

- (void) generateWavePhasesFromPositive: (GLFunction *) G_plus negative: (GLFunction *) G_minus
{
    GLFunction *k = [[GLFunction functionOfRealTypeFromDimension: self.kDim withDimensions: self.spectralDimensions forEquation:self.equation] scalarMultiply: 2*M_PI];
	GLFunction *l = [[GLFunction functionOfRealTypeFromDimension: self.lDim withDimensions: self.spectralDimensions forEquation:self.equation] scalarMultiply: 2*M_PI];
    GLFunction *K_H = [[[k multiply: k] plus: [l multiply: l]] sqrt];
	K_H = [K_H setValue: 1 atIndices: @":,0,0"]; // prevent divide by zero.
    GLFunction *sqrtH = [self.eigendepths sqrt];
    
	GLScalar *rho0 = [self.rho mean];
    
	self.zeta_plus = [G_plus multiply: [[K_H multiply: sqrtH] dividedBy: self.eigenfrequencies]];
    self.zeta_minus = [G_minus multiply: [[K_H multiply: sqrtH] dividedBy: self.eigenfrequencies]];
	
	self.rho_plus = [[self.zeta_plus times: rho0] times: @(1/g)];
    self.rho_minus = [[self.zeta_minus times: rho0] times: @(1/g)];
	
	self.w_plus = [[[G_plus multiply: [K_H multiply: sqrtH]] swapComplex] makeHermitian];
    self.w_minus = [[[[G_minus multiply: [K_H multiply: sqrtH]] swapComplex] negate] makeHermitian];
    
    GLFunction *denominator = [[self.eigenfrequencies multiply: K_H] multiply: sqrtH];
	self.u_plus = [[G_plus multiply: [[[[k multiply: self.eigenfrequencies] minus: [[l times: @(self.f0)] swapComplex]] dividedBy: denominator] negate]] makeHermitian];
    self.u_minus = [[G_minus multiply: [[[k multiply: self.eigenfrequencies] plus: [[l times: @(self.f0)] swapComplex]] dividedBy: denominator]] makeHermitian];
	
	self.v_plus = [[G_plus multiply: [[[[l multiply: self.eigenfrequencies] plus: [[k times: @(self.f0)] swapComplex]] dividedBy: denominator] negate]] makeHermitian];
    self.v_minus = [[G_minus multiply: [[[l multiply: self.eigenfrequencies] minus: [[k times: @(self.f0)] swapComplex]] dividedBy: denominator]] makeHermitian];
}


@end
