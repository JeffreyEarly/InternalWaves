//
//  GLInternalWaveInitialization.h
//  InternalWaves
//
//  Created by Jeffrey J. Early on 1/14/14.
//  Copyright (c) 2014 Jeffrey J. Early. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>

@interface GLInternalWaveInitialization : NSObject

/** Compute the internal wave modes a set of wavenumbers.
 @param rho Density profile given as a function of z only.
 @param dimensions A set of vertical and horizontal dimensions. The same vertical dimensions (z) must be included, and at least one horizontal dimension.
 @param latitude Latitude at which the modes will be used (for the Coriolis frequency).
 @param equation The GLEquation object that sould be used for all calculations.
 @returns A GLInternalWaveInitialization object with N2, eigenfrequencies, S, and Sprime and all the variable phases populated.
 */
- (GLInternalWaveInitialization *) initWithDensityProfile: (GLFunction *) rho fullDimensions: (NSArray *) dimensions latitude: (GLFloat) latitude equation: (GLEquation *) equation;

/// Optional maximum number of modes to be included in the transformation matrices. Setting this to 0 will use all modes.
@property NSUInteger maximumModes;

/** Initializes all variables with the Garrett-Munk spectrum.
 @param energyLevel A multiplicative factor, use 1.0 for default settings.
 */
- (void) createGarrettMunkSpectrumWithEnergy: (GLFloat) energyLevel;

/// The equation used for all computations.
@property(strong) GLEquation *equation;

/// The spatial dimensions, e.g., (x, y, z), although they will be in the order given during initialization.
@property(strong) NSArray *fullDimensions;

/// The associated spectral dimensions, e.g., (k, l, mode).
@property(strong) NSArray *spectralDimensions;

/// The latitude used for the mode and phase computation.
@property GLFloat latitude;

/// The Coriolis frequency at the given latitude.
@property GLFloat f0;

/// Density profile (function of z only).
@property(strong) GLFunction *rho;

/// Stratification profile used in the calculations.
@property(strong) GLFunction *N2;

@property(strong) GLFunction *eigenfrequencies;
@property(strong) GLLinearTransform *S;
@property(strong) GLLinearTransform *Sprime;

@property(strong) GLFunction *zeta_plus;
@property(strong) GLFunction *zeta_minus;
@property(strong) GLFunction *rho_plus;
@property(strong) GLFunction *rho_minus;
@property(strong) GLFunction *u_plus;
@property(strong) GLFunction *u_minus;
@property(strong) GLFunction *v_plus;
@property(strong) GLFunction *v_minus;
@property(strong) GLFunction *w_plus;
@property(strong) GLFunction *w_minus;

@end
