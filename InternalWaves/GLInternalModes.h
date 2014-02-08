//
//  GLInternalModes.h
//  InternalWaves
//
//  Created by Jeffrey J. Early on 1/14/14.
//  Copyright (c) 2014 Jeffrey J. Early. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>

@interface GLInternalModes : NSObject

/** Compute the internal geostrophic (omega=0) modes.
 @param rho Density profile given as a function of z only.
 @returns A GLInternalModes object with N2, eigenfrequencies, S, and Sprime populated.
 */
- (NSArray *) internalModesFromDensityProfile: (GLFunction *) rho;

/** Compute the internal wave modes for a given wavenumber.
 @param rho Density profile given as a function of z only.
 @param k Wavenumber given in cycles/meter.
 @param latitude Latitude at which the modes will be used (for the Coriolis frequency).
 @returns A GLInternalModes object with N2, eigenfrequencies, S, and Sprime populated.
 */
- (NSArray *) internalModesFromDensityProfile: (GLFunction *) rho wavenumber: (GLFloat) k latitude: (GLFloat) latitude;

/** Compute the internal wave modes a set of wavenumbers.
 @param rho Density profile given as a function of z only.
 @param dimensions A set of vertical and horizontal dimensions. The same vertical dimensions (z) must be included, and at least one horizontal dimension.
 @param latitude Latitude at which the modes will be used (for the Coriolis frequency).
 @returns A GLInternalModes object with N2, eigenfrequencies, S, and Sprime populated.
 */
- (NSArray *) internalWaveModesFromDensityProfile: (GLFunction *) rho withFullDimensions: (NSArray *) dimensions forLatitude: (GLFloat) latitude;
- (NSArray *) internalWaveModesGIPFromDensityProfile: (GLFunction *) rho withFullDimensions: (NSArray *) dimensions forLatitude: (GLFloat) latitude;

/// Set the maximum number of modes to be used in the transformation matrix.
@property NSUInteger maximumModes;

/// Stratification profile as compute for the mode calculation. N^2(z) = -g/mean(rho) * d/dz(rho)
@property(strong) GLFunction *N2;

/// Wavenumber function associated with x. This may be nil.
@property(strong) GLFunction *k;

/// Wavenumber function associated with y. This may be nil.
@property(strong) GLFunction *l;

// Not used yet.
@property(strong) GLFunction *eigenfrequencies;
@property(strong) GLFunction *eigendepths;
@property(strong) GLLinearTransform *S;
@property(strong) GLLinearTransform *Sprime;


@end
