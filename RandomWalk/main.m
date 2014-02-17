//
//  main.m
//  RandomWalk
//
//  Created by Jeffrey Early on 2/16/14.
//  Copyright (c) 2014 Jeffrey J. Early. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>

int main(int argc, const char * argv[])
{

    @autoreleasepool {
        NSUInteger Nx = 10;
        NSUInteger Ny = 10;
        GLFloat Dx = 10;
        GLFloat Dy = 10;
        
        /************************************************************************************************/
        /*		Define the problem dimensions															*/
        /************************************************************************************************/
        GLEquation *equation = [[GLEquation alloc] init];
        GLDimension *xFloatDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: Nx domainMin: -(Dx-1)/2 length: Dx-1];
        xFloatDim.name = @"x";
        GLDimension *yFloatDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: Ny domainMin: -(Dy-1)/2 length: Dy-1];
        yFloatDim.name = @"y";
        
        NSArray *floatDimensions = @[xFloatDim, yFloatDim];
		GLFunction *xFloat = [GLFunction functionOfRealTypeFromDimension: xFloatDim withDimensions: floatDimensions forEquation: equation];
		GLFunction *yFloat = [GLFunction functionOfRealTypeFromDimension: yFloatDim withDimensions: floatDimensions forEquation: equation];
        
		GLFunction *xPosition = [GLFunction functionFromFunction: xFloat];
		GLFunction *yPosition = [GLFunction functionFromFunction: yFloat];
        
        GLFloat timeStep = 10;
        GLFloat kappa = 1; // m^2/s
        GLFloat norm = sqrt(timeStep*2*kappa);
        norm = sqrt(4)*norm/timeStep; // the integrator multiplies by deltaT, so we account for that here.
        // There's something else I'm not getting. Some reason why we need that sqrt(4) there. It's differen for RK5.
        // RK4: dt/3 f(0) + dt/6 f(1) + dt/6 *f(4) + dt/3*f(3)
        // Mean of 1/3 and 1/6? 1/4. It's like the geometric mean to get the same norm, yup. Hence, sqrt 4.
        // RK5 will depend on those C coefficients.
        GLAdaptiveRungeKuttaOperation *integrator = [GLAdaptiveRungeKuttaOperation rungeKutta4AdvanceY: @[xPosition, yPosition] stepSize: timeStep fFromTY:^(GLScalar *time, NSArray *yNew) {
            GLFunction *xStep = [GLFunction functionWithNormallyDistributedValueWithDimensions: floatDimensions forEquation: equation];
            GLFunction *yStep = [GLFunction functionWithNormallyDistributedValueWithDimensions: floatDimensions forEquation: equation];
            
            xStep = [xStep times: @(norm)];
            yStep = [yStep times: @(norm)];
            
			return @[xStep, yStep];
		}];
        
        GLFloat maxTime = 10000;
        NSArray *newPositions = [integrator stepForwardToTime: maxTime];
        GLFunction *x = newPositions[0];
        GLFunction *y = newPositions[1];
        
        GLScalar *meanSquareSeparation0 = [[[yPosition times: yPosition] plus: [yPosition times: yPosition]] mean];
        GLScalar *meanSquareSeparation = [[[x times: x] plus: [y times: y]] mean];
        
        GLFloat a = *(meanSquareSeparation0.pointerValue);
        GLFloat b = *(meanSquareSeparation.pointerValue);
        [meanSquareSeparation0 dumpToConsole];
        [meanSquareSeparation dumpToConsole];
        
        GLFloat kappaDeduced = (0.25)*(b-a)/maxTime;
        NSLog(@"kappa: %f, actual kappa: %f", kappa, kappaDeduced);
    }
    return 0;
}

