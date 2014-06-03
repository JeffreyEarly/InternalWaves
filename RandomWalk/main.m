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
        NSUInteger N = 100;

        GLEquation *equation = [[GLEquation alloc] init];
        GLDimension *floatDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: N domainMin: 1 length: N];
        floatDim.name = @"x";
        
        NSArray *floatDimensions = @[floatDim];
		GLFunction *xPosition = [GLFunction functionOfRealTypeWithDimensions: floatDimensions forEquation: equation];
		GLFunction *yPosition = [GLFunction functionOfRealTypeWithDimensions: floatDimensions forEquation: equation];
        
        GLFloat timeStep = 10;
		GLFloat maxTime = 1000;
		GLDimension *tDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: 1+round(maxTime/timeStep) domainMin: 0 length: maxTime];
		tDim.name = @"t";
        GLFloat kappa = 1; // m^2/s
        GLFloat norm = sqrt(timeStep*2*kappa);
        norm = sqrt(36./10.)*norm/timeStep; // the integrator multiplies by deltaT, so we account for that here.
        // RK4: dt/3 f(0) + dt/6 f(1) + dt/6 *f(4) + dt/3*f(3)
        // sqrt of ( (1/3)^2 + (1/6)^ + (1/6)^2 + (1/3)^2 )
        GLAdaptiveRungeKuttaOperation *integrator = [GLAdaptiveRungeKuttaOperation rungeKutta4AdvanceY: @[xPosition, yPosition] stepSize: timeStep fFromTY:^(GLScalar *time, NSArray *yNew) {
            GLFunction *xStep = [GLFunction functionWithNormallyDistributedValueWithDimensions: floatDimensions forEquation: equation];
            GLFunction *yStep = [GLFunction functionWithNormallyDistributedValueWithDimensions: floatDimensions forEquation: equation];
            
            xStep = [xStep times: @(norm)];
            yStep = [yStep times: @(norm)];
            
			return @[xStep, yStep];
		}];
        
		NSArray *newPositions = [integrator integrateAlongDimension: tDim];
        
        //NSArray *newPositions = [integrator stepForwardToTime: maxTime];
        GLFunction *x = newPositions[0];
        GLFunction *y = newPositions[1];
        
        GLScalar *meanSquareSeparation = [[[x times: x] plus: [y times: y]] mean: 1];
        
		// kappa = 1/2 d/dt X^2 in *each* direction, hence 1/4.
        GLFloat kappaDeduced = (0.25)*(meanSquareSeparation.pointerValue[tDim.nPoints-1]-meanSquareSeparation.pointerValue[0])/maxTime;
        NSLog(@"kappa: %f, actual kappa: %f", kappa, kappaDeduced);
		
		GLFunction *u = [x diff: @"t"];
		GLFunction *v = [y diff: @"t"];
		
		GLScalar *uu = [[u times: u] mean];
		GLScalar *vv = [[v times: v] mean];
		GLScalar *uv = [[u times: v] mean];
		
		// The decorrelation time is the interval itself, dt, so kappa = v^2 dt.
		NSLog(@"(uu, vv, uv)=(%g, %g, %g)", (uu.pointerValue[0])*(timeStep), (vv.pointerValue[0])*(timeStep), (uv.pointerValue[0])*(timeStep));
    }
    return 0;
}

