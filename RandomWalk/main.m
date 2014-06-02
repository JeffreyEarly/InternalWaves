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
        
        GLScalar *meanSquareSeparation = [[[[x times: x] plus: [y times: y]] mean: 2] mean: 1];
        
        GLFloat kappaDeduced = (0.25)*(meanSquareSeparation.pointerValue[tDim.nPoints-1]-meanSquareSeparation.pointerValue[0])/maxTime;
        NSLog(@"kappa: %f, actual kappa: %f", kappa, kappaDeduced);
		
		GLFunction *u = [x diff: @"t"];
		GLFunction *v = [y diff: @"t"];
		
		GLScalar *uu = [[u times: u] mean];
		GLScalar *vv = [[v times: v] mean];
		GLScalar *uv = [[u times: v] mean];
		
		NSLog(@"(uu, vv, uv)=(%g, %g, %g)", (uu.pointerValue[0])*(timeStep/2.), (vv.pointerValue[0])*(timeStep/2.), (uv.pointerValue[0])*(timeStep/2.));
		
//		GLFunction *vx = [x diff:@"t"];
//		GLFunction *vy = [y diff:@"t"];
//		GLScalar *Vxx = [[vx times: vx] mean];
//		GLScalar *Vyy = [[vy times: vy] mean];
//		GLScalar *Vxy = [[vx times: vy] mean];
//		
//		NSLog(@"(kappa_xx, kappa_yy, kappa_xy)=(%g, %g, %g)", Vxx.pointerValue[0]*timeStep/2., Vyy.pointerValue[0]*timeStep/2., Vxy.pointerValue[0]*timeStep/2.);
		
    }
    return 0;
}

