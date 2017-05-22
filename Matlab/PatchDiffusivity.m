function [r2, kappa_r] = PatchDiffusivity(t,x,y,xStride,yStride)
% size(x) = [length(t) numDrifters]
% xStride tells us the index spacing between consecutive drifters in x
% yStride tells us the index spacing between consecutive drifters in y
% By assumption, one of these strides is 1.

nDrifters = size(x,2);
if xStride == 1
    Nx = yStride;
    Ny = nDrifters/Nx;
elseif yStride == 1
    Ny = xStride;
    Nx = nDrifters/Ny;
else
    error('One of the strides must be 1, by assumption');
end

divisors = 2.^((1:floor(log2(min(Nx,Ny))))-1);

nReps = sum(divisors.^2);
r2 = zeros(nReps,1);
kappa_r = zeros(nReps,1);
iteration = 1;
for divisor=divisors
    nPatchX = floor(Nx/divisor); % number of particles in x-direction in the patch
    nPatchY = floor(Ny/divisor);
    
    for i=1:divisor
       for j=1:divisor
           xInd = (i-1)*nPatchX + (1:nPatchX);
           yInd = (j-1)*nPatchY + (1:nPatchY);
           linearIndices = zeros(nPatchX*nPatchY,1);
           iIndex = 1;
           for ix = xInd
              for iy = yInd
                  linearIndices(iIndex) = (ix-1)*xStride + (iy-1)*yStride + 1;
                  iIndex = iIndex+1;
              end
           end
           
           if (max(linearIndices)>nDrifters)
               fprintf('too big.');
           end
           
           [X2,Y2] = MeanSquareSeparation(x(:,linearIndices),y(:,linearIndices));
           D2 = X2+Y2;
           [p,~,mu]=polyfit(t,D2,1);
           
           kappa_r(iteration) = (p(1)/mu(2))/4;
           r2(iteration) = mean(D2,1);
           iteration = iteration+1;
       end
    end
end

end