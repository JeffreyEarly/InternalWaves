% Given a vector y where size(y)=[N 8] such that we have N floats in three
% different experiments. They are ordered [x_float y_float z_float x_diff
% y_diff z_diff x_drifter y_drifter], and should therefore return [u v w u
% v w u v].
function flux = FluxForFloatDiffusiveDrifter(t,y,z_drifter,wavemodel,method)
    nFloats = size(y,1);
    floatIndices = 1:nFloats;
    diffusiveFloatIndices = floatIndices + nFloats;
    drifterIndices = diffusiveFloatIndices + nFloats;
    
    [u,v] = wavemodel.VelocityAtTimePosition(t,cat(1,y(:,1),y(:,4),y(:,7)), cat(1,y(:,2),y(:,5),y(:,8)), cat(1,y(:,3),y(:,6),z_drifter), method);
    w = wavemodel.VerticalVelocityAtTimePosition(t,cat(1,y(:,1),y(:,4)), cat(1,y(:,2),y(:,5)), cat(1,y(:,3),y(:,6)), method);
    
    flux = cat(2,u(floatIndices),v(floatIndices),w(floatIndices),u(diffusiveFloatIndices),v(diffusiveFloatIndices),w(diffusiveFloatIndices),u(drifterIndices),v(drifterIndices));
end