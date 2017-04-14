% Given a vector y where size(y)=[N 8] such that we have N floats in three
% different experiments. They are ordered [x_float y_float z_float x_diff
% y_diff z_diff x_drifter y_drifter], and should therefore return [u v w u
% v w u v].
function flux = FluxForFloatDiffusiveDrifter(t,y,z_drifter,wavemodel)
    [u_float,v_float,w_float] = wavemodel.VelocityAtTimePosition(t,y(:,1),y(:,2),y(:,3));
    [u_diffusive_float,v_diffusive_float,w_diffusive_float] = wavemodel.VelocityAtTimePosition(t,y(:,4),y(:,5),y(:,6));
    [u_drifter,v_drifter] = wavemodel.VelocityAtTimePosition(t,y(:,7),y(:,8),z_drifter);
    
    flux = cat(2,u_float,v_float,w_float,u_diffusive_float,v_diffusive_float,w_diffusive_float,u_drifter,v_drifter);
end