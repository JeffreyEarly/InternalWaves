% Given a vector y where size(y)=[N 8] such that we have N floats in three
% different experiments. They are ordered [x_float y_float z_float x_diff
% y_diff z_diff x_drifter y_drifter], and should therefore return [u v w u
% v w u v].
function flux = FluxForLagrangianErrorExperiment(t,y,wavemodel)
    [u_exact,v_exact,w_exact] = wavemodel.VelocityAtTimePosition(t,y(:,1),y(:,2),y(:,3),'exact');
    [u_linear,v_linear,w_linear] = wavemodel.VelocityAtTimePosition(t,y(:,4),y(:,5),y(:,6),'linear');
    [u_spline,v_spline,w_spline] = wavemodel.VelocityAtTimePosition(t,y(:,7),y(:,8),y(:,9),'spline');
    
    flux = cat(2,u_exact,v_exact,w_exact,u_linear,v_linear,w_linear,u_spline,v_spline, w_spline);
end