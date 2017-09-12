figure

file = cell(2,1);
file{1} = 'particle_data_linear.mat';
file{2} = 'particle_data_nonlinear.mat';

for i = 1:2
    load(file{i});
    
    % tindices=1:floor(length(timevar)/2);
    % tindices=floor(length(timevar)/2):length(timevar);
    tindices = 1:length(timevar);
    
    t = timevar(tindices);
    x = xvar(:,tindices)';
    y = yvar(:,tindices)';
    z = zvar(:,tindices)';
    
    [r2, kappa_r_slope ] = PairwiseRelativeDiffusivity( t, x, y, 'slope' );
    
%     scatter(sqrt(r2), kappa_r_slope)
    
    theBins = 100:500:7000;
    theBins = 10.^(linspace(log10(100),log10(7000),15));
    hold on
    [xMean, yMean, yStdErr] = histogramWithErrorbars(sqrt(r2),kappa_r_slope,theBins);
end

legend('linear','nonlinear')