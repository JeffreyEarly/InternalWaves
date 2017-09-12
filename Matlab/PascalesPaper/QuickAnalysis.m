load('particle_data_nonlinear.mat');

% tindices=1:floor(length(timevar)/2);
% tindices=floor(length(timevar)/2):length(timevar);
tindices = 1:length(timevar);

t = timevar(tindices);
x = xvar(:,tindices)';
y = yvar(:,tindices)';
z = zvar(:,tindices)';

[~,~,q,r] = CenterOfMass(x,y);
figure
plot(q,r)

figure
scatter(q(1,:),r(1,:))

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First let's look at z
[Mqq, Mrr] = SecondMomentMatrix( z, z );

D2 = Mqq;
[p,~,mu]=polyfit(t,D2,1);
kappa_fit = (p(1)/mu(2))/2; % Factor 2 this time, b/c 1D
intercept = p(2)-p(1)*mu(1)/mu(2);

kappa_endpoint = 0.5*(D2(end)-D2(1))/(t(end)-t(1));

figure
subplot(1,2,1)
plot(t, D2)
hold on
plot(t,2*kappa_fit*t + intercept)
title(sprintf('Second moment (kappa_{fit}, kappa_{endpoint})=(%.2g,%.2g) m^2/s', kappa_fit, kappa_endpoint))
xlim([min(t) max(t)])

[minD, maxD, theta] = SecondMomentMatrix( x, y, 'eigen' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now let's check out x and y
subplot(1,2,2)
plot(t,[minD maxD])

D2 = minD+maxD;

hold on
plot(t, D2)

% Now let's fit a line through the second moment
[p,~,mu]=polyfit(t,D2,1);
kappa_fit = (p(1)/mu(2))/4;
intercept = p(2)-p(1)*mu(1)/mu(2);

% Alternatively, we can just use the endpoints
kappa_endpoint = 0.25*(D2(end)-D2(1))/(t(end)-t(1));

hold on
plot(t,4*kappa_fit*t + intercept)
title(sprintf('Second moment (kappa_{fit}, kappa_{endpoint})=(%.2g,%.2g) m^2/s', kappa_fit, kappa_endpoint))
xlim([min(t) max(t)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Okay, what does the pairwise diffusivity look like for x & y

[r2, kappa_r_slope ] = PairwiseRelativeDiffusivity( t, x, y, 'slope' );
figure
scatter(sqrt(r2), kappa_r_slope)

theBins = 100:500:7000;
figure
[xMean, yMean, yStdErr] = histogramWithErrorbars(sqrt(r2),kappa_r_slope,theBins);