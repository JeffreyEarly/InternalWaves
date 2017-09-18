method = 'slope';

file = cell(2,1);
file{1} = 'particle_data_linear.mat';
file{2} = 'particle_data_nonlinear.mat';

load(file{1});

% tindices=1:floor(length(timevar)/2);
% tindices=floor(length(timevar)/2):length(timevar);
tindices = 1:length(timevar);

t = timevar(tindices);
x = xvar(:,tindices)';
y = yvar(:,tindices)';
z = zvar(:,tindices)';

[r2_h_lin, kappa_r_h_lin, r2_z_lin, kappa_r_z_lin, uv_correlation_lin, w_correlation_lin] = PairwiseRelativeDiffusivityWithZ( t, x, y, z, 'slope' );

load(file{2});

t = timevar(tindices);
x = xvar(:,tindices)';
y = yvar(:,tindices)';
z = zvar(:,tindices)';

[r2_h_nl, kappa_r_h_nl, r2_z_nl, kappa_r_z_nl, uv_correlation_nl, w_correlation_nl] = PairwiseRelativeDiffusivityWithZ( t, x, y, z, 'slope' );


horizontalBins = 10.^(linspace(log10(100),log10(7000),15));
verticalBins = 10.^(linspace(log10(0.01),log10(20),15));
verticalDiffusivityBins = 10.^(linspace(log10(1e-7),log10(5e-5),15));

figure

% First plot horizontal diffusivity vs distance
subplot(2,2,1)
[xMean, yMean, yStdErr] = histogramWithErrorbars(sqrt(r2_h_lin),kappa_r_h_lin,horizontalBins);
hold on
[xMean, yMean, yStdErr] = histogramWithErrorbars(sqrt(r2_h_nl),kappa_r_h_nl,horizontalBins);
legend('linear','nonlinear')
xlabel('lateral separation (m)')
ylabel('lateral diffusivity (m^2/s)')

subplot(2,2,2)
[xMean, yMean, yStdErr] = histogramWithErrorbars(sqrt(r2_h_lin),kappa_r_z_lin,horizontalBins);
hold on
[xMean, yMean, yStdErr] = histogramWithErrorbars(sqrt(r2_h_nl),kappa_r_z_nl,horizontalBins);
legend('linear','nonlinear')
xlabel('horizontal separation (m)')
ylabel('vertical diffusivity (m^2/s)')

subplot(2,2,3)
[xMean, yMean, yStdErr] = histogramWithErrorbars(sqrt(r2_z_lin),kappa_r_h_lin,verticalBins);
hold on
[xMean, yMean, yStdErr] = histogramWithErrorbars(sqrt(r2_z_nl),kappa_r_h_nl,verticalBins);
legend('linear','nonlinear')
xlabel('vertical separation (m)')
ylabel('lateral diffusivity (m^2/s)')

subplot(2,2,4)
[xMean, yMean, yStdErr] = histogramWithErrorbars(sqrt(r2_z_lin),kappa_r_z_lin,verticalBins);
hold on
[xMean, yMean, yStdErr] = histogramWithErrorbars(sqrt(r2_z_nl),kappa_r_z_nl,verticalBins);
legend('linear','nonlinear')
xlabel('vertical separation (m)')
ylabel('vertical diffusivity (m^2/s)')


figure
histogramWithErrorbars(kappa_r_z_lin, kappa_r_h_lin,verticalDiffusivityBins);
hold on
histogramWithErrorbars(kappa_r_z_nl, kappa_r_h_nl, verticalDiffusivityBins);
legend('linear','nonlinear')
xlabel('vertical diffusivity (m^2/s)')
ylabel('lateral diffusivity (m^2/s)')
xlog
