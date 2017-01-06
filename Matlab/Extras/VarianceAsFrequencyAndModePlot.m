latitude = 33;
f0 = 2 * 7.2921E-5 * sin( latitude*pi/180 );

j_star = 3;
L_gm = 1.3e3; % thermocline exponential scale, meters
invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
E_gm = 6.3e-5; % non-dimensional energy parameter
E = L_gm*L_gm*L_gm*invT_gm*invT_gm*E_gm;
N0 = invT_gm;

H = (j_star+(1:3000)).^(-5/2);
H_norm = 1/sum(H);
B_norm = 1/atan(sqrt(N0*N0/(f0*f0)-1));

GM2D_int = @(omega0,omega1,j) E*H_norm*B_norm*((j+j_star).^(-5/2))*(atan(f0/sqrt(omega0*omega0-f0*f0)) - atan(f0/sqrt(omega1*omega1-f0*f0)));
GM2D_uv_int = @(omega0,omega1,j) E*H_norm*B_norm*((j+j_star).^(-5/2))*( f0*sqrt(omega1*omega1-f0*f0)/(2*omega1*omega1) - (3/2)*atan(f0/sqrt(omega1*omega1-f0*f0)) - f0*sqrt(omega0*omega0-f0*f0)/(2*omega0*omega0) + (3/2)*atan(f0/sqrt(omega0*omega0-f0*f0)));
GM2D_w_int = @(omega0,omega1,j) E*H_norm*B_norm*((j+j_star).^(-5/2))*( f0*sqrt(omega1*omega1-f0*f0) - f0*f0*atan(f0/sqrt(omega1*omega1-f0*f0)) - f0*sqrt(omega0*omega0-f0*f0) + f0*f0*atan(f0/sqrt(omega0*omega0-f0*f0)));
GM2D_zeta_int = @(omega0,omega1,j) E*H_norm*B_norm*((j+j_star).^(-5/2))*( ((omega1*omega1-f0*f0)^(3/2))/(2*f0*omega1*omega1) - (1/2)*atan(f0/sqrt(omega1*omega1-f0*f0)) - sqrt(omega1*omega1-f0*f0)/(2*f0) - ((omega0*omega0-f0*f0)^(3/2))/(2*f0*omega0*omega0) + (1/2)*atan(f0/sqrt(omega0*omega0-f0*f0)) + sqrt(omega0*omega0-f0*f0)/(2*f0) );

omegaAxis = linspace(f0,N0,1000)';
modeAxis = (1:35)';

HKE = zeros(length(omegaAxis),length(modeAxis));
for j = modeAxis
   for iOmega = 1:(length(omegaAxis)-1)
       HKE(iOmega,j) = GM2D_uv_int(omegaAxis(iOmega),omegaAxis(iOmega+1),j);
   end
end

figure
pcolor(omegaAxis,modeAxis,log10(HKE/E)')
shading flat