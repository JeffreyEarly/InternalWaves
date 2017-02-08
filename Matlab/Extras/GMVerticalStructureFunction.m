latitude = 33;
f0 = 2 * 7.2921E-5 * sin( latitude*pi/180 );

Nz = 256;
L = 5000;
z = abs(L/2)*(cos(((0:Nz-1)')*pi/(Nz-1))+1) - L;

j_star = 3;
L_gm = 1.3e3; % thermocline exponential scale, meters
invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
E_gm = 6.3e-5; % non-dimensional energy parameter
E = L_gm*L_gm*L_gm*invT_gm*invT_gm*E_gm;
N0 = invT_gm;
g = 9.81;

H = (j_star+(1:3000)).^(-5/2);
H_norm = 1/sum(H);

j_max = Nz/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Exact structure function
%
rho0 = 1025;
rho = rho0*(1 + L_gm*N0*N0/(2*g)*(1 - exp(2*z/L_gm)));
[F, G, h, N2] = InternalWaveModesFromDensityProfile_Spectral( rho, z, z, 0.0, latitude, 'const_G_norm' );

PhiExact = 0*z;
GammaExact = 0*z;
for j = 1:j_max
    PhiExact = PhiExact + (1/h(j))*F(:,j).^2*H_norm*(j_star+j).^(-5/2); % the 1/h converts it to the const_F_norm
    GammaExact = GammaExact + G(:,j).^2*H_norm*(j_star+j).^(-5/2);
end
PhiExact = (L_gm*N0./sqrt(N2)).*PhiExact;
GammaExact = (1/g) * L_gm*N0*sqrt(N2).*GammaExact;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% WKB structure function
%
wkb_integral = @(arg) L_gm*( sqrt( N0*N0*exp(2*arg/L_gm) - f0*f0) - f0*atan(sqrt( N0*N0*exp(2*arg/L_gm) - f0*f0)/f0));
d = wkb_integral(0) - wkb_integral(-L);
xi = (wkb_integral(z) - wkb_integral(-L)) / d; 

PhiWKB = 0*z;
GammaWKB = 0*z;
for j = 1:j_max
    PhiWKB = PhiWKB + 2 * (cos(j*pi*xi)).^2*H_norm*(j_star+j).^(-5/2);
    GammaWKB = GammaWKB + 2*(sin(j*pi*xi)).^2*H_norm*(j_star+j).^(-5/2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Figure!
%
figure
subplot(1,2,1)
plot(PhiWKB,z), hold on, plot(PhiExact,z)
xlim([0 1.1*max(PhiExact)])
subplot(1,2,2)
plot(GammaWKB,z), hold on, plot(GammaExact,z)
xlim([0 1.1*max(GammaWKB)])