latitude = 33;
f0 = 2 * 7.2921E-5 * sin( latitude*pi/180 );

j_star = 3;
L_gm = 1.3e3; % thermocline exponential scale, meters
invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
E_gm = 6.3e-5; % non-dimensional energy parameter
E = L_gm*L_gm*L_gm*invT_gm*invT_gm*E_gm;
N0 = invT_gm;

H = (j_star+(1:1024)).^(-5/2);
H_norm = 1/sum(H);
B_norm = 1/atan(sqrt(N0*N0/(f0*f0)-1));


GM2D_uv_int = @(omega0,omega1,j) B_norm*H_norm*E*((j+j_star).^(-5/2))*( f0*sqrt(omega1*omega1-f0*f0)/(2*omega1*omega1) - (3/2)*atan(f0/sqrt(omega1*omega1-f0*f0)) - f0*sqrt(omega0*omega0-f0*f0)/(2*omega0*omega0) + (3/2)*atan(f0/sqrt(omega0*omega0-f0*f0)));


L = 5000;
z = linspace(0,L,65)';

H = H_norm*(j_star+(1:1024)).^(-5/2);

j_max = 300;

F2 = 0*z;
G2 = 0*z;
for j = 1:j_max
    F2 = F2 + (cos(z*j*pi/L)).^2*H_norm*(j_star+j).^(-5/2);
    G2 = G2 + (sin(z*j*pi/L)).^2*H_norm*(j_star+j).^(-5/2);
end

figure
subplot(1,2,1)
plot(F2,z)
xlim([0 1.1*max(F2)])
subplot(1,2,2)
plot(G2,z)
xlim([0 1.1*max(G2)])





return

F2 = 0*z;
G2 = 0*z;
for j = 1:j_max
    F2 = F2 + cos(z*j*pi/L)*H_norm*(j_star+j).^(-5/2);
    G2 = G2 + sin(z*j*pi/L)*H_norm*(j_star+j).^(-5/2);
end
F2 = F2.^2;
G2 = G2.^2;

figure
subplot(1,2,1)
plot(F2,z)
xlim([0 1.1*max(F2)])
subplot(1,2,2)
plot(G2,z)
xlim([0 1.1*max(G2)])