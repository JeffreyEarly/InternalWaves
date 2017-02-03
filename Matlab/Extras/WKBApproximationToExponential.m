D = 4000;
z = linspace(-D,0,100)';
b = 1300;

N = exp(z/b);
% figure, plot(N,z)

omega = 0.5;
p = sqrt(N.*N-omega*omega);
z_turning_point = b*log(omega);

z_upper = linspace(z_turning_point,0,1000)';
N_upper = exp(z_upper/b);
p_upper = sqrt(N_upper.*N_upper-omega*omega);

psi_upper = sin(-cumtrapz(p_upper)/10 + pi/4)./sqrt(p_upper);

z_lower = linspace(-D,z_turning_point-500,1000)';
N_lower = exp(z_lower/b);
p_lower = sqrt(omega*omega - N_lower.*N_lower);

psi_lower = exp(cumtrapz(p_lower))./sqrt(p_lower);

figure, plot(psi_upper,z_upper)
hold on
plot(psi_lower,z_lower)
