j_star = 3;
H = (j_star+(1:1024)).^(-5/2);
H_norm = 1/sum(H);
H = H_norm*(j_star+(1:1024)).^(-5/2);

L = 100;
z = linspace(0,L,65)';

j_max = 3;

F2 = 0*z;
G2 = 0*z;
for j = 1:j_max
    F2 = F2 + (cos(z*j*pi/L)*H_norm*(j_star+j).^(-5/2)).^2;
    G2 = G2 + (sin(z*j*pi/L)*H_norm*(j_star+j).^(-5/2)).^2;
end

figure
subplot(1,2,1)
plot(F2,z)
xlim([0 1.1*max(F2)])
subplot(1,2,2)
plot(G2,z)
xlim([0 1.1*max(G2)])


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