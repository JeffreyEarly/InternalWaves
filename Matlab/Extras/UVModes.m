N = 16;
f_bar = zeros(N,1);
j = 15;
f_bar(j+1) = 1;


f_tilde = cat(1, f_bar(1), 0.5*f_bar(2:N), 0, 0.5*f_bar(N:-1:2));


% f = ifct(f_bar);
% x = linspace(0,N,N);
% f_test = cos(2*pi*x/N);

f = real(fft(f_tilde));
f = f(1:N);
f_alt = 2*sqrt(2)*mirt_idctn(f_bar); % mirt_idctn uses DCT-II, we use DCT-I
% f_alt = ifft(f_tilde,'symmetric');
% f_alt = 2*N*f_alt(1:N);
Lz = N;
dz = Lz/N;
x = dz*(0:N-1)';
f_test = cos(j*pi*x/N);


figure
subplot(2,1,1)
plot(x,f_test, 'LineWidth', 3.0, 'Color', 'blue'), hold on
plot(x,f,'LineWidth', 1.0, 'Color', 'white', 'Marker', 'o', 'MarkerSize',6, 'MarkerFaceColor', 'black')
plot(x,f_alt)
title('(u,v)-mode')

f = real(fft(f_tilde));
f = f(1:N+1);
plot(dz*(0:N)',f)


% Now let's test W!
f_bar = zeros(N,1);
f_bar(j+1) = 1;

f_tilde = sqrt(-1)*cat(1, f_bar(1), 0.5*f_bar(2:N), 0, -0.5*f_bar(N:-1:2));

f = real(fft(f_tilde));
f = f(1:N);
Lz = N;
dz = Lz/N;
x = dz*(0:N-1)';
f_test = sin(j*pi*x/N);


subplot(2,1,2),
plot(x,f_test, 'LineWidth', 3.0, 'Color', 'blue'), hold on
plot(x,f,'LineWidth', 1.0, 'Color', 'white', 'Marker', 'o', 'MarkerSize',6, 'MarkerFaceColor', 'black')
title('(w)-mode')

g = cat(1,f,0,-f(N:-1:2));
g_bar = ifft(g);
g_bar_tilde = 2*imag(g_bar(2:N+1));
