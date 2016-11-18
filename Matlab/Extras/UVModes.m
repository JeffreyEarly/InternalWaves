N = 16;
f_bar = zeros(N,1);
j = 3;
f_bar(j+1) = 1;


f_tilde = cat(1, f_bar(1), 0.5*f_bar(2:N), 0, 0.5*f_bar(N:-1:2));
f_tilde = cat(1, f_bar(1), 0.5*f_bar(2:N), 0, 0.5*f_bar(N:-1:2));


% f = ifct(f_bar);
% x = linspace(0,N,N);
% f_test = cos(2*pi*x/N);

f = real(fft(f_tilde));
f = f(1:N);
Lz = N;
dz = Lz/N;
x = dz*(0:N-1)';
f_test = cos(j*pi*x/N);


  figure, plot(x,f,'LineWidth', 2)
  hold on, plot(x,f_test)


% Now let's test W!
j = 3;
f_bar = zeros(N,1);
f_bar(j+1) = 1;

f_tilde = sqrt(-1)*cat(1, f_bar(1), 0.5*f_bar(2:N), 0, -0.5*f_bar(N:-1:2));

f = real(fft(f_tilde));
f = f(1:N);
Lz = N;
dz = Lz/N;
x = dz*(0:N-1)';
f_test = sin(j*pi*x/N);


figure, plot(x,f,'LineWidth', 2)
hold on, plot(x,f_test)