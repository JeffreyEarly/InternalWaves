L = 1;
N = 8;
delta = L/N;
n = (0:(N-1))';
x = n*delta + delta/2;
f = n/(2*N*delta);

a = 1 + 3*cos(4*pi*x);

mirt_dctn(a)