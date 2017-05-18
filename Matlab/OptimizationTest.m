clear

N = 1e6; % N points
M = 64; % M waves

x = rand(N,1);
k = rand(1,M);
A = rand(1,M);

tic
a = sum(A.*cos(x*k),2);
toc

tic
a1 = cos(x*k)*A'; %  [N M] * [M 1]
toc

tic
b = sum(A.*cos(bsxfun(@times,x,k)),2);
toc

c = zeros(size(x));
tic
for i=1:M
   c = c + A(i)*cos(k(i)*x); 
end
toc