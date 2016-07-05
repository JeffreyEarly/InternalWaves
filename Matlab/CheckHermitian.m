function CheckHermitian(A)
M = size(A,1);
N = size(A,2);
K = size(A,3);

for k=1:K
   for i=M:-1:1
       for j=N:-1:1
           ii = mod(M-i+1, M) + 1;
           jj = mod(N-j+1, N) + 1;
           if ( A(i,j,k) ~= conj(A(ii,jj,k)) )
              fprintf('i,j,k = (%d,%d,%d)\n',i,j,k);  
           end
       end
   end
end

end