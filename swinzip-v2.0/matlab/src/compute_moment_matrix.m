function [M, U]=compute_moment_matrix(pos,k,KK,d,seli,k2)

% compute_moment_matrix - computes moment matrices based on a polynomial basis.
  
% build moment matrix       
kk = length(seli);                   % nbr of point in this bin, ie. 2k^2 in regular case
M = zeros( kk, 2*k2 );
for ii=1:kk    
    for j=1:k2
       M(ii,j)=1;
       for dd=1:d
          M(ii,j)=M(ii,j)*pos(dd,seli(ii))^KK(j,d-dd+1);
       end
    end
end
   
% orthogonalize
[U,R] = qr(M);
U = transpose(U);