function [L,H] = Gene_Hessenberg(A,v,m)
% The generalized Hessenberg (GH) procedure with pivoting strategy for 
% generate the Krylov subspace.
% Input:
%   A    square matrix (n by n)
%   v    initial vector
%   m    number of iterations
% Output: 
%   L    non-orthonormal basis of Krylov space (n by m+1)
%   H    upper Hessenberg matrix, A*L(:,1:m) = L*H (the size of H is (m+1)-by m)
% References:
% 1. M. Heyouni, H. Sadok, A new implementation of the CMRH method for 
%    solving dense linear systems, J. Comput. Appl. Math., Volume 213, 
%    Issue 2, 2008, pp. 387-399. (See Algorithm 2)
% 2. H. Sadok, CMRH: A new method for solving nonsymmetric linear systems 
%    based on the Hessenberg reduction algorithm, Numer. Algorithms, 
%    Volume 20, Issue 4, 1999, pp. 303-321.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The codes are developed and written by Dr. Xian-Ming Gu, who is an 
% Assitant Professor of SWUFE.
% Contact information:
% e-mail: guxianming@live.cn, x.m.gu@rug.nl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Date: 03-07-2018. 17:43 (The University of Macau)

% Begin to implement the GH-process with pivoting strategy for A

n = length(v);       % the size of test matrix.
L = zeros(n,m + 1);  % initialize the  
H = zeros(m + 1,m);  % the size of test matrix.
p = 1:n;             % Auxiliary vector for order.
[~,i0] = max(abs(v));
alpha = v(i0);
L(:,1) = v/alpha;
%swap p(1) and p(i0)
tmp = p(1);
p(1) = p(i0);
p(i0) = tmp;
%s = r(i0)*e1;
for k = 1:m
%     u = A*L(:,k);   % Matrix-vector product.
   if isa(A, 'double')
      u = A*(L(:,k));
   elseif isa(A, 'function_handle')
      u = A(L(:,k));
   else
      error('A is neither a matrix or a function handle!!!');
   end % end if
%    matvec = matvec + 1;
   for j = 1:k
       H(j,k) = u(p(j));
       u = u - H(j,k)*L(:,j);
   end
   tp = zeros(n-k,1);
   if (k < n && ~isequal(u,zeros(n,1)))
      for i = k + 1:n
          tp(i-k) = u(p(i)); % Vector of u(p(k+1)),...,u(p(n))
      end
      [~,ind_p]= max(abs(tp)); % Maximum value of the vector abs(tp)
      i0 = k + ind_p;
      H(k + 1,k) = u(p(i0));
      L(:,k + 1) = u/H(k + 1,k);
      % Swap contents
      tmp = p(k+1);
      p(k + 1) = p(i0);
      p(i0) = tmp;
   else
      H(k + 1,k) = 0;  % Stop.
   end
end