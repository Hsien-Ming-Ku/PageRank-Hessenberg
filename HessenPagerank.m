function [q,iter,res] = HessenPagerank(A,q,dim,itmax,tol)
% The restarted and refined Hessenberg process for computing PageRank,
% which will be proposed by X.-M. Gu, S.-L. Lei, Z.-L. Shen, K. Zhang, C.
% Wen (2018, in preparation). 
% -----------------------------------------------
% written by Dr. Xian-Ming Gu, who recently works at the School of Economic
% Mathematics, SWUFE.
% Date: 24 July, 2018, the Visting Scholar at University of Macau
%       4 January, 2019 (modified some parts), at SWUFE
%       7 January, 2019 (modified some parts), at SWUFE
% Email: guxianming@live.cn, guxm@swufe.edu.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------- Check the input parameter A ----------------
if isa(A, 'double')
   matrix_or_function = 1;
elseif isa(A, 'function_handle')
   matrix_or_function = 2;
else
   error('A is neither a matrix or a function handle !!');
end % end if
% end
% -------------------------------------------------
n = length(q); it = 0; 
m = dim;
p = zeros(n,1); res = zeros(itmax+1,1);
L = zeros(n,m+1);
H = zeros(m + 1,m);
%%%%%%%%%%%%%%%%%%%%%%%%
[~,i0] = max(abs(q)); 
beta = q(i0); 
p(1) = i0; L(:,1) = q/beta; 
% qt = q - A(q);
if matrix_or_function == 1
   qt = q - A*q;
else
   qt = q - A(q);   % Other cases will result in the error.
end % end if
nres = norm(qt,1)/norm(q,1); res(1) = nres;
%%%%%%%%%%%%%%%%%%%
It = zeros(m + 1,m);
It(1:m,1:m) = eye(m);
%-----------------
while (it < itmax) && (nres > tol)
%     fprintf(' %3.0f %6.4e  \n',it,nres);
    it = it + 1;
    for k = 1:m
%         u = A*L(:,k);
        if matrix_or_function == 1
           u = A*(L(:,k));
        else
           u = A(L(:,k));   % Other cases will result in the error.
        end % end if
%----------------- the inner for-loop: LU-like process ------------------
        for j = 1 : k
            pj = p(j); H(j,k) = u(pj); 
            u = u - H(j,k)*L(:,j);
        end 
        if ( k < n )
            [ibeta,i0] = max(abs(u)); H(k+1,k) = u(i0); p(k+1) = i0; 
            if ( abs(ibeta) > 0 )
                L(:,k+1) = u/H(k+1,k);
            else
                L(:,k+1) = zeros(n,1);
                disp('Happy break-down');
            end
        end
    end
    H1 = H - It;   % H_{m+1,m} - [I_m;0] = U*\Sigma*(V^T);
    [U,S,V] = svd(H1);
    q = (L(:,1:m))*V(:,m);
    deta = diag(S);
    deta = deta(end);
    qt = deta*(L*U(:,m)); 
%   qt = q - A(q);
    nres = norm(qt,1)/norm(q,1); res(it + 1) = nres;
    if nres < tol
       break;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (it < itmax) && (nres > tol)
        [~,i0] = max(abs(q)); 
        beta = q(i0);
        p(1) = i0; L(:,1) = q/beta;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
end
iter = it;
% fprintf(' %3.0f %6.4e  \n',it,nres);
res = res(1:iter+1);
end
