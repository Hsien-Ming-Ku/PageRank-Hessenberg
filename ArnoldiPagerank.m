function [q,iter,res] = ArnoldiPagerank(A,q,dim,itmax,tol,it_pow)
% The restarted refined Arnoldi process for computing PageRank, which is 
% proposed by G. H. Golub, C. Greif (2006), BIT Numer. Math., 46: 759-771.
% -----------------------------------------------------------------------
% written by Dr. Xian-Ming Gu, who recently works at the School of Economic
% Mathematics, SWUFE.
% Date: 24 July, 2018, the Visting Scholar at University of Macau
%       4 January, 2019 (modified some parts), at SWUFE
%       7 January, 2019 (modified some parts), at SWUFE
% Email: guxianming@live.cn, guxm@swufe.edu.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------- Check the input parameter A ----------------
if isa(A, 'double')
   matrix_or_function = 1;
elseif isa(A, 'function_handle')
   matrix_or_function = 2;
else
   error('A is neither a matrix or a function handle !!');
end % end if
% -------------------------------------------------
n = length(q); it = 0; 
m = dim;
res = zeros(itmax+1,1);
L = zeros(n,m + 1);
H = zeros(m + 1,m);
% -----------------------------------------------
beta = norm(q); 
L(:,1) = q/beta; 
% qt = q - A(q);
if matrix_or_function == 1
   qt = q - A*q;
else
   qt = q - A(q);   % Other cases will result in the error.
end % end if
nres = norm(qt,1)/norm(q,1); res(1) = nres;
% --------------------------------------------------------
It = zeros(m+1,m);
It(1:m,1:m) = eye(m);
% ---------------------------------------------------------
while (it < itmax) && (nres > tol)
%     fprintf(' %3.0f %6.4e  \n',it,nres);
    it = it + 1;
    for k = 1:m
%       u = A*L(:,k);
        if matrix_or_function == 1
           u = A*(L(:,k));
        else
           u = A(L(:,k));  % Other cases will result in the error.
        end % end if
%----------------- the inner for-loop: QR-MGS process ------------------
        for j = 1 : k
%             temLj = L(:,j);
            H(j,k) = u'*L(:,j); 
            u = u - H(j,k)*L(:,j);
        end
        hlast = norm(u);
        if hlast < eps
           disp('Happy break-down');
           break;
        end
        H(k + 1,k) = hlast; 
        L(:,k + 1) = u/H(k + 1,k);
    end
    H1 = H - It;   % H_{m+1,m} - [I_m;0] = U*\Sigma*(V^T);
    [U,S,V] = svd(H1);
    q = L(:,1:m)*(V(:,m));
    deta = diag(S);
    deta = deta(end);
    qt = deta*(L*U(:,m));
%     qt = q - A(q);
%     qt = qt/norm(qt,1);
    nres = norm(qt,1)/norm(q,1); res(it + 1) = nres;
    if nres < tol
       break;
    end
%     norm(q),
%%%%%%%%%%%%%%%%%%%%%
    if (it < itmax) && (nres > tol)
        for k = 1:it_pow
%             mv = mv + 1;
            q = q/norm(q,1);
%         w0 = q; w0 = w0/norm(w0);  % The first choice
            q = A(q);
        end
        beta = norm(q);
        L(:,1) = q/beta;
    end
end
iter = it;
res = res(1:iter + 1);
% fprintf(' %3.0f %6.4e  \n',it,nres);
end