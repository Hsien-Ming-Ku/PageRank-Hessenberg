function [q,iter,res,mv] = IDRsPagerank(A,Atrace,s,dim,q,itmax,it_pow,tol)
%%%%%
%%%%%%%%
%%%%%%%%%%%
n = length(q);
w0 = q; w0 = w0/norm(w0);
m = dim;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = rand(s,n);
j = floor(m/(s+1));
mv = 0;
oms_i = ((Atrace/n) + .0001*randn(j,1));
%%%%%%%%%%%%%%%%%%%
It = zeros(m+1,m);
It(1:m,1:m) = eye(m);
it = 0;
qt = q - A(q);
nres = norm(qt,1);
res = zeros(itmax+1,1);
res(1) = nres;
% rev
while (it < itmax) && (nres > tol)
%     fprintf(' %3.0f %6.4e  \n',it,nres);
    it = it + 1;
    %%% Compute the Hessenberg relation: A*L_{m} = L_{m+1}*H;
    [W,H,oms] = idr_process_standard(A,s,m,w0,P,oms_i);
%     cond(W),
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    H1 = H - It;   % H_{m+1,m} - [I_m;0] = U*\Sigma*(V^T);
    [U,S,V] = svd(H1);
    q = (W(:,1:m))*V(:,m);
    rt = S(m,m)*(W*U(:,m)); 
%     r = A(q1) - q1;
    nres = norm(rt,1)/norm(q,1);
    res(it+1) = nres;
    if (it < itmax) && (nres > tol)
%         w0 = q; w0 = w0/norm(w0);  % The first choice
        %%% -----------------------------------------
%         P = rand(s,n);
        oms_i = ((Atrace/n) + .0001*randn(j,1));
%        oms_i = oms;
        for k = 1:it_pow
            mv = mv + 1;
            q = q/norm(q,1);
%         w0 = q; w0 = w0/norm(w0);  % The first choice
            q = A(q);
        end
        w0 = q; w0 = w0/norm(w0);
    end
iter = it;
res = res(1:iter + 1);
q = q/norm(q,1);
% fprintf(' %3.0f %6.4e  \n',it,nres);
end