function [xm,iter,res] = weighted_ArnoldiPagerank(A,q,tol,dim,itmax)
% The restarted adpative refined weighted Arnoldi process for computing 
% PageRank, which is proposed by J.-F. Yin, G.-J. Yin, M. Ng (2012), NLAA,
% 19: 73-85.
% ---------------------------------------------------------------------
% Written by Dr. Xian-Ming Gu, SWUFE, 2019-1-5: 14:17;
% via the help from Dr. Guojian Yin (HKUST) and Dr. Zhao-Li Shen (SICAU)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------- check the input parameter A ----------------------------
if isa(A, 'double')
   matrix_or_function = 1;
elseif isa(A, 'function_handle')
   matrix_or_function = 2;
else
   error('A is neither a matrix or a function handle !!');
end % end if
%%%%----------------------------------------------------------------
% initial settings
n = length(q); it = 0;
wt = ones(n,1); % weight
% sqrtn = sqrt(n);
x0 = q;         % initial guess
%----------- variables for the weighted Arnoldi process ------------------
res = zeros(itmax+1,1);
m = dim;
V = zeros(n,m+1);
H = zeros(m+1,m);
MI = zeros(m+1,m);
MI(1:m,1:m) = eye(m);
%----------- start the restart weighted Arnoldi process ------------------
if matrix_or_function == 1
   r0 = x0 - A*(x0);
else
   r0 = x0 - A(x0);
end 
%------------------------------------ 
nres = norm(r0,1)/norm(x0,1);
res(1) = nres;
V(:,1) = x0/normw(x0,wt); %%% -- indicate that weighted techniques is used
while (it < itmax) && (nres > tol) 
% beta=norm(r0,2);
      it = it + 1;
      for j = 1:m  
%           temV(j).v = V(:,j);
%           wttemV(j).v = wt.*temV(j).v;
          if matrix_or_function == 1
              w = A*V(:,j);
          else
              w = A(V(:,j));
          end
%---------------- the inner for-loop: QR-MGS process -------------------   
          for i = 1:j
%               temVi = wt.*V(:,i);
              H(i,j) = w'*(wt.*V(:,i));
              w = w - H(i,j)*V(:,i);
          end
          hlast = normw(w,wt);
          if hlast < eps
             disp('Happy break-down');
             break;       
          end
          H(j + 1,j) = hlast;   % *****important
          V(:,j + 1) = w/H(j+1,j);
      end
% temHm=H(1:m,:);
     [U,S,W] = svd(H - MI);
% U: m+1 * m+1
% S: m+1 * m
% W: m * m
% s=diag(S);
% deta=min(s);
% um=U(:,m+1);
     deta = diag(S);
     deta = deta(end);
%rm=deta*(V*um(1:m)+um(end)*vm1);
     rm = deta*(V*U(:,m));
     xm = V(:,1:m)*(W(:,m));
% rm=r_original;
     nres = norm(rm,1)/norm(xm,1);
     res(it + 1) = nres;
     if nres < tol
        break;
     end
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if (it < itmax) && (nres > tol)        
        tem =(abs(rm));
        wt = tem/norm(tem,1);
        beta = normw(xm,wt); %%% -- show that weighted techniques is used
        V(:,1) = xm/beta;
     end
end
iter = it;
res = res(1:iter + 1);
end
