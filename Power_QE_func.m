% Power method accelerated by the quadratic extrapolation 
% for PageRank computations: A*x = x;
% written by Dr. Xian-Ming Gu (SWUFE), 
%% related paper: 
%% convergence criterion ||Ax-x ||_1<tol, A=alpha*(P+v*d')+(1-alpha)*v*e';

%% if the input all_res==1, the residual for every single iteration is recorded
%% else only the residual after every restart cycle is recorded 
%% output: it: the number of restarted cycles
function [x0, mv,iter, res] = Power_QE_func(A, q, tol, maxit)
%--------------- check the input A ------------------------------
if isa(A, 'double')
   matrix_or_function=1;
elseif isa(A, 'function_handle')
   matrix_or_function=2;
else
  error('A is neither a matrix or a function handle!!!');
end % end if
%%%%----------------------------------------------------------------

%%%%%%%%%%%%% begin power iteration %%%%%%%%%%%%%%%%
m = 50;  % the window size: using one QE per m steps of Power iteration.
x0 = q; 
if matrix_or_function==1
   xappro=A*(x0);
else
   xappro= A(x0);
end
%     xappro = alpha * (A' * x0 + (d' * x0) * v) + (1 - alpha) * v;
%-------------------- 
n = length(x0);
iter = 0; 
mv = 1;
res = zeros(maxit+1,1);
relres = norm(xappro - x0,1)/norm(x0,1);
res(1) = relres; 
% ---------------------
V1(:,1) = x0;
if matrix_or_function == 1
   v  = A*(V1(:,1));
else
   v  = A(V1(:,1));
end
V1(:,2) = v/norm(v);
if matrix_or_function == 1
   v = A*V1(:,2);
else
   v = A(V1(:,2));
end
% v  = A(V1(:,2));
V1(:,3) = v/norm(v);
while (iter < maxit) && (relres > tol)
%-------replace the original mv_pr used in PWFOM project by the new one,
%for myself
    %w=mv_pr(temV(j).v);
     if matrix_or_function == 1
        w = A*V1(:,3);
     else
        w = A(V1(:,3));
     end
%     xappro = alpha * (A' * x0 + (d' * x0) * v) + (1 - alpha) * v;
%-------------------- 
      mv = mv + 1;
      iter = iter + 1;
%     lamda = V1(:,3)'*w;
      lambda = 1;
      r = w - lambda*V1(:,3);
      relres = norm(r,1);
      xappro = w/norm(w);
      norm1_x = norm(xappro,1);
%       relres = norm(w - xappro,1)/norm1_x;
      res(iter+1) = relres/norm1_x;
      if relres < tol*norm1_x
         x0 = xappro/norm(xappro,1);
         break;
      end
      if (mod(iter,m))==0 
    % 'quadratic extrapolation'
         w = QE(V1(:,1), V1(:,2), V1(:,3),xappro,n);
         xappro = w/norm(w);
      %'end quadratic extrapolation'
      end
      V1(:,1) = V1(:,2);
      V1(:,2) = V1(:,3);
      V1(:,3) = xappro;
end        
res = res(1:iter+1);
end