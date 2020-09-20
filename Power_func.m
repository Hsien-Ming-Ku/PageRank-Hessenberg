%%%% Power method for PageRank computation, Ax = x;
% written by Dr. Z.-L. Shen (SICAU) and Dr. C. Wen (UESTC)
% Modified by Dr. Xian-Ming Gu (SWUFE)
% The related paper: 
% convergence criterion: ||Ax-x ||_1/||x||_1 < tol, 
% where A = alpha*(P + v*d') + (1 - alpha)*v*e';
function [x0, iter, mv, res] = Power_func(A,q,tol,maxit)
%--------------- check the input A ------------------------------
        if isa(A, 'double')
           matrix_or_function = 1;
        elseif isa(A, 'function_handle')
           matrix_or_function = 2;
        else
          error('A is neither a matrix or a function handle!!!');
        end % end if
%%%%----------------------------------------------------------------
%%% --------------- begin the power iteration --------------------------
x0 = q; 
%------------------------------ initial process ------------------------
if matrix_or_function == 1
   xappro = A*(x0);
else
   xappro = A(x0);
end
% in fact:    xappro = alpha * (A' * x0 + (d' * x0) * v) + (1 - alpha) * v;
%-------------------- 
iter = 0; 
mv = 1;
res = zeros(maxit+1,1);
relres = norm(xappro - x0,1)/norm(x0,1);
res(1) = relres;
x0 = xappro;
while (iter < maxit) && (relres > tol)
%-----------------------------------------------------------------------
      if matrix_or_function == 1
         xappro = A*(x0);
      else
         xappro = A(x0);
      end
%     xappro = alpha * (A' * x0 + (d' * x0) * v) + (1 - alpha) * v;
%-------------------- 
      mv = mv + 1;
      iter = iter + 1; 
      relres = norm(xappro - x0,1)/norm(x0,1);
      res(iter + 1) = relres;
      if relres < tol
         break;
      end        
      x0 = xappro;
end
res = res(1:iter+1);
end