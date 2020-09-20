% Power extrapolation method for PageRank computation, Ax = x
% written by Z.-L. Shen (SAU)
% related paper: Dr. Xueyuan Tan (JCAM 2017), extrapolation techniques
% X. Tan, A new extrapolation method for PageRank computations, J. Comput. 
% Appl. Math., Vol. 313, 2017, pp. 383-392.
% ---------------------------------------------------------------
% Convergence criterion: ||Ax - x||_1/||x||_1 < tol, 
% where A = alpha*(P + v*d') + (1 - alpha)*v*e';
function [x,iter,mv,res] = Power_extrapolation_func(A,q,alpha,l,tol,maxit)
%--------------- check the input A ------------------------------
        if isa(A, 'double')
           matrix_or_function=1;
        elseif isa(A, 'function_handle')
           matrix_or_function=2;
        else
           error('A is neither a matrix or a function handle!!!');
        end % end if
%%%%----------------------------------------------------------------
% set parameters for power-extrapolation method
n = length(q);
m = 50;   % the window size suggested by Dr. Tan in JCAM's paper
mu = (1 - alpha) + alpha*l/n;                          % refer to the paper [1]
%  run the power-extrapolation method
x = q;                                                                  % initial guess
% normr=norm(P'*x+(f*x)*h-x,1);              % initial resudial norm
if matrix_or_function == 1
   qt = x - A*x;
else
   qt = x - A(x);
end
normr = norm(qt,1)/norm(x,1);
res = zeros(maxit+1,1);
res(1) = normr;
iter = 0;                                                                   % the number of iterations
mv = 1;                                                               % the number of mat-vec operations
while(normr>tol && iter<maxit)
%      x_tem=P'*x;
%      x_tem=x_tem+(f*x)*h;
%      x_tem=A(x);
     if matrix_or_function == 1
       x_tem = A*x;
     else
       x_tem = A(x);
     end
     normr = norm(x_tem - x,1)/norm(x,1);
     res(iter + 1) = normr;
     if normr < tol
         break;
     end
     if (mod(iter,m)==0)
        x_pre = x_tem; 
     end
     if (mod(iter,m)==1 && (iter~=1))   
   % ---begin extrapolation process ---------------
        x_tem = x_tem - (mu - 1)*x_pre;
        x_tem = x_tem/norm(x_tem,1);
     end
     x = x_tem;
     iter = iter+1;    
     mv = mv + 1;
end
res = res(1:iter + 1);               %resudial 
end