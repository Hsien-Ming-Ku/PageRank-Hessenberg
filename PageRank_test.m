% Numerical tests for computing PageRank by means of 
% 1. Arnoldi-type method
% 4. Power method;
% 5. Power-LE method; from Dr. Xueyuan Tan (JCAM 2017)
% ---------------------------------------------------------------
% written by Dr. Xian-Ming Gu, who recently works at the School of Economic
% Mathematics, SWUFE.
% Date: 24 July, 2018, the Visting Scholar at University of Macau
%       7 January, 2019 (modified some parts), at SWUFE
%       24 August, 2020 (finalized the numerical results) at SWUFE.
% Email: guxianming@live.cn, guxm@swufe.edu.cn
%
%% --------------- Setting some parameters -----------------------
clear;
% clc;
close;
tol = 1e-8;
dim = 8;
maxit = 800;
%% ------------ generate the PageRank problems ------------
load amazon-2008; % the web adjacency matrix 
% cnr-2000 fail
% in-2000 fail
% eu-2005
%load('D:\sx-stackoverflow.mat');  % web adjacency matrix file
G = Problem.A;
% G = CS;
clear Problem;
%%%        initialization of the PageRank problem 
%%%        fomating A into A = P' + h*f
alpha = 0.99;                    % Damping factor for pagerank problem
n = size(G,2);                   % dimension of the web adjacency matrix G
nnzG = nnz(G);
I = speye(n);
v = 1/n*ones(n,1);               % the personalization vector
sumrow = full(sum(G,2));
dangling = find(sumrow==0);      % set of dangling nodes
l = length(dangling);            % amount of dangling nodes
mu = (1-alpha) + alpha*l/n;  
d = zeros(1,n);
d(dangling) = 1;                 % the index vector representing dangling nodes 
sumrow(dangling) = 1;
P = diag(sparse(1./sumrow))*G;   % the transition matrix of the web link graph
clear subrow dangling G;
% --------------------------------
P = alpha*P;                           
h = v;            
f = alpha*d + (1 - alpha)*ones(1,n);
P = P';
Atrace = trace(P) + sum(f.*h');
A = @(x)mv_pagerank(P,h,f,x);
%%%%%%%%%%%%%%%%%%%
q = ones(n,1);
q = q/norm(q,1);
% matvec = 0;
%%---- Test the Arnoldi-type algorithm --------------------
tic;
[q1,iter1,res1] = ArnoldiPagerank(A,q,dim,maxit,tol);
t1 = toc;
fprintf('- CPU time for Arnoldi-type algorithm: %.4f.\n ', t1);
fprintf('- Outer iters for Arnoldi-type algorithm: %.4g.\n ', iter1);
% %% --------------------------------------------------------------------
% ---- Test the weighted Arnoldi-type algorithm --------------------
tic;
[q2,iter2,res2] = weighted_ArnoldiPagerank(A,q,tol,dim,maxit);
t2 = toc;
fprintf('- CPU time for GArnoldi-type algorithm: %.4f.\n ', t2);
fprintf('- Outer iters for GArnoldi-type algorithm: %.4g.\n ', iter2);
 % ---- Test the Power method -----------------------
 maxit = dim*maxit;
 tic;
 [q3,iter3,mv1,res3] = Power_func(A,q,tol,maxit);
 t3 = toc;
 fprintf('- CPU time for the Power-type method: %.4f.\n ', t3);
 fprintf('- Number of MatVecs for the Power-type method: %.4g.\n ', iter3);
 %%% --------------------------------------------------------------------
 % ----- Test the Power + extrapolation method (X. Tan, 2017) --------
 maxit = dim*maxit;
 tic;
 [q4,iter4,mv2,res4] = Power_extrapolation_func(A,q,alpha,l,tol,maxit);
 t4 = toc;
 fprintf('- CPU time for the Power-Ext method: %.4f.\n ', t4);
 fprintf('- Number of MatVecs for the Power-Ext method: %.4g.\n ', iter4);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% ---- Test the Power + QE algorithm --------------------
 maxit = dim*maxit;
 tic;
 [x0, iter, mv, res] = Power_QE_func(A, q, tol, maxit);
 t6 = toc;
 fprintf('- CPU time for Power + QE method: %.4f.\n ', t6);
 fprintf('- Outer iters for Power + QE method: %.4g.\n ', iter);
 %%%%%%%%%%% Test the Hessenberg algorithm ------------------
tic;
[q5,iter5,res5] = HessenPagerank(A,q,dim,maxit,tol);
t5 = toc;
fprintf('- CPU time for Hessenberg-type algorithm: %.4f.\n ', t5);
fprintf('- Outer iters for Hessenberg-type algorithm: %.4g.\n ', iter5);

disp('+++++++++++++++++++++++ end test ++++++++++++++++++++++++++++++++');
