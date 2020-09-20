function w = QE(v1, v2, v3, v,n)
% The subroutine is to implement the Quadratic Extrapolation for 
% accelerating vector iterations. Then it will be used to accelerate the 
% Power iteration;
% Written by Dr. Guojian Yin (HKUST) 
% and modified by Dr. Xian-Ming Gu (SWUFE)
% Reference:
% S.D. Kamvar, T.H. Haveliwala, C.D. Manning, G.H. Golub, Extrapolations 
% methods for accelerating PageRank computations, WWW'2003, May 20-24, 2003, 
% Budapest, Hungary.
% ---------------------------------------------------------------------
% v1, v2, v3, v is the recent iteration vector u^{k-3}, u^{k-2}, u^{k-1},
% u^{k}
% ------------------------------------------------
V3 = [v2,v3,v];
V4 = zeros(n,3);
% Y = zeros(n,2);
% gamma = zeros(2,1);
for k = 1:3
    V4(:,k) = V3(:,k) - v1;
end
Y = [V4(:,1),V4(:,2)];
gamma_3 = 1;
gamma = -pinv(Y)*V4(:,3);
% size(Y)
% gamma = (-Y)\V4(:,3);
gamma_1 = gamma(1);
gamma_2 = gamma(2);
gamma_0 = -(gamma_1 + gamma_2 + gamma_3);
beta_0 = -gamma_0;
beta_1 = gamma_2 + gamma_3;
beta_2 = gamma_3;
w = beta_0*V3(:,1) + beta_1*V3(:,2) + beta_2*V3(:,3);