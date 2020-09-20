function y = normw(x,wt)
% To compute the weighted norm of x with repect to 
% the weighting vector wt, this subroutine will be used in 
% the weighted Arnoldi process;
% Written by Dr. Zhao-Li Shen (SCIAU)
y = sqrt(x'*(x.*wt));
end