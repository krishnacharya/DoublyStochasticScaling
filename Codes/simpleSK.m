function [c,r] = simpleSK(A,maxIters,tol)
iter = 0;
[m, n] = size(A)
r = ones(m,1)
while(iter < maxIters)
    disp("iteration number "+iter)
    c = 1./( A'*r );
    r = 1./( A *c );
    iter = iter + 1;
    disp(r')
    disp(c')
    disp(diag(r)*A*diag(c))
end