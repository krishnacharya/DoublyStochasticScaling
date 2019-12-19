function [c, r, iter, err] = knight_sk(B, maxIters, tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [c, r, iter, err] = knight_sk(B, maxIters, tol)
%
% implements P A Knights's presentation of Sinkhorn-Knopp algorithm
%      for scaline A in 1-norm
% 
%Upon exit B = sparse(1:size(A, 1), 1:size(A, 1), r) * A * sparse(1:size(A, 1), 1:size(A, 1), c)
%      is the soubly stochastic matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by  Daniel Ruiz and Bora Ucar, 
% references: A symmetry preserving algorithm for matrix scaling
%             Daniel Ruiz and Bora Ucar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ii,jj,vv] = find(B);
[m, n] = size(B);

A = sparse(ii,jj, abs(vv), m, n);

if(nnz(B) ~= nnz(A)),
    error('lost some nonzeros');
end
r = ones(m,1);
res = 1;
iter = 0;
while ((res > tol) && (iter<maxIters)),
    c = 1./( A'*r );
    r = 1./( A *c );
    Dr = spdiags(r,0,m,m);
    Dc = spdiags(c,0,n,n);
    res = max(norm(Dr*A*c-1,inf), norm(Dc*A'*r-1,inf));
    iter = iter + 1;
    err = res;
    if(mod(iter, 100) == 0),
        fprintf('iter %d err %1.6e\n', iter, res);
    end
end

%err = res;
