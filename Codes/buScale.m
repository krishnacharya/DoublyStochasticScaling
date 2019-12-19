function As = buScale(A, Dr, Dc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% As = buScale(A, Dr, Dc)
% Scales A with
%  As = spdiags(ones(size(D))./D, 0, n,n)* A *spdiags(ones(size(D))./D, 0, n, n);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by  Daniel Ruiz and Bora Ucar, 
% references: A symmetry preserving algorithm for matrix scaling
%             Daniel Ruiz and Bora Ucar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m,n] = size(A);
As = spdiags(ones(size(Dr))./Dr, 0, m, m) * A * spdiags(ones(size(Dc))./Dc, 0, n, n);

