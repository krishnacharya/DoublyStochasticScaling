function [r, c, As, numIters, err ] = buScaleSK(A, maxIters, tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% usage  [r, c, As, numIters, err] = buScaleSK( A, maxIters, tol)
%        Dr and Dc are optional inputs.
% Scales the matrix 1./Dr*A*1/Dc with 1./DDr and 1./DDc to have 
% all rows and cols in the scaled matrix As to have 1-norm equal to 1
%
% Implements P A Knights description of Sinkhorn-Knopp algorithm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by  Daniel Ruiz and Bora Ucar, 
% references: A symmetry preserving algorithm for matrix scaling
%             Daniel Ruiz and Bora Ucar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[c, r, numIters, err] = knight_sk( A, maxIters, tol);

As = buScale(A, 1./r, 1./c);


