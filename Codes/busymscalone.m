function [DDr, DDc, As, numIters, Err] = busymscalone( A, NbIter, butol, Dr, Dc )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% usage  [DDr, DDc, As, numIters, Err] = busymscalone( A, NbIter, butol, Dr, Dc )
%        Dr and Dc are optional inputs.
% Scales the matrix 1./Dr*A*1/Dc with 1./DDr and 1./DDc to have
% all rows and cols in the scaled matrix As to have 1-norm equal to 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by  Daniel Ruiz and Bora Ucar,
% references: A symmetry preserving algorithm for matrix scaling
%             Daniel Ruiz and Bora Ucar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numIters = 0;
if(nargin<3),
    error('give at least three args');
end
[m,n] = size(A);
if (nargin == 3)
    Dr = ones(1,m);  Dc = ones(1,n);
end
DDr = Dr(:)';  DDc = Dc(:)';  As = A;

d = max(abs(A'));  I = find(d==0); I1 = ones(size(I));
d = max(abs(A));   J = find(d==0); J1 = ones(size(J));

if ( sum([size(I1) size(J1)]) ~= 0 ),
    if ( issparse(A) )
        % If A (input) is sparse, we treat everything accordingly
        for k = 1:NbIter;
            As = spdiags(ones(m,1)./Dr(:), 0, m, m) * As * spdiags(ones(n,1)./Dc(:), 0, n, n);
            Dr = full(sqrt(sum(abs(As'))));  Dc = full(sqrt(sum(abs(As))));
            Dr(I) = I1;  Dc(J) = J1; %see this handles zero rows cols
            DDr = DDr.*Dr; DDc = DDc.*Dc;
            
            ErrR = max(abs(Dr.^2-ones(1,m)));
            ErrC = max(abs(Dc.^2-ones(1,n)));
            Err = max(ErrR, ErrC);
            if(mod(k,10000) == 1),
             %   fprintf('iter =%d err = %1.6f\n', k, Err);
            end
            if(Err<=butol),
                break;
            end
        end
        As = spdiags(ones(m,1)./Dr(:), 0, m, m) * As * spdiags(ones(n,1)./Dc(:), 0, n, n);
        %  End case when A is sparse
    else
        % If A (input) is full, things are more simple ...
        for k = 1:NbIter;
            As = diag(ones(1,m)./Dr) * As * diag(ones(1,n)./Dc);
            Dr = sqrt(sum(abs(As')));  Dc = sqrt(sum(abs(As)));
            Dr(I) = I1;  Dc(J) = J1;
            DDr = DDr.*Dr; DDc = DDc.*Dc;
            
            ErrR = max(abs(Dr.^2-ones(1,m)));
            ErrC = max(abs(Dc.^2-ones(1,n)));
            Err = max(ErrR, ErrC);
            if(Err<=butol),
                break;
            end
        end
        As = diag(ones(1,m)./Dr) * As * diag(ones(1,n)./Dc);
        %  End case when A is full
    end
else
    if ( issparse(A) )
        % If A (input) is sparse, we treat everything accordingly
        for k = 1:NbIter;
            As = spdiags(ones(m,1)./Dr(:), 0, m, m) * As * spdiags(ones(n,1)./Dc(:), 0, n, n);
            Dr = full(sqrt(sum(abs(As'))));  Dc = full(sqrt(sum(abs(As))));
            DDr = DDr.*Dr; DDc = DDc.*Dc;
            ErrR = max(abs(Dr.^2-ones(1,m)));
            ErrC = max(abs(Dc.^2-ones(1,n)));
            Err = max(ErrR, ErrC);
            if(mod(k,10000) == 1),
                fprintf('iter =%d err = %1.6f\n', k, Err);
            end
            if(Err<=butol),
                break;
            end
        end
        As = spdiags(ones(m,1)./Dr(:), 0, m, m) * As * spdiags(ones(n,1)./Dc(:), 0, n, n);
        %  End case when A is sparse
    else
        % If A (input) is full, things are more simple ...
        error('Full case is not fully covered');
        for k = 1:NbIter;
            As = diag(ones(1,m)./Dr) * As * diag(ones(1,n)./Dc);
            Dr = sqrt(sum(abs(As')));  Dc = sqrt(sum(abs(As)));
            DDr = DDr.*Dr; DDc = DDc.*Dc;
            ErrR = max(abs(Dr.^2-ones(1,m)));
            ErrC = max(abs(Dc.^2-ones(1,n)));
            Err = max(ErrR, ErrC);
            if(Err<=butol),
                break;
            end
        end
        As = diag(ones(1,m)./Dr) * As * diag(ones(1,n)./Dc);
        %  End case when A is full
    end
end
numIters = k;

rsums = abs(abs(As) * ones(size(As,2),1)-1);
csums = abs(abs(As') * ones(size(As,1),1)-1);


Err = max(max(rsums), max(csums));

