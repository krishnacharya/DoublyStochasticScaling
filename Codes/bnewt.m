function [x, res, MVP] = bnewt(A,tol,x0,delta,fl,Delta, outerIterLmt)
% BNEWT A balancing algorithm for symmetric matrices
%
% X = BNEWT(A) attempts to find a vector X such that
% diag(X)*A*diag(X) is close to doubly stochastic. A must
% be symmetric and nonnegative.
%
% X = BNEWT(A,TOL) specifies the tolerance of the method.
% If TOL = [] then BNEWT uses the default, 1e-6.
%
% X = BNEWT(A,TOL,X0) specifies the initial guess. If none
% is given then BNEWT uses the default, vectors of ones.
%
% X = BNEWT(A,TOL,X0,DEL) determines how close we are
% willing to allow our balancing vectors can get to the edge of the
% positive cone. We use a relative measure on the size of elements. The
% default value for DEL is 0.1
%
% X = BNEWT(A,TOL,X0,DEL,FL) will output intermediate convergence
% statistics if FL = 1. The default value for FL is 1.
%
%
% [X, RES] = BNEWT(A,TOL,X0,DEL,FL) returns the residual error, too.

% BNEWT is based on a Newton step to solve
% diag(X)*A*diag(X) - 1 = 0 as an outer iteration. This can be written as
% x_new = (A + diag((A*x_old)./(xold)))\(A*x_old + 1./x_old).
% We solve this step approximately by using CG as an iteration and with
% preconditioner D = diag(xold) applied to both sides of the system.

% The iteration continues until the residual is smaller than TOL.
% The residual is measured by norm(diag(x)*A*x - e,2).

% Initialise
[n, ~]=size(A); e = ones(n,1); res=[];
if nargin < 7, outerIterLmt = 10; end

if nargin < 6, Delta = 3; end
if nargin < 5,  fl = 0; end
if nargin < 4,  delta= 0.1; end
if nargin < 3,  x0 = e; end
if nargin < 2,  tol = 1e-6; end

g=0.9; etamax = 0.1; % Parameters used in determining inner iteration stopping criterion. Default g=0.9; eta = 0.1;
eta = etamax; stop_tol = tol*.5;  %Relative tolerance?

x = x0; rt = tol^2; 
v = x.*(A*x); rk = e - v; 
rho_km1 = rk'*rk; rout = rho_km1; rold = rout;    

MVP = 0;  % We'll count matrix vector products.
i = 0; % Outer iteration count.
 
if fl == 1, fprintf('it    in. it    res\n'), end
while rout > rt  % Outer iteration
    i = i + 1; k = 0; y = e; 
    innertol = max([eta^2*rout,rt]);
    while  rho_km1 > innertol %Inner iteration by CG
        k = k + 1;
        if k == 1
            Z = rk./v;  p=Z; rho_km1 = rk'*Z; 
        else
            beta=rho_km1/rho_km2;
            p=Z + beta*p;
        end        
        % w is not just A*p as we are solving a system with matrix
        % B = diag(x)*A*diag(x)+diag(diag(x)*A*diag(x)e). 
        % This can be done efficiently without forming B.
        w = x.*(A*(x.*p)) + v.*p;
        alpha = rho_km1/(p'*w);
        ap = alpha*p;
        % We want to avoid -ve components. To achieve this we limit how
        % much closer we can get to the boundary in each inner iteration.
        ynew = y + ap;
        if min(ynew) <= delta
            if delta == 0, break, end
            ind = find(ap < 0);
            gamma = min((delta  - y(ind))./ap(ind));
            y = y + gamma*ap;
            break
        end  
        if max(ynew) >= Delta
            ind = find(ynew > Delta);
            gamma = min((Delta-y(ind))./ap(ind));
            y = y + gamma*ap;
            break
        end 
        y = ynew; 
        rk = rk - alpha*w;  rho_km2 = rho_km1;
        Z = rk./v;  rho_km1 = rk'*Z;
    end
     x = x.*y; 
    v = x.*(A*x); 
    rk = 1 - v; rho_km1 = rk'*rk; rout = rho_km1;
    
    MVP = MVP + k + 1;
    
    % Update inner iteration stopping criterion.
    rat = rout/rold;  rold = rout; res_norm = sqrt(rout);
    eta_o = eta;  eta = g*rat;
    if g*eta_o^2 > 0.1
        eta = max([eta,g*eta_o^2]);
    end
    eta = min([eta,etamax]);
    eta = max([eta,stop_tol/ res_norm]);
    
    if fl == 1
        fprintf('%3d %6d   %.3e %.3e %.3e %.3e\n', i,k, res_norm,min(y),max(y),min(x));
        res=[res; res_norm];
    end
    if (i >= outerIterLmt)
        fprintf('outer iterlimit (%d) attained, retiurning\n', outerIterLmt);
        break;
    end
end

 fprintf('Matrix-vector products = %6d\n', MVP)