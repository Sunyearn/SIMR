function [x, isfailed, residual, Niter, resvec]=solver(A,b,x,method,MaxIter,Tol,...
                                        Nedge)
% a simple driver function to solve 'almost' any (sparse matrix based)
% linear system, provided by A, b and x
% 
% available methods for solving the problem include 'direct', 'bicg',
% 'qmr' and 'gmres'
if nargin<5
    MaxIter=200;
    Tol=1e-10;
    Nedge=[];
elseif nargin <6
    Tol=1e-10;
    Nedge=[];
elseif nargin <7
    Nedge=[];
end
isfailed=0;
switch lower(method)
    case 'direct' 
%   direct solver here - just use the matlab default '\' 
%   USE AT YOUR OWN RISK
%   could take years to compute (if you do not run out of memory first)
%   preconditioning with a simple (jacobian) pre-conditioner ...
        n = size(A,1);
        T = spdiags(1./spdiags(A,0),0,n,n); % T1=P^(-1)
        A1 = T*A;
        b1 = T*b;
        x = A1\b1;
        residual = norm(A*x- b)/norm(b);
        Niter = 1;
        resvec = residual;
    case 'qmr'
%   QMR here - uses less memory and faster for larger system
%   slower than bicg, but seems to yield more reasonable response
        if ~isempty(Nedge) % we are solving the system function
            [L,U] = precond(A,'BILU',Nedge);
        else
            [L,U] = precond(A,'SOR');
        end
        [x, isfailed, residual, Niter, resvec] = qmr(A,b,Tol,MaxIter,L,U,x);
    case 'bicg'
%   BiCG here - uses less memory and much faster for larger system
%   seems to be the most efficient so far
%   probably will be used in future version of ModEM 
%   the bad thing is the 'unconverged' result seems worse than QMR
%   even if the relative residual is at the same level - strange
        if ~isempty(Nedge) % we are solving the system function
            [L,U] = precond(A,'BILU',Nedge);
        else
            [L,U] = precond(A,'ILU');
        end
%        [x isfailed residual Niter] = bcgs(A,b,Tol,MaxIter,L,U,x);
        [x, isfailed, residual, Niter, resvec] = bicgstab(A,b,Tol,MaxIter,L,U,x);
    case 'bcgs'
%   home-made BiCGstab - have not tested that (yet)
        if ~isempty(Nedge) % we are solving the system function
            [L,U] = precond(A,'BILU',Nedge);
        else
            [L,U] = precond(A,'ILU');
        end
    
        [x, isfailed, residual, Niter, resvec] = bcgs(A,b,Tol,MaxIter,L,U,x);
    case 'gmres'
%   GMRES here - uses fewer (outer) iterations
%   and hence less divergence corrections
%   takes about twice the time as of bicg
        if ~isempty(Nedge) % we are solving the system function
            [L,U] = precond(A,'BILU',Nedge);
        else
            [L,U] = precond(A,'SSOR');
        end
        [x, isfailed, residual, iter, resvec] = gmres(A,b,5,Tol,MaxIter,L,U,x);
        Niter=iter(1);
    case 'pcg'
%   PCG here - seems to be not stable (at all!)
%       The only good part is, well, the result looks somewhat reasonable
%       even if it does not converge well.
%       maybe I got something wrong...
%       but this has never been completely tested. 
%       if you have time reading this, test it!
        if ~isempty(Nedge) % we are solving the system function
            [L,U] = precond(A,'BILU',Nedge);
        else
            [L,U] = precond(A,'IC');
        end
        [x, isfailed, residual, Niter, resvec] = pcg(A,b,Tol,MaxIter,L,U);
    otherwise
        error('solver method not recognized, must be direct, qmr, bicg, pcg or gmres')
end
return