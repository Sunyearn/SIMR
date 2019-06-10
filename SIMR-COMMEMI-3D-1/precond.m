function [L,U]=precond(A,method,Nedge)
% a simple driver function to switch between different preconditioners here
% this is supposed to be a generized function (but they never are, for my
% functions)
% 
% available methods include 'NONE', 'JAC', 'BILU', 'SOR', 'SSOR' and 'ILU'
% 
% it should be noted that the matlab iLU method would take like forever to
% make the L/U. 
switch upper(method)
    case 'NONE'
        % as the name implies. probably won't converge using this, but
        % isn't it a good practice?
        L = spdiags(ones(size(A,1),1),0,size(A,1),size(A,2));
        U = L;
    case 'JAC' 
        % Jacobian. This could make smaller systems solveable, but may be
        % as slow as direct method.
        U = spdiags(spdiags(A,0),0,size(A,1),size(A,1));
        L = spdiags(ones(size(A,1),1),0,size(A,1),size(A,2));
    case 'SOR' % DLU
        % successive over-relaxation, simple and suprisingly efficient
        % comparing with Jacobian
        L = tril(A);
        U = (diag(diag(A)))\triu(A);
    case 'SSOR' % DLLT
        % symmetric successive over-relaxation, not (overly) complex
        n = size(A,1);
        % w = 2/(1+sin(pi/sqrt(size(A,1)))); % w=1.8
        w = 1.75;
        D = spdiags(diag(A)./w,0,n,n);
        L = 1/(2-w)*(D + tril(A,-1));
        U = D\(D + triu(A,1));
    case 'IC'
        % Incomplete Cholesky factorization 
        % could complain about nonpositive pivot.
        % L = ichol(A, struct('type','nofill','michol','on'));
        L = ichol(A, struct('type','nofill','michol','off'));
%         L = ichol(A);
        U = L';
    case 'ILU'
        % Incomplete LU factorization
        [L, U] =ilu(A, struct('type','ilutp','droptol',1e-2));
    case 'BILU' % block ILU - x, y, z are setup for different blocks
        % seems to be the most efficient so far
        % probably will be used in future version of ModEM 
        Nedgex = Nedge(1);
        Nedgey = Nedge(2);
        Nedgez = Nedge(3);
        ix = 1:Nedgex;
        iy = Nedgex+1:Nedgex+Nedgey;
        iz = Nedgex+Nedgey+1:Nedgex+Nedgey+Nedgez;

%         [Lx,Ux] = ilu(A(ix,ix),struct('type','ilutp','droptol',1e-1,'udiag',0));
%         [Ly,Uy] = ilu(A(iy,iy),struct('type','ilutp','droptol',1e-1,'udiag',0));
%         [Lz,Uz] = ilu(A(iz,iz),struct('type','ilutp','droptol',1e-1,'udiag',0));
        [Lx,Ux] = ilu(A(ix,ix),struct('type','nofill'));
        [Ly,Uy] = ilu(A(iy,iy),struct('type','nofill'));
        [Lz,Uz] = ilu(A(iz,iz),struct('type','nofill'));
        L=blkdiag(Lx,Ly,Lz);
        U=blkdiag(Ux,Uy,Uz);
    case 'MILU' % ILU with more blocks...
        nblock= 6;
        edges=sum(Nedge);
        nlen = floor(edges/nblock);
        L=[];
        U=[];
        for iblock=1:nblock
            ix = (iblock-1)*nlen+1:nlen*iblock;            
            if iblock == nblock
                ix = ix(1):edges;
            end
            [Li,Ui] = ilu(A(ix,ix),struct('type','nofill'));
            L=blkdiag(L,Li);
            U=blkdiag(U,Ui);
        end
    otherwise 
        error('method not recognized, must be SOR, SSOR, IC, ILU or BILU')
end