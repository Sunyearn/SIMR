function [field3d]=mtfwd3d_initial(block,level,blo,per,method,MaxIter,Tol)
%=========================================================================%
% input and output parametres:
% Dx:       steps of grid in X direction
% Dy:       steps of grid in Y direction
% Dz:       steps of grid in Z direction
% res:      resistivity model
% per:      array of periods 
% pol:      polarization mode of plane wave         
% Tol:      target error tolerance for iteration
% Maxiter:  maximum number of iteration
%=========================================================================%
% other parametres
% Nx:       number of grid in X direction
% Ny:       number of grid in Y direction
% Nz:       number of grid in Z direction
% Nza:      number of Air layers 
% Nper:     the amount of periods
% Ra:       resistivity value of Air layers
% sigma:    conductivity model (=1./res)
%============================ prepare paras ==============================%
Dx = block(level).Dx(blo).vec;
Dy = block(level).Dy(blo).vec;
Dz = block(level).Dz(blo).vec;
res = block(level).res(blo).matrix;
%=========================================================================%
if nargin<5
    error('not enough input parameters')
elseif nargin<6
    method='bicg';
    MaxIter=400;
    Tol=1e-10;
elseif nargin<7
    MaxIter=400;
    Tol=1e-10;
elseif nargin<8
    Tol=1e-10;
end
%=========================================================================%
%                             some constants
%=========================================================================%
mju0=4.0*pi*1E-7; % unit is Wb/(A¡¤m) in SI. 
%=========================================================================%
%                      setting up some basic parameters
%=========================================================================%
polstr={'XY','YX'};   % polarization string
Nper=length(per);     % the amount of peroids
field3d = struct;     % data structure for output
Nx=length(Dx);
Ny=length(Dy);
Nz=length(Dz);

Nza=block(level).Nza(blo);% number of air layers of "blo" block at current "level" level
if size(res,1)~=Nx||size(res,2)~=Ny||size(res,3)~=Nz
    error('size of resistivity not matching sizes of column and layer')
end
sigma=1./res;   % conductivity model
Nnode=(Nx+1)*(Ny+1)*(Nz+1); % number of nodes
% number of edges in three directions
Nedgex=Nx*(Ny+1)*(Nz+1);    % the number of edges in X direction
Nedgey=(Nx+1)*Ny*(Nz+1);    % the number of edges in Y direction
Nedgez=(Nx+1)*(Ny+1)*Nz;    % the number of edges in Z direction
Nedge=[Nedgex,Nedgey,Nedgez];
Ne=sum(Nedge);              % sum of edges in three directions

Nfacex=(Nx+1)*Ny*Nz;        % the number of faces in X direction
Nfacey=Nx*(Ny+1)*Nz;        % the number of faces in Y direction
Nfacez=Nx*Ny*(Nz+1);        % the number of faces in Z direction
Nface=[Nfacex,Nfacey,Nfacez];
Nf=sum(Nface);              % sum of faces in three directions

% set up indices for boundary and inner elements
[idxeb,idxei]=get_bdidx3d(Nx,Ny,Nz,'edge'); % idxeb(index of boundary edges)£¬idxei(index of inside edges)
[idxnb,idxni]=get_bdidx3d(Nx,Ny,Nz,'node'); % idxnb(index of boundary nodes)£¬idxni(index of inside nodes)

% Nedgei=[sum of the inner edges in X direction, sum of the inner edges in Y direction,sum of the inner edges in Z direction]
Nedgei=[length(idxei(idxei<=Nedgex)),...
        length(idxei(idxei>Nedgex&idxei<=Nedgex+Nedgey)),...
        length(idxei(idxei>Nedgex+Nedgey))];

Ledge = cal_edgelength(Dx,Dy,Dz,'edge'); % array of length of edges
Lface = cal_edgelength(Dx,Dy,Dz,'face'); % array of average-length of edges at nodes
Aface = cal_facearea(Dx,Dy,Dz,'face'); % array of area of faces
Aedge = cal_facearea(Dx,Dy,Dz,'edge'); % array of average-area of faces at edges
Vcell = cal_cellvol(Dx,Dy,Dz,'cell'); % array of volumn of cells
sigma_avg=ave_sigma(sigma,Vcell,Nedge,'edge');% calculate the averaged sigma (on edges)

Ex=zeros(Nx,Ny+1,Nz+1); % electric field in X direction
Ey=zeros(Nx+1,Ny,Nz+1); % electric field in Y direction
Ez=zeros(Nx+1,Ny+1,Nz); % electric field in Z direction
eprev=zeros(Ne,2);  % a matrix contains the temperary electric field
%=========================================================================%
%                    now build up the matrix of operators
%========================== the grad operator  ===========================%
[irow,icol,ivalue]=mk_grad3d(Nx,Ny,Nz,Nedge,'N2E'); % [N2E: number of nodes] by [number of edges]
Tg = sparse(irow,icol,ivalue); % the 'topology' part of gradient operator
G=spdiags(1./Ledge,0,Ne,Ne)*Tg; % the gradient operator
%======================== the curl(-curl) operator  ======================%
Tc=mk_curlT(Nx,Ny,Nz,Nedge,Nface); % the 'topology' part of curl operator
C = spdiags(1./Aface,0,Nf,Nf)*Tc*spdiags(Ledge,0,Ne,Ne); % the curl operator
CC = spdiags(Ledge,0,Ne,Ne)*Tc'*spdiags(Lface./Aface,0,Nf,Nf)*Tc*...
     spdiags(Ledge,0,Ne,Ne); % CC = (Curl'*Volumn)*Curl, Volumn's scale is (amount of edges) by (amount of faces)
%============================ the div operator  ==========================%
Td=Tg';
D = Td*spdiags(Aedge,0,Ne,Ne);
Dts = D*spdiags(sigma_avg,0,Ne,Ne); % divergence operator
DtstG = Dts*G; % the matrix used to do divergence correction
%%
%=========================================================================%
%                       and loop through the periods
%=========================================================================%
for iper=1:Nper
    omega=2*pi./per(iper);
    disp(['computing period #' num2str(iper) ' @ ' num2str(per(iper)) 's'])
%=========================================================================%
%                    setting up frequency and sigma
%=========================================================================%
% take exp(iwu) as time harmonic dependence
    iwus = 1j * omega * mju0 * sigma_avg .* Ledge .* Aedge;
%===================== assemble the system equation  =====================%
    A = CC + spdiags(iwus,0,Ne,Ne);
    Dts(:,idxeb)=0;
    Dts(idxnb,:)=0;
%% 
%=========================================================================%
%                 setup source - right hand side or "b"
%     Note to self: an interface should be reserved for future (CSEM) use
%=========================================================================%
    for ipol=1:2
        disp(['for polarization ' num2str(ipol)])
        [Ex,Ey,Ez]=setbc3dmt(Ex,Ey,Ez,Dx,Dy,Dz,per(iper),res,'1D',...
            polstr{ipol},'E');% setup the boundary condition
        e=[reshape(Ex,Nedgex,1);...
           reshape(Ey,Nedgey,1);...
           reshape(Ez,Nedgez,1)];
        if iper~=1
           e(idxei)=eprev(idxei,ipol);
        end
        b = zeros(length(idxei),1)-A(idxei,idxeb)*e(idxeb); % construct the right hand side
%=========================================================================%
%                     now solve the linear system
%=========================================================================%
        isfailed=1;
        MaxsIter=40;
        MaxdIter=20;
        fcount=0;
        dIter=0;
        sIter=0;
        resall=zeros(MaxIter*2,1);
%%
        while sIter<MaxIter
%=========================================================================%
%                               solver
%=========================================================================%        
            [e(idxei), isfailed, residual, iter, resvec]=solver(A(idxei,idxei),...
                b,e(idxei),method,MaxsIter,Tol,Nedgei);
            if length(resvec)==1
                resall=resvec;
            else
                resall(sIter*2+1:sIter*2+iter*2)=resvec(1:iter*2);
            end
            disp([' A iteration# ', num2str(iter), ...
                  ' A current residual : ' num2str(residual)])
            phi = Dts(idxni,idxei)*e(idxei);
            nphi = norm(phi);
            disp([' A norm of Phi is: ', num2str(nphi)])
            sIter=sIter+iter;
            if isfailed==0
                disp(' A desired misfit reached, finishing...')
                break;
            elseif iter == 1
                if fcount==0 % give it a second chance...
                    disp(' A WARNING: THE SOLVER FAILED AFTER 1 ITER.')
                    fcount=1; 
                else % stop the iteration
                    disp(' A WARNING: THE SOLVER FAILED (AGAIN), aborting...')
                    disp(' A check your solver...')
                    break                    
                end
            end
%=========================================================================%
%                        Divergence correction 1
%=========================================================================%
            if strcmp(method,'direct')
                break % skip the divergence correction
            end
            e(idxei) = divcor(e(idxei),...
            Dts(idxni,idxei),DtstG(idxni,idxni),G(idxei,:),...
                idxni,idxnb);
            dIter=dIter+1;
            disp(['divergence correction # ' num2str(dIter)])
            if sIter>=MaxIter||dIter>=MaxdIter
                disp('Maximum iteration reached, exiting...')
                break
            elseif (isfailed == 0) 
                break
            end
        end
%%
        disp(['after totally ' num2str(sIter) ' iterations'])
        residual = A(idxei,idxei)*e(idxei) - b;
        residual = norm(residual) / norm(b);
        resall = resall/norm(b);
%         figure;
%         semilogy(0.5:0.5:sIter,resall(1:sIter*2));
        disp(['with a relative residual of ' num2str(residual)])
%=========================================================================%
%                   store the e for the next starting guess
%=========================================================================%
        eprev(:,ipol)=e;

%=========================================================================%
%                      now calculate H from E
%=========================================================================%
        h = (1j./mju0/omega)*C*e; % Ne E fields --> Nf Ffields

%=========================================================================%
%               save these fields in the storage structure
%=========================================================================%
        Hx=reshape(h(1:Nfacex),Nx+1,Ny,Nz);
        Hy=reshape(h(Nfacex+1:Nfacex+Nfacey),Nx,Ny+1,Nz);
        Hz=reshape(h(Nfacex+Nfacey+1:Nf),Nx,Ny,Nz+1);
        Ex=reshape(e(1:Nedgex),Nx,Ny+1,Nz+1);
        Ey=reshape(e(Nedgex+1:Nedgex+Nedgey),Nx+1,Ny,Nz+1);
        Ez=reshape(e(Nedgex+Nedgey+1:Ne),Nx+1,Ny+1,Nz);
        field3d(iper,ipol).Dx=Dx;
        field3d(iper,ipol).Dy=Dy;
        field3d(iper,ipol).Dz=Dz;
        field3d(iper,ipol).Ex=Ex;
        field3d(iper,ipol).Ey=Ey;
        field3d(iper,ipol).Ez=Ez;
        field3d(iper,ipol).Hx=Hx;
        field3d(iper,ipol).Hy=Hy;
        field3d(iper,ipol).Hz=Hz;        
        field3d(iper,ipol).Nza=Nza;
        field3d(iper,ipol).per=per(iper);
        field3d(iper,ipol).isfailed=isfailed;
    end
end
return 
