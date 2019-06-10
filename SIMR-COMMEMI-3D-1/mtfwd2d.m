function [Ux,Vy,Vz]=mtfwd2d(Dy,Dz,res,per,pol) %, Tol, Niter)
% a simple driver function for MT 2d forward script using finite difference 
% method... to be used for 3D boundary conditions
% NOTE this is a silent version, which does not generate any annoying
% messages 
% this calculate Nper periods and *one* polarization
%
% function mtfwd1d is called to generate the boundary condition
% 
% essentially this solves the Maxwell equations based on E/B
% for E, 
%   curl(curl(E))+i*omega*mju*sigma*E=0
% as we realy just have one E component (x), the equation can be simplified
% to be:
%   d(d(Ex)/dy)/dy + d(d(Ex)/dz)/dz - i*omega*mju*sigma*Ex = 0
% for B, similarly
%   curl(rho*curl(B))+i*omega*mju*B=0
% which can also be simplified: 
%   d(rho*d(Bx)/dy)/dy + d(rho*d(Bx)/dz)/dz - i*omega*mju*Bx = 0
% *this* could be helpful
%       |  i  j  k  |
%       | dx dy dz  |
%       | Ux Uy Uz  |
% DONG Hao 2011.06.30 in Qumarleb
%=========================================================================%
% input and output parametres:
% Dy:       Ny array of column width of each cell column of the model
% Dz:       Nz array of layer depth of each cell layer of the model 
% res:      Ny by Nz matrix of resistivity of the model 
% per:     array of calculated data period 
% pol:      polarization mode of planar wave (can be 'E' or 'B')
%           in 'E'(TE), Ux=Ex, Vy=Bx, Vz=Bz;
%           in 'B'(TM), Ux=Bx, Vy=Ex, Vz=Ez;
%
% the following is currently not used as we now use a direct solver
%           
% Tol:      target error tolerance for iteration
% Niter:    maximum number of iteration (for each linear solver)
%=========================================================================%
% other parametres that might be useful:
% Ny:       number of model columns (cells)
% Nz:       number of model layers (cells)
% Nper:     number of periods
% Nza:       number of Air layers 
% Ra:       resistivity of Air layers
% sigma:    Nx by Ny by Nz cell conductivity, 1/res
%=========================================================================%
if nargin<4
    error('not enough input parameters')
elseif nargin<5
    pol='E';
end
% switch upper(pol)
%     case 'E'
%         disp('calculating E Polarization')
%     case 'B'
%         disp('calculating B Polarization')
%     otherwise
%         error('polarization mode must be E or B')
% end

%=========================================================================%
%                          setting up some parameters
%=========================================================================%
mju0=4.0*pi*1E-7; % unit is Wb/(A��m) in SI. 
Ra=1E10; % Air resisitivity
isign=-1; % -i*omega*mju
Nper=length(per);
Ny=length(Dy);
Nz=length(Dz);
if size(res,1)~=Ny||size(res,2)~=Nz
    error('size of resistivity not matching sizes of column and layer')
end
Nza=find(res(1,:)<Ra,1)-1;
sigma=1./res;
% average edge length 
% i.e. distances of mid points of adjacent cell, used for 2nd order FD
DDy=([Dy(1:end);0]+[0;Dy(1:end)])/2;
DDz=([Dz(1:end);0]+[0;Dz(1:end)])/2;

% allocating memory
% please note that U should be located on nodes and V should be located on
% cell edges. 
Ux=zeros(Ny+1,Nz+1,Nper);
Vy=Ux;
Vz=Ux;
%=========================================================================%
%                    Setting up the matrix operators
%=========================================================================%
% compute number of edges and nodes
Nnode = (Ny+1)*(Nz+1);
NedgeY = Ny*(Nz+1);
NedgeZ = (Ny+1)*Nz;
Nedge = NedgeY+NedgeZ;
% create the differential operator
G=mk_grad2d(Ny,Nz);
% now build up the cell area vector (will be used as diag matrix later)
[tdz,tdy]=meshgrid(Dz,Dy);
Acell = tdy.*tdz;
sta = sigma.*Acell; % sigma-times-area 
rta = 1./sigma.*Acell; % res-times-area
% ...the average area vector of (inner) nodes
[tdz,tdy]=meshgrid(DDz,DDy);
AAnode = reshape(tdy.*tdz,(Ny+1)*(Nz+1),1);
% ...the edge length vector of horizontal and vertical edges
Ledge = reshape(Dy*ones(1,Nz+1),NedgeY,1);% y
Ledge = [Ledge;reshape(ones(Ny+1,1)*Dz',NedgeZ,1)];% z
% and the average edge length vector of horizontal and vertical edges
ALedge = reshape(ones(Ny,1)*DDz',NedgeY,1); 
ALedge = [ALedge;reshape(DDy*ones(1,Nz),NedgeZ,1)];
% get the boundary and inner field (at each node) indexing...
[idxb,idxi]=get_bdidx2d(Ny+1,Nz+1);
%=========================================================================%
%                       Now loop through the periods
%=========================================================================%
for iper=1:Nper
    omega=2*pi./per(iper);
%     disp(['calculating period ' num2str(per(iper))])
    switch upper(pol)
        case 'E'
%=========================================================================%
%                         building the operators
%=========================================================================%
            A =  spdiags(1./AAnode,0,Nnode,Nnode)*G'...
                *spdiags(ALedge./Ledge,0,Nedge,Nedge)*G;
            % start calculating averaged sigma over the nodes
            % 
            sigma_avg = zeros(Ny+1,Nz+1);
            sigma_avg(2:end-1,2:end-1)=(sta(1:end-1,1:end-1)+...
                sta(2:end,1:end-1)+sta(1:end-1,2:end)+sta(2:end,2:end))/4;
            sigma_avg = sigma_avg./reshape(AAnode,Ny+1,Nz+1);
            %   i*omega*mu*sigma
            iwus = isign*1j*omega*mju0*reshape(sigma_avg,(Ny+1)*(Nz+1),1);
            A = A + spdiags(iwus,0,Nnode,Nnode);
%=========================================================================%
%            setup source - currently a uniform source is used 
%     Note to self: an interface should be reserved for future (CSEM) use
%=========================================================================%
% changed to use 1D boundary here
            Ux(:,:,iper)=setbc2dmt(Ux(:,:,iper),Dy,Dz,per(iper),res,pol); 
            x = reshape(Ux(:,:,iper),Nnode,1);
            b = -A(idxi,idxb)*x(idxb);
            A = A(idxi,idxi);
        case 'B'
%=========================================================================%
%                         building the operators
%=========================================================================%
            % start calculating averaged sigma over the Y/Z edges
            % respectively 
            % averaged between two adjacent cells along Y edge
            res_avgy = zeros(Ny,Nz+1);
            res_avgy(:,2:end-1) = rta(:,1:end-1)+rta(:,2:end);
            % no, you don't really want to average the air with the earth
            % at the air-earth boundary ...
            res_avgy(:,Nza+1) = rta(:,Nza+1)*2;
            res_avgy(:,2:end-1) = res_avgy(:,2:end-1)./(Acell(:,1:end-1)+ ...
                        Acell(:,2:end)); %Y
            % averaged between two adjacent cells along Z edge
            res_avgz = zeros(Ny+1,Nz);
            res_avgz(2:end-1,:) = rta(1:end-1,:)+rta(2:end,:);
            res_avgz(2:end-1,:) = res_avgz(2:end-1,:)./(Acell(1:end-1,:)+ ...
                        Acell(2:end,:)); %Z
            res_avg = [reshape(res_avgy,NedgeY,1);reshape(res_avgz,NedgeZ,1)];
            % 
            A =  spdiags(1./reshape(AAnode,Nnode,1),0,Nnode,Nnode)*G'...
                *spdiags(res_avg.*ALedge./Ledge,0,Nedge,Nedge)*G;
            %   i*omega*mu
            iwu = ones(Ny+1,Nz+1); 
            iwu(2:end-1,2:end-1) = isign*1j*omega*mju0;
            A = A + spdiags(reshape(iwu,(Ny+1)*(Nz+1),1),0,Nnode,Nnode);
%=========================================================================%
%            setup source - currently a uniform source is used 
%     Note to self: an interface should be reserved for future (CSEM) use
%=========================================================================%
% changed to use 1D boundary here
            Ux(:,:,iper)=setbc1dmt(Ux(:,:,iper),Dy,Dz,per(iper),res,pol); 
            x = reshape(Ux(:,:,iper),Nnode,1);
            b = -A(idxi,idxb)*x(idxb);
            A = A(idxi,idxi);
        otherwise
            error('polarization option not recognized, must be "E" or "B"')
    end
%=========================================================================%
%                     now solve the linear system
%=========================================================================%
%   direct solve with a simple (jacobian) pre-conditioner ...
    T = spdiags(1./spdiags(A,0),0,length(idxi),length(idxi)); % T=P^(-1)
    A = T*A;
    b = T*b;
    x(idxi) = A\b; 
    Ux(:,:,iper) = reshape(x,Ny+1,Nz+1);
%=========================================================================%
%                          now E to B/ B to E
%=========================================================================%
    switch upper(pol)
        case 'E'
            [Vy(:,:,iper),Vz(:,:,iper)] = E2B(Ux(:,:,iper),Dy,Dz,G,...
                Ledge,isign,omega);
        case 'B'
            [Vy(:,:,iper),Vz(:,:,iper)] = B2E(Ux(:,:,iper),Dy,Dz,G,...
                Ledge,res_avg);
    end
end
return