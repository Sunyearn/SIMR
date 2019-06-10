function sigma_avg=ave_sigma(sigma,Vcell,Nelem,mode)
% a simple script to average sigma onto the edges or the faces of the mesh
% used for building the i*omega*mju*sigma vector 
% Nelem could be 3 by 1 Nface or Nedge 
% mode could be 'EDGE' or 'FACE'
if nargin<3
    error('not enough input parameters')
elseif nargin<4
    mode='EDGE';
end
[Nx,Ny,Nz]=size(sigma);
Vcell=reshape(Vcell,Nx,Ny,Nz);
stv = sigma.*Vcell;
switch upper(mode)
    case 'EDGE' 
        % sigma averaged on x edges
        stv_temp = zeros(Nx,Ny+2,Nz+2);
        stv_temp(:,2:(end-1),2:(end-1)) = stv;
        stv_temp(:,1,2:(end-1)) = stv_temp(:,2,2:(end-1));
        stv_temp(:,end,2:(end-1)) = stv_temp(:,end-1,2:(end-1));
        stv_temp(:,:,1) = stv_temp(:,:,2);
        stv_temp(:,:,end) = stv_temp(:,:,end-1);
        Vcell_temp = zeros(Nx,Ny+2,Nz+2);
        Vcell_temp(:,2:(end-1),2:(end-1)) = Vcell;
        Vcell_temp(:,1,2:(end-1)) = Vcell_temp(:,2,2:(end-1));
        Vcell_temp(:,end,2:(end-1)) = Vcell_temp(:,end-1,2:(end-1));
        Vcell_temp(:,:,1) = Vcell_temp(:,:,2);
        Vcell_temp(:,:,end) = Vcell_temp(:,:,end-1);
        sx=zeros(Nx,Ny+1,Nz+1);
        sx(:,:,:) = stv_temp(:,1:end-1,1:end-1)+stv_temp(:,2:end,1:end-1)+...
            stv_temp(:,1:end-1,2:end)+stv_temp(:,2:end,2:end);
        sx(:,:,:) = sx(:,:,:)./(Vcell_temp(:,1:end-1,1:end-1)+ ...
            Vcell_temp(:,2:end,1:end-1)+Vcell_temp(:,1:end-1,2:end)+Vcell_temp(:,2:end,2:end));
        % sigma averaged on y edges
        stv_temp = zeros(Nx+2,Ny,Nz+2);
        stv_temp(2:(end-1),:,2:(end-1)) = stv;
        stv_temp(1,:,2:(end-1)) = stv_temp(2,:,2:(end-1));
        stv_temp(end,:,2:(end-1)) = stv_temp(end-1,:,2:(end-1));
        stv_temp(:,:,1) = stv_temp(:,:,2);
        stv_temp(:,:,end) = stv_temp(:,:,end-1);
        Vcell_temp = zeros(Nx+2,Ny,Nz+2);
        Vcell_temp(2:(end-1),:,2:(end-1)) = Vcell;
        Vcell_temp(1,:,2:(end-1)) = Vcell_temp(2,:,2:(end-1));
        Vcell_temp(end,:,2:(end-1)) = Vcell_temp(end-1,:,2:(end-1));
        Vcell_temp(:,:,1) = Vcell_temp(:,:,2);
        Vcell_temp(:,:,end) = Vcell_temp(:,:,end-1);
        sy=zeros(Nx+1,Ny,Nz+1);
        sy(:,:,:) = stv_temp(1:end-1,:,1:end-1)+stv_temp(2:end,:,1:end-1)+...
            stv_temp(1:end-1,:,2:end)+stv_temp(2:end,:,2:end);
        sy(:,:,:) = sy(:,:,:)./(Vcell_temp(1:end-1,:,1:end-1)+ ...
            Vcell_temp(2:end,:,1:end-1)+Vcell_temp(1:end-1,:,2:end)+Vcell_temp(2:end,:,2:end));
        % sigma averaged on z edges
        stv_temp = zeros(Nx+2,Ny+2,Nz);
        stv_temp(2:(end-1),2:(end-1),:) = stv;
        stv_temp(1,2:(end-1),:) = stv_temp(2,2:(end-1),:);
        stv_temp(end,2:(end-1),:) = stv_temp(end-1,2:(end-1),:);
        stv_temp(:,1,:) = stv_temp(:,2,:);
        stv_temp(:,end,:) = stv_temp(:,end-1,:);
        Vcell_temp = zeros(Nx+2,Ny+2,Nz);
        Vcell_temp(2:(end-1),2:(end-1),:) = Vcell;
        Vcell_temp(1,2:(end-1),:) = Vcell_temp(2,2:(end-1),:);
        Vcell_temp(end,2:(end-1),:) = Vcell_temp(end-1,2:(end-1),:);
        Vcell_temp(:,1,:) = Vcell_temp(:,2,:);
        Vcell_temp(:,end,:) = Vcell_temp(:,end-1,:);
        sz=zeros(Nx+1,Ny+1,Nz);
        sz(:,:,:) = stv_temp(1:end-1,1:end-1,:)+stv_temp(2:end,1:end-1,:)+...
            stv_temp(1:end-1,2:end,:)+stv_temp(2:end,2:end,:);
        sz(:,:,:) = sz(:,:,:)./(Vcell_temp(1:end-1,1:end-1,:)+ ...
            Vcell_temp(2:end,1:end-1,:)+Vcell_temp(1:end-1,2:end,:)+Vcell_temp(2:end,2:end,:));
    case 'FACE' % not completely tested yet, but probably will work...
        % sigma averaged on x faces
        sx=zeros(Nx+1,Ny,Nz);
        sx(2:end-1,:,:) = stv(1:end-1,:,:)+stv(2:end,:,:);
        sx(2:end-1,:,:) = sx(2:end-1,:,:)...
            ./(Vcell(1:end-1,:,:)+ Vcell(2:end,:,:));
        % sigma averaged on y faces
        sy=zeros(Nx,Ny+1,Nz);
        sy(2:end-1,:,:) = stv(:,1:end-1,:)+stv(:,2:end,:);
        sy(2:end-1,:,:) = sy(:,2:end-1,:)...
            ./(Vcell(:,1:end-1,:)+ Vcell(:,2:end,:));
        % sigma averaged on z faces
        sz=zeros(Nx,Ny,Nz+1);
        sz(2:end-1,:,:) = stv(:,:,1:end-1)+stv(:,:,2:end);
        sx(2:end-1,:,:) = sx(2:end-1,:,:)...
            ./(Vcell(:,:,1:end-1)+ Vcell(:,:,2:end));
    otherwise
        error('mode not recognized, must be FACE or EDGE')
end
% assemble and convert sigma into a vector to build 
% i*omega*mju*sigma
sigma_avg=[reshape(sx,Nelem(1),1);...
           reshape(sy,Nelem(2),1);...
           reshape(sz,Nelem(3),1)];
return