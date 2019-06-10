function Aface=cal_facearea(Dx,Dy,Dz,mode)
% a simple function to build the face area vector for differential
% operators 
%  Ledge   |Lface|
% |-----|-----|-----|
% |     |     |     |
% |     |     |     |
% |-----|-----|-----|
% |     |     |     |
% |     |     |     |
% |-----|-----|-----|
% mode could be 'EDGE' or 'NODE'
if nargin<3
    error('not enough input parameters')
elseif nargin<4
    mode='FACE';
end
Nx=size(Dx,1);
Ny=size(Dy,1);
Nz=size(Dz,1);
switch upper(mode)
    case 'FACE'
        Nfacex=(Nx+1)*Ny*Nz;
        Nfacey=Nx*(Ny+1)*Nz;
        Nfacez=Nx*Ny*(Nz+1);                
        [Ly,Lx,Lz]=meshgrid(Dy,Dx,Dz);
        Ax=Ly.*Lz;
        Ay=Lx.*Lz;
        Az=Lx.*Ly;
        % padding the extra rows,columns or layers
        Ax(Nx+1,:,:)=Ax(Nx,:,:);
        Ay(:,Ny+1,:)=Ay(:,Ny,:);
        Az(:,:,Nz+1)=Az(:,:,Nz);
    case 'EDGE'
        Nfacex=(Nx)*(Ny+1)*(Nz+1);
        Nfacey=(Nx+1)*(Ny)*(Nz+1);
        Nfacez=(Nx+1)*(Ny+1)*(Nz);        
        DDx=([Dx(1:end);0]+[0;Dx(1:end)])/2;
        DDy=([Dy(1:end);0]+[0;Dy(1:end)])/2;
        DDz=([Dz(1:end);0]+[0;Dz(1:end)])/2;        
        [Ly,Lx,Lz]=meshgrid(DDy,DDx,DDz);
        Ax=Ly.*Lz;
        Ay=Lx.*Lz;
        Az=Lx.*Ly;
        % subtracting the extra rows,columns or layers
        Ax(Nx+1,:,:)=[];
        Ay(:,Ny+1,:)=[];
        Az(:,:,Nz+1)=[];       
    otherwise
        error('mode not recognized, must be EDGE or FACE.')
end
Aface = [reshape(Ax,Nfacex,1);...
         reshape(Ay,Nfacey,1);...
         reshape(Az,Nfacez,1)];
return