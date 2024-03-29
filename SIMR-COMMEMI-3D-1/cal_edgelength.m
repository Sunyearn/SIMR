function Ledge=cal_edgelength(Dx,Dy,Dz,mode)
% a simple function to build the edge length vector for differential
% operators 
%  Ledge   |Lface|
% |-----|-----|-----|
% |     |     |     |
% |     |     |     |
% |-----|-----|-----|
% |     |     |     |
% |     |     |     |
% |-----|-----|-----|
% mode could be 'EDGE' or 'NODE' 这里应该是"EDGE"和"FACE"吧？
if nargin<3
    error('not enough input parameters')
elseif nargin<4
    mode='EDGE';
end
Nx=size(Dx,1);
Ny=size(Dy,1);
Nz=size(Dz,1);
switch upper(mode)
    case 'EDGE'
        Nedgex=(Nx)*(Ny+1)*(Nz+1);
        Nedgey=(Nx+1)*(Ny)*(Nz+1);
        Nedgez=(Nx+1)*(Ny+1)*(Nz);
        [Ly,Lx,Lz]=meshgrid(Dy,Dx,Dz);%之所以把第一个矩阵对应dy，是因为在调用矩阵时我们，我们习惯用第一个变量位置代表行，第二个变量位置代表列
        % padding the extra rows,columns or layers %
        % 为什么要额外加入行、列和层》？因为这里的额外行、列、层对应着三个方向上棱边元本应具有的量，
        Lx(:,Ny+1,:)=Lx(:,Ny,:); 
        Lx(:,:,Nz+1)=Lx(:,:,Nz);
        
        Ly(Nx+1,:,:)=Ly(Nx,:,:);
        Ly(:,:,Nz+1)=Ly(:,:,Nz);
        
        Lz(Nx+1,:,:)=Lz(Nx,:,:);
        Lz(:,Ny+1,:)=Lz(:,Ny,:);
    case 'FACE'
        Nedgex=(Nx+1)*Ny*Nz;
        Nedgey=Nx*(Ny+1)*Nz;
        Nedgez=Nx*Ny*(Nz+1);        
%         DDx=([Dx;0]+[0;Dx])/2;
%         DDy=([Dy;0]+[0;Dy])/2;  %%%  DONG
%         DDz=([Dz;0]+[0;Dz])/2;
        DDx=([Dx;0]+[0;Dx])/2;
        DDy=([Dy;0]+[0;Dy])/2;
        DDz=([Dz;0]+[0;Dz])/2;
        [Ly,Lx,Lz]=meshgrid(DDy,DDx,DDz);
        % subtracting the extra rows,columns or layers 删去额外的行、列和层
        Lx(:,Ny+1,:)=[];
        Lx(:,:,Nz+1)=[];
        Ly(Nx+1,:,:)=[];
        Ly(:,:,Nz+1)=[];
        Lz(Nx+1,:,:)=[];
        Lz(:,Ny+1,:)=[];
    otherwise
        error('mode not recognized, must be EDGE or FACE.')
end
Ledge = [reshape(Lx,Nedgex,1);...
         reshape(Ly,Nedgey,1);...
         reshape(Lz,Nedgez,1)];
return