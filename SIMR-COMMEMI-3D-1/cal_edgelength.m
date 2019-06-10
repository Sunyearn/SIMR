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
% mode could be 'EDGE' or 'NODE' ����Ӧ����"EDGE"��"FACE"�ɣ�
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
        [Ly,Lx,Lz]=meshgrid(Dy,Dx,Dz);%֮���԰ѵ�һ�������Ӧdy������Ϊ�ڵ��þ���ʱ���ǣ�����ϰ���õ�һ������λ�ô����У��ڶ�������λ�ô�����
        % padding the extra rows,columns or layers %
        % ΪʲôҪ��������С��кͲ㡷����Ϊ����Ķ����С��С����Ӧ���������������Ԫ��Ӧ���е�����
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
        % subtracting the extra rows,columns or layers ɾȥ������С��кͲ�
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