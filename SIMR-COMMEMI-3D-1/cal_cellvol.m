function Vcell=cal_cellvol(Dx,Dy,Dz,mode)
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
% mode could be 'CELL' or 'NODE'
if nargin<3
    error('not enough input parameters')
elseif nargin<4
    mode='CELL';
end
Nx=size(Dx,1);
Ny=size(Dy,1);
Nz=size(Dz,1);
switch upper(mode)
    case 'CELL'
        Ncell=Nx*Ny*Nz;
        [Ly,Lx,Lz]=meshgrid(Dy,Dx,Dz);
        V=Ly.*Lx.*Lz;
    case 'NODE'
        Ncell=(Nx+1)*(Ny+1)*(Nz+1);
        DDx=([Dx(1:end);0]+[0;Dx(1:end)])/2;
        DDy=([Dy(1:end);0]+[0;Dy(1:end)])/2;
        DDz=([Dz(1:end);0]+[0;Dz(1:end)])/2;        
        [Ly,Lx,Lz]=meshgrid(DDy,DDx,DDz);
        V=Ly.*Lx.*Lz;
    otherwise
        error('mode not recognized, must be CELL or NODE.')
end
Vcell = reshape(V,Ncell,1);
return