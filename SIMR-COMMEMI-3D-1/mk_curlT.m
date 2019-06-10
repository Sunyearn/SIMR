function T=mk_curlT(Nx,Ny,Nz,Nedge,Nface)
%  a simple function to create the "topology" part of 3D curl operator
%  matrix (hence called "T") for FD method 
%  Both the curl from edge to face, and the curl from (averaged) 
%  face to (averaged) edge involves the same 'topology' part of T 
%                           1         2
%       *                    edge y(+)
%       |\                  *---->----*   1
%    Ez v ^ Ex              |         |
%       |  \                |         |
% Hy <----  *---  edge z(-) ^    x    v edge z(+)
%       *   |               |  face x |
%        \  |               |         |
%     -Ex v ^ -Ez           *----<----*   2
%          \|                edge y(-)
%           *
%  basicly this just build the 'sign' part of 
%  curl(U) = ex(dUz/dy-dUy/dz)+ey(dUx/dz-dUz/dx)+ez(dUy/dx-dUx/dy)
%          = ex((Uz(k+1)-Uz(k))/Dy(j)-(Uy(j+1)-Uy(j))/Dz(k)) 
%          + ey((Ux(i+1)-Ux(i))/Dz(k)-(Uz(k+1)-Uz(k))/Dx(i))
%          + ez((Uy(j+1)-Uy(j))/Dx(i)-(Ux(i+1)-Ux(i))/Dy(j))
%  
%  essentially, a face has four edges around it - 
%  take face x for example, there will be two edge ys and two edge zs.
%  the x component of curl operator will involve 
%
%   ex((Uz(k+1)-Uz(k))/Dy(j)-(Uy(j+1)-Uy(j))/Dz(k)) 
% 
%  which happens to reflect the two "plus" edges (Uz(k+1),Uy(k)) that along
%  the "integration path" and two "minus" edges (Uz(k),Uy(k+1))that against 
%  the path.
%  
%  modified to use sparse matrix (row, column, value format) to deal with
%  larger (not so large, actually) problems
%
%  DONG Hao 2012.07.11 in Beijing
%=========================================================================%
% first make the indices for edges - 
% x
Nedgex=Nedge(1);
iedgex = zeros(Nx,Ny+1,Nz+1);
iedgex(1:Nedgex) = 1:Nedgex;
% y
Nedgey=Nedge(2);
iedgey = zeros(Nx+1,Ny,Nz+1);
iedgey(1:Nedgey) =Nedgex+1:Nedgex+Nedgey;
% z
Nedgez=Nedge(3);
iedgez = zeros(Nx+1,Ny+1,Nz);
iedgez(1:Nedgez) = Nedgex+Nedgey+1:Nedgex+Nedgey+Nedgez;

% for x-faces, which only involves y and z edges (4*Nfacex edges in total)
Nfacex = Nface(1); %number of face -> (Nx+1)*Ny*Nz
% determine rows
irowx = reshape(ones(4,1)*(1:Nfacex),Nfacex*4,1);
% determine columns 
I = mod(1:Nfacex,Nx+1);
I(I==0) = Nx+1;
J  = mod(ceil((1:Nfacex)/(Nx+1)),Ny);
J(J==0) = Ny;
K = ceil((1:Nfacex)/((Nx+1)*Ny));
icolx = reshape( ...
    [ iedgey((K-1)*(Nx+1)*Ny+(J-1)*(Nx+1)+I); ... % iedgey - Nx+1 by Ny
      iedgez((K-1)*(Nx+1)*(Ny+1)+(J-1)*(Nx+1)+I); ... % iedgez - Nx+1 by Ny
      iedgey((K)*(Nx+1)*Ny+(J-1)*(Nx+1)+I); ... % iedgey (j+1)
      iedgez((K-1)*(Nx+1)*(Ny+1)+(J)*(Nx+1)+I)],... % iedgez (k+1)
      Nfacex*4,1);
% determine the topology - 
value = [ones(1,Nfacex) ; -ones(1,Nfacex) ; -ones(1,Nfacex) ; ones(1,Nfacex)];
valuex = reshape(value,4*Nfacex,1);

%  for y-faces, which only involves x and z edges (4*Nfacey edges in total)
Nfacey = Nface(2); %number of face -> Nx*(Ny+1)*Nz
irowy = reshape(ones(4,1)*(Nfacex+1:Nfacex+Nfacey),Nfacey*4,1);
% determine columns 
I = mod(1:Nfacey,Nx);
I(I==0) = Nx;
J  = mod(ceil((1:Nfacey)/Nx),Ny+1);
J(J==0) = Ny+1;
K = ceil((1:Nfacey)/(Nx*(Ny+1)));
icoly = reshape( ...
    [ iedgex((K-1)*Nx*(Ny+1)+(J-1)*Nx+I); ...           % iedgex - Nx by Ny+1
      iedgez((K-1)*(Nx+1)*(Ny+1)+(J-1)*(Nx+1)+I); ...   % iedgez - Nx by Ny+1
      iedgex((K)*Nx*(Ny+1)+(J-1)*Nx+I); ...             % iedgex (j+1)
      iedgez((K-1)*(Nx+1)*(Ny+1)+(J-1)*(Nx+1)+I+1);], Nfacey*4,1);  % iedgez (i+1)
% determine the topology -   
value = [-ones(1,Nfacey) ; ones(1,Nfacey) ; ones(1,Nfacey) ; -ones(1,Nfacey)];
valuey = reshape(value,4*Nfacey,1);

%  for z-faces, which only involves x and y edges (4*Nfacez edges in total)
Nfacez = Nface(3); %number of face -> Nx*Ny*(Nz+1)
irowz=reshape(ones(4,1)*(Nfacex+Nfacey+1:Nfacex+Nfacey+Nfacez),Nfacez*4,1);
% determine columns 
I = mod(1:Nfacez,Nx);
I(I==0) = Nx;
J  = mod(ceil((1:Nfacez)/Nx),Ny);
J(J==0) = Ny;
K = ceil((1:Nfacez)/(Nx*Ny));
icolz = reshape( ...
    [ iedgex((K-1)*Nx*(Ny+1)+(J-1)*Nx+I); ...        % iedgex - Nx by Ny
      iedgey((K-1)*(Nx+1)*Ny+(J-1)*(Nx+1)+I); ...    % iedgey - Nx by Ny
      iedgex((K-1)*Nx*(Ny+1)+(J)*Nx+I); ...       % iedgex (j+1)
      iedgey((K-1)*(Nx+1)*Ny+(J-1)*(Nx+1)+I+1)], Nfacez*4,1); %iedgey (i+1)
% determine the topology -   
value = [ones(1,Nfacez) ; -ones(1,Nfacez) ; -ones(1,Nfacez) ; ones(1,Nfacez)];
valuez = reshape(value, Nfacez*4,1);
% now construct sparse matrix
T = sparse([irowx;irowy;irowz],[icolx;icoly;icolz],[valuex;valuey;valuez]);
% T here is Nface by Nedge matrix 
return