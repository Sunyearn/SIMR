function [irow,icol,ivalue]=mk_grad3d(Nx,Ny,Nz,Nelem,mode)
%  a simple function to create the 'topology' part of 3D gradient operator 
%  (matrix) for Finite difference method
%  basicly this just build the 'sign' part of 
%
%  grad(U) = ex(dUx/dx)+ey(dUy/dy)+ez(dUz/dz)
%          = ex((Ux(i+1)-Ux(i))/Dx(i))
%          + ey((Uy(j+1)-Uy(j))/Dy(j))
%          + ez((Uz(k+1)-Uz(k))/Dz(k))
%  which is essentially differencial operators along each axis 
%          1    2    3
%        3 *----*----* 1
%         /    /    /|
%      2 *----*----* |
%    Ex ^    /    /| * 2
%    1 *-->-*----* |/|
%      | Ey |    | * |
%   Ez V    |    |/| * 3
%      *----*----* |/
%      |    |    | *
%      |    |    |/
%      *----*----*
%  
%  modified to use sparse matrix (row, column, value format) to deal with
%  larger (not so large, actually) problems - this was MUCH easier to setup
%  when using full matrix...
%  Nelem could be 3 by 1 Nface or Nedge 
%  mode could be 'N2E' or 'C2F'   N2E 和 C2F都什么模式呢？答：node to edge 和 cell to (inner) face 

%  DONG Hao 2012.07.11 in Beijing
%=========================================================================%
if nargin<4
    error('not enough input parameters')
elseif nargin<5
    mode='N2E';
end
switch upper(mode)
   case 'N2E' % node to edge
        Nedge = Nelem;
        Nelem2 = (Nx+1)*(Ny+1)*(Nz+1); % Nnode here
		% for x-edges (2*Nedgex edges in total)
		Nedgex = Nedge(1);
        % determine rows
		irowx = reshape(ones(2,1)*(1:Nedgex),Nedgex*2,1);%irowx的意思是row Ux的序号i？？irowx拿来做什么用，与后面的icolx有关联么？
		% determine columns 
        % just a little trick to convert from 3d to 1d index
        I = mod(1:Nedgex,Nx);
        I(I==0) = Nx;   % 我去，这操作犀利！！！！！！！！！！！！！！！！！！！！！！
        J  = mod(ceil((1:Nedgex)/(Nx)),(Ny+1));
        J(J==0) = Ny+1;
        K = ceil((1:Nedgex)/((Nx)*(Ny+1)));
		icolx =reshape( [(K-1)*(Nx+1)*(Ny+1)+(J-1)*(Nx+1)+I ; ... %为什么上一行提示的K是((Nx)*(Ny+1))的倍数取整，这里却是(Nx+1)*(Ny+1)的倍数
                        (K-1)*(Nx+1)*(Ny+1)+(J-1)*(Nx+1)+I+1 ],Nedgex*2,1); %???????????????????????????????
		% determine the topology - 
		valuex = reshape([-ones(1,Nedgex); ones(1,Nedgex)] , Nedgex*2,1);
        
		% for y--edges (2*Nedgey edges in total)
		Nedgey = Nedge(2);
		% determine rows
		irowy = reshape(ones(2,1)*(Nedgex+1:Nedgex+Nedgey),Nedgey*2,1);
		% determine columns 
        I = mod(1:Nedgey,Nx+1);
        I(I==0) = Nx+1;
        J  = mod(ceil((1:Nedgey)/(Nx+1)),(Ny));
        J(J==0) = Ny;
        K = ceil((1:Nedgey)/((Nx+1)*(Ny)));
		icoly =reshape([(K-1)*(Nx+1)*(Ny+1)+(J-1)*(Nx+1)+I; ... %?????????????????????????????????
                        (K-1)*(Nx+1)*(Ny+1)+(J)*(Nx+1)+I], Nedgey*2,1); %???????????????????????????????
		% determine the topology - 
		valuey = reshape([-ones(1,Nedgey); ones(1,Nedgey)],Nedgey*2,1);
		
		% for z-edges (2*Nedgez edges in total)
		Nedgez = Nedge(3);
        % determine rows
		irowz = reshape(ones(2,1)*(Nedgex+Nedgey+1:Nedgex+Nedgey+Nedgez),Nedgez*2,1);
		% determine columns 
		I = mod(1:Nedgez,Nx+1);
        I(I==0) = Nx+1;
        J  = mod(ceil((1:Nedgez)/(Nx+1)),(Ny+1));
        J(J==0) = Ny+1;
        K = ceil((1:Nedgez)/((Nx+1)*(Ny+1)));
		icolz =reshape([(K-1)*(Nx+1)*(Ny+1)+(J-1)*(Nx+1)+I;... %??????????????????????????????????
                        (K)*(Nx+1)*(Ny+1)+(J-1)*(Nx+1)+I], Nedgez*2,1); %????????????????????????????????
		% determine the topology + -  
		valuez = reshape([-ones(1,Nedgez); ones(1,Nedgez)],Nedgez*2,1);
	case 'C2F' % cell to (inner) face 
        % this is a little bit complicated as the 1d index for inner faces
        % is no longer continuous 
        % 
		% for x inner faces ((Nx-1)*Ny*Nz faces) 
        % (number of all x faces is (Nx+1)*Ny*Nz)
        Nface=Nelem;
        Nelem2 = Nx*Ny*Nz; % Ncell here
        ifacex = ones(Nx+1,Ny,Nz);
        % mask out the boundary faces...
        ifacex([1 Nx+1],:,:)=0;
        ifacex=find(ifacex==1);
		Nfacex = length(ifacex);
        % determine rows
		irowx = reshape(ones(2,1)*ifacex',Nfacex*2,1);
		% determine columns
        I = mod(ifacex,Nx+1);
        I(I==0) = Nx+1;
        J  = mod(ceil((ifacex)/(Nx+1)),Ny);
        J(J==0) = Ny;
        K = ceil((ifacex)/((Nx+1)*Ny));
		icolx = reshape([((K-1)*(Nx)*(Ny)+(J-1)*(Nx)+I-1)' ;...
                         ((K-1)*(Nx)*(Ny)+(J-1)*(Nx)+I)'], Nfacex*2,1); 
		% determine the topology - 
		valuex = reshape([-ones(1,Nfacex); ones(1,Nfacex)],Nfacex*2,1);
		
		% for y inner faces (Nx*(Ny-1)*Nz faces) 
        % (number of all y faces is Nx*(Ny+1)*Nz)
        ifacey = ones(Nx,Ny+1,Nz);
        % mask out the boundary faces...
        ifacey(:,[1, Ny+1],:)=0;
        ifacey=find(ifacey==1);
		Nfacey = length(ifacey);% y-faces
        % determine rows
		irowy = reshape(ones(2,1)*(Nface(1)+ifacey)',Nfacey*2,1);
		% determine columns
        I = mod(ifacey,Nx);
        I(I==0) = Nx;
        J  = mod(ceil((ifacey)/(Nx)),Ny+1);
        J(J==0) = Ny+1;
        K = ceil((ifacey)/((Nx)*(Ny+1)));
		icoly =reshape([((K-1)*(Nx)*(Ny)+(J-2)*(Nx)+I)'; ...
                        ((K-1)*(Nx)*(Ny)+(J-1)*(Nx)+I)'], Nfacey*2,1); 
		% determine the topology - 
		valuey = reshape([-ones(1,Nfacey); ones(1,Nfacey)],Nfacey*2,1);
		
		% for z inner faces (Nx*Ny*(Nz-1) faces) 
        % (number of all z faces is Nx*Ny*(Nz+1))
        ifacez = ones(Nx,Ny,Nz+1);
        % mask out the boundary faces...
        ifacez(:,:,[1 Nz+1])=0;
        ifacez=find(ifacez==1);
		Nfacez = length(ifacez);
        % determine rows
		irowz = reshape(ones(2,1)*(Nface(1)+Nface(2)+ifacez)',Nfacez*2,1);
		% determine columns
        I = mod(ifacez,Nx);
        I(I==0) = Nx;
        J  = mod(ceil((ifacez)/Nx),Ny);
        J(J==0) = Ny;
        K = ceil((ifacez)/(Nx*Ny));
		icolz =reshape([((K-2)*(Nx)*(Ny)+(J-1)*(Nx)+I)';...
                        ((K-1)*(Nx)*(Ny)+(J-1)*(Nx)+I)'], Nfacez*2,1); 
		% determine the topology + - 
		valuez = reshape([-ones(1,Nfacez); ones(1,Nfacez)],Nfacez*2,1);    
    otherwise
        error('mode not recognized, must be N2E or C2F')
end 
%  now construct sparse matrix for all three directions
%  G should be Nedge by Nnode for NTE
%              Nface by Ncell for CTF
G = sparse([irowx;irowy;irowz],[icolx;icoly;icolz],...
    [valuex;valuey;valuez],sum(Nelem),Nelem2);
irow = [irowx;irowy;irowz];
icol = [icolx;icoly;icolz];
ivalue = [valuex;valuey;valuez];
return