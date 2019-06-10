function [idxb,idxi]=get_bdidx3d(Nx,Ny,Nz,mode)
% a simple script to get the boundary (and internal index) index for
% edges, faces and nodes for vectorized matrix for 3D FD grid
% set here counters for Edges, Faces and Nodes of the grid


% �������Ӧ�ú�������X��Y��Z���������й�ϵ
switch upper(mode)
    case 'EDGE'  %���Ԫ�ǰ�����ߵķ������У�����x�򣨼�������x������صı߽����⣩������y���������z����ÿһ�����������ڲ�ѭ��Ϊ�У���һ�����꣩�����ѭ��Ϊ�У��ڶ������꣩��������Ǹߣ����������꣩
        x = ones(Nx,Ny+1,Nz+1); %��x����ƽ�е���߾���
        y = ones(Nx+1,Ny,Nz+1); %��y����ƽ�е���߾���
        z = ones(Nx+1,Ny+1,Nz); %��z����ƽ�е���߾���
        % set the boundaries to be zeros...
        % upper and lower boundary:
        x(:,:,[1,Nz+1])=0;      %���ַ����ܼ�㣬�ܺ�
        y(:,:,[1,Nz+1])=0;
        % left and right boundary:
        x(:,[1,Ny+1],:)=0;
        z(:,[1,Ny+1],:)=0;
        % front and back boundary:
        y([1, Nx+1],:,:)=0;
        z([1, Nx+1],:,:)=0;
        t = [reshape(x,Nx*(Ny+1)*(Nz+1),1); ...
             reshape(y,(Nx+1)*Ny*(Nz+1),1); ...
             reshape(z,(Nx+1)*(Ny+1)*Nz,1)];
    case 'FACE'  %��Ԫ�ǰ�����ķ��������У�����x�򣨼�������x������صı߽����⣩������y���������z����ÿһ�����������ڲ�ѭ��Ϊ�У���һ�����꣩�����ѭ��Ϊ�У��ڶ������꣩��������Ǹߣ����������꣩
        x = ones(Nx+1,Ny,Nz);
        y = ones(Nx,Ny+1,Nz);
        z = ones(Nx,Ny,Nz+1);
        % set the boundaries to be zeros...
        % upper and lower boundary:
        z(:,:,[1,Nz+1])=0;
        % left and right boundary:
        y(:,[1,Ny+1],:)=0;
        % front and back boundary:
        x([1, Nx+1],:,:)=0;
        t = [reshape(x,(Nx+1)*Ny*Nz,1); ...
             reshape(y,Nx*(Ny+1)*Nz,1); ...
             reshape(z,Nx*Ny*(Nz+1),1)];
    case 'NODE'  %�ڵ�Ԫ�ǰ��սڵ㣬�ڸ���ͬ��������û�з���ĸ�����ڲ�ѭ��Ϊ�У���һ�����꣩�����ѭ��Ϊ�У��ڶ������꣩��������Ǹߣ����������꣩
        x = ones(Nx+1,Ny+1,Nz+1);
        % set the boundaries to be zeros...
        % upper and lower boundary:
        x(:,:,[1,Nz+1])=0;
        % left and right boundary:
        x(:,[1,Ny+1],:)=0;
        % front and back boundary:
        x([1, Nx+1],:,:)=0;
        t = reshape(x,(Nx+1)*(Ny+1)*(Nz+1),1);
    otherwise
        error('mode not recognized, must be EDGE, FACE or NODE')
end
% interior edges/faces/nodes
idxi = find(t==1);
% boundary edges/faces/nodes
idxb = find(t==0);
return