function [Ux,Uy,Uz]=setbc3dmt(Ux,Uy,Uz,Dx,Dy,Dz,per,res,mode,pol,type)
% setup boundary condition
Nx=size(Dx,1);
Ny=size(Dy,1);
Nz=size(Dz,1);
xi = [0; cumsum(Dx)];
yi = [0; cumsum(Dy)];
Ux(:)=0;
Uy(:)=0;
Uz(:)=0;
switch mode
    case '1D'
        switch pol
            case 'XY' % N-S source 
%=========================================================================%
%                         source on the top
%=========================================================================%
                Ux(:,:,1) = 1;  % can be either E or B field
                Uy(:,:,1) = 0;  % this is reduntant, actually
%=========================================================================%
%                         boundary on the left
%=========================================================================%
                % yes, run *MANY* 1d fwd here
                % loop through x edges...
                for i1d=1:Nx
                    res1d=reshape(res(i1d,1,:),Nz,1);
                    Ux(i1d,1,:)= mtfwd1d(per,res1d,Dz,type);
                end
                Uz(:,1,:)=0;
%=========================================================================%
%                         boundary on the right
%=========================================================================%
                % yes, run *MANY* 1d fwd here
                % loop through x edges...
                for i1d=1:Nx
                    res1d=reshape(res(i1d,Ny,:),Nz,1);
                    Ux(i1d,Ny+1,:)= mtfwd1d(per,res1d,Dz,type);
                end
                Uz(:,Ny+1,:)=0;
%=========================================================================%
%                         boundary on the bottom
%=========================================================================%
%   interpolate - you lier! you said you should just use the zeros!
                xi = (xi(1:Nx)+xi(2:Nx+1))/2;
                u = [Ux(:,1,end) Ux(:,end,end)];
                x = [xi xi];
                y = [ones(Nx,1)*yi(1) ones(Nx,1)*yi(end)];
                [yyi,xxi]=meshgrid(yi,xi);
                Ux(:,:,Nz+1) = interp2(y,x,u,yyi,xxi);
                Uy(:,:,Nz+1) = 0;
%=========================================================================%
%           now try to make a initial guess of the inner fields by
%                     interpolating the boundary
%=========================================================================%
% for now only do it for Ux 
                for iz=2:Nz
                    u = [Ux(:,1,iz) Ux(:,end,iz)];
                    Ux(:,:,iz) = interp2(y,x,u,yyi,xxi);
                end                
            case 'YX' % E-W source
%=========================================================================%
%                         source on the top
%=========================================================================%
                Uy(:,:,1) = 1;  % can be either E or B field
                Ux(:,:,1) = 0;  % this is reduntant, actually
%=========================================================================%
%                         boundary on the front
%=========================================================================%
                % yes, run *MANY* 1d fwd here
                % loop through y edges...
                for i1d=1:Ny
                    res1d=reshape(res(1,i1d,:),Nz,1);
                    Uy(1,i1d,:)= mtfwd1d(per,res1d,Dz,type);
                end
                Uz(1,:,:) = 0;  % this is reduntant, actually
%=========================================================================%
%                         boundary on the back
%=========================================================================%
                % yes, run *MANY* 1d fwd here
                % loop through y edges...
                for i1d=1:Ny
                    res1d=reshape(res(Nx,i1d,:),Nz,1);
                    Uy(Nx+1,i1d,:)= mtfwd1d(per,res1d,Dz,type);
                end
                Uz(Nx+1,:,:) = 0;  % this is reduntant, actually
%=========================================================================%
%                         boundary on the bottom
%=========================================================================%
%   interpolate
                yi = (yi(1:Ny)+yi(2:Ny+1))/2;
                u = [Uy(1,:,end); Uy(end,:,end)];
                y = [yi'; yi'];
                x = [ones(1,Ny)*xi(1); ones(1,Ny)*xi(end)];
                [yyi,xxi]=meshgrid(yi,xi);
                Uy(:,:,Nz+1) = interp2(y,x,u,yyi,xxi);
                Ux(:,:,Nz+1) = 0;
%=========================================================================%
%           now try to make a initial guess of the inner fields by
%                     interpolating the boundary
%=========================================================================%
% for now only do it in y               
                for iz=2:Nz
                    u = [Uy(1,:,iz); Uy(end,:,iz)];
                    Uy(:,:,iz) = interp2(y,x,u,yyi,xxi);
                end
        end
    case '2D'
        switch pol
            case 'XY' % N-S source 
%=========================================================================%
%               boundary on the X edges - should arange all BCs
%=========================================================================%
                for islice=1:Nx
                    res2d=reshape(res(islice,:,:),Ny,Nz,1);
                    Ux2d= mtfwd2d(Dy,Dz,res2d,per,type);
                    Ux(islice,:,:)=Ux2d;
                end
            case 'YX' % E-W source 
%=========================================================================%
%               boundary on the Y edges - should arange all BCs
%=========================================================================%
                for islice=1:Ny
                    res2d=reshape(res(:,islice,:),Nx,Nz,1);
                    Uy2d= mtfwd2d(Dx,Dz,res2d,per,type);
                    Uy(:,islice,:)=Uy2d;
                end           
        end
    case 'dirichlet'
        switch pol
            case 'XY' % N-S source 
%=========================================================================%
%                         source on the top
%=========================================================================%
                Ux(:,:,1) = 1;  % can be either E or B field
            case 'YX' % E-W source
%=========================================================================%
%                         source on the top
%=========================================================================%
                Uy(:,:,1) = 1;  % can be either E or B field
        end
end
return 