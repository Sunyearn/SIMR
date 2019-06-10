function [Nedge,Nface,Nnode] = cal_fen(Nx,Ny,Nz)
    Nedge_x = Nx*(Ny+1)*(Nz+1);
    Nedge_y = (Nx+1)*Ny*(Nz+1);
    Nedge_z = (Nx+1)*(Ny+1)*Nz;

    Nface_x = (Nx+1)*Ny*Nz;
    Nface_y = Nx*(Ny+1)*Nz;
    Nface_z = Nx*Ny*(Nz+1);

    Nedge = [Nedge_x,Nedge_y,Nedge_z;...
             Nedge_x,(Nedge_x+Nedge_y),(Nedge_x+Nedge_y+Nedge_z)];
    Nface = [Nface_x,Nface_y,Nface_z;...
             Nface_x,(Nface_x+Nface_y),(Nface_x+Nface_y+Nface_z)];
    Nnode = (Nx+1)*(Ny+1)*(Nz+1);
return