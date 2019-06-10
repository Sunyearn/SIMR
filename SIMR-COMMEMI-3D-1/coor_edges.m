function [coor_edge_X,coor_edge_Y,coor_edge_Z] = coor_edges(Cx,Cy,Cz)

    coor_x = Cx;
    coor_y = Cy;
    coor_z = Cz;
    
    Nx = length(Cx)-1;
    Ny = length(Cy)-1;
    Nz = length(Cz)-1;
    
    middle_x = (coor_x(1:(end-1)) + coor_x(2:end) )/2;
    middle_y = (coor_y(1:(end-1)) + coor_y(2:end) )/2;
    middle_z = (coor_z(1:(end-1)) + coor_z(2:end) )/2;
    
    [coor_edge_X_y,coor_edge_X_x,coor_edge_X_z] = meshgrid(coor_y,middle_x,coor_z);
    [coor_edge_Y_y,coor_edge_Y_x,coor_edge_Y_z] = meshgrid(middle_y,coor_x,coor_z);
    [coor_edge_Z_y,coor_edge_Z_x,coor_edge_Z_z] = meshgrid(coor_y,coor_x,middle_z);
    
    coor_edge_X_x = reshape(coor_edge_X_x,Nx*(Ny+1)*(Nz+1),1);
    coor_edge_X_y = reshape(coor_edge_X_y,Nx*(Ny+1)*(Nz+1),1);
    coor_edge_X_z = reshape(coor_edge_X_z,Nx*(Ny+1)*(Nz+1),1);
    
    coor_edge_Y_x = reshape(coor_edge_Y_x,(Nx+1)*Ny*(Nz+1),1);
    coor_edge_Y_y = reshape(coor_edge_Y_y,(Nx+1)*Ny*(Nz+1),1);
    coor_edge_Y_z = reshape(coor_edge_Y_z,(Nx+1)*Ny*(Nz+1),1);
    
    coor_edge_Z_x = reshape(coor_edge_Z_x,(Nx+1)*(Ny+1)*Nz,1);
    coor_edge_Z_y = reshape(coor_edge_Z_y,(Nx+1)*(Ny+1)*Nz,1);
    coor_edge_Z_z = reshape(coor_edge_Z_z,(Nx+1)*(Ny+1)*Nz,1);
    
    coor_edge_X = [coor_edge_X_x,coor_edge_X_y,coor_edge_X_z];
    coor_edge_Y = [coor_edge_Y_x,coor_edge_Y_y,coor_edge_Y_z];
    coor_edge_Z = [coor_edge_Z_x,coor_edge_Z_y,coor_edge_Z_z];
end