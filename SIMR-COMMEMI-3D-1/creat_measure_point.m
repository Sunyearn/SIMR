function [point,line_x,line_y] = creat_measure_point(block,ilevel,iblock,point,line_x,line_y)
% this code is for setting the coordinates of synthetic data points
%% Prepare
    Nza = block(ilevel).Nza(iblock);% the number of air layers in current block
	Cx = block(ilevel).Cx(iblock).vec;% node coordinates in x direction
	Cy = block(ilevel).Cy(iblock).vec;% node coordinates in y direction
	Cz = block(ilevel).Cz(iblock).vec((Nza+1):end);% eliminate the node coordinates of air layers in z direction
	Cxm = ( Cx(1:(end-1)) + Cx(2:end) )/2; 
	Cym = ( Cy(1:(end-1)) + Cy(2:end) )/2;
	Czm = ( Cz(1:(end-1)) + Cz(2:end) )/2;
    Nx = length(Cxm);
    Ny = length(Cym);
%% Struct "point"
    % the point sits at the approximate centre of whole domain
    pos_x = Cxm(ceil(Nx/2));
    pos_y = Cym(ceil(Ny/2));
    pos_z = Cz(1:10);
    [My,Mx,Mz] = meshgrid(pos_y,pos_x,pos_z);
    point(ilevel,iblock).x = Mx;
    point(ilevel,iblock).y = My;
    point(ilevel,iblock).z = Mz;
%% Struct "line_x"
    Ny = length(Cym);
    Cym_l = Cym(ceil(Ny/2))*ones(size(Cxm));
    Czm_l = Cz(1)*ones(size(Cxm));
    line_x(ilevel,iblock).x = Cxm;
    line_x(ilevel,iblock).y = Cym_l;
    line_x(ilevel,iblock).z = Czm_l;
%% Struct "line_y"
    Nx = length(Cxm);
    Cxm_l = Cxm(ceil(Nx/2)+2)*ones(size(Cym));
    Czm_l = Cz(1)*ones(size(Cym));
    line_y(ilevel,iblock).x = Cxm_l;
    line_y(ilevel,iblock).y = Cym;
    line_y(ilevel,iblock).z = Czm_l;
end