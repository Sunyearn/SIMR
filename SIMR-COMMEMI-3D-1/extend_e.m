function [model,block] = extend_e(model,block,ilevel,iblo,iper)
% The "_now" is about finer-grid, "_pre" is about coarser-grid. The finer-grid's fields are
% interpolated of replaced by the coarser-grid's fields.
block_source = block(ilevel).source(iblo,:);
if block_source(end) ~= 0
%% prepare
    Ex_pre = model(ilevel-1).Ex(:,:,iper);
    Ey_pre = model(ilevel-1).Ey(:,:,iper);
    Ez_pre = model(ilevel-1).Ez(:,:,iper);
    
    Nx_pre = length(model(ilevel-1).Dx);
    Ny_pre = length(model(ilevel-1).Dy);
    Nz_pre = length(model(ilevel-1).Dz);
    Cx_pre = model(ilevel-1).Cx;
    Cy_pre = model(ilevel-1).Cy;
    Cz_pre = model(ilevel-1).Cz;
    Dx_pre = model(ilevel-1).Dx;
    Dy_pre = model(ilevel-1).Dy;
    Dz_pre = model(ilevel-1).Dz;
    Res_pre = model(ilevel-1).res;
    Nedge_pre = model(ilevel-1).fen.edge(1,:);
    Sum_edge_pre = sum(Nedge_pre);
    
    Nx_now = length(model(ilevel).Dx);
    Ny_now = length(model(ilevel).Dy);
    Nz_now = length(model(ilevel).Dz);
    Cx_now = model(ilevel).Cx;
    Cy_now = model(ilevel).Cy;
    Cz_now = model(ilevel).Cz;
    Dx_now = model(ilevel).Dx;
    Dy_now = model(ilevel).Dy;
    Dz_now = model(ilevel).Dz;
    Res_now = model(ilevel).res;
    Nedge_now = model(ilevel).fen.edge(1,:);
    Sum_edge_now = sum(Nedge_now); 
    
    ratio = (model(ilevel).ratio) ./ (model(ilevel-1).ratio);
    overlap_edge = find_div_overlap(Nx_pre,Ny_pre,Nz_pre, Nx_now,Ny_now,Nz_now,ratio,'EDGE');
    [ghost.edge2,ghost.edge4,ghost.edge8,ghost_control] = find_div_ghost(Nx_pre,Ny_pre,Nz_pre, ...
                                                                         Nx_now,Ny_now,Nz_now,ratio,'edge');
    [coor_edge_X_pre,coor_edge_Y_pre,coor_edge_Z_pre] = coor_edges(Cx_pre,Cy_pre,Cz_pre);
    coor_edge_pre = [coor_edge_X_pre;coor_edge_Y_pre;coor_edge_Z_pre];
    [coor_edge_X_now,coor_edge_Y_now,coor_edge_Z_now] = coor_edges(Cx_now,Cy_now,Cz_now);
    coor_edge_now = [coor_edge_X_now;coor_edge_Y_now;coor_edge_Z_now];
    % the average sigma
    Sigma_pre = 1./Res_pre;
    Vcell_pre = cal_cellvol(Dx_pre,Dy_pre,Dz_pre,'cell');
    sigma_avg_pre = ave_sigma(Sigma_pre,Vcell_pre,Nedge_pre,'edge'); %构建棱边元上的模型参数矩阵
    Sigma_now = 1./Res_now;
    Vcell_now = cal_cellvol(Dx_now,Dy_now,Dz_now,'cell');
    sigma_avg_now = ave_sigma(Sigma_now,Vcell_now,Nedge_now,'edge'); %构建棱边元上的模型参数矩阵
    
for ib = 1:(block_source(end)) % number of blocks in upper level grid
%% replacement and interpretation
    for ipol = 1:2
        Ex = Ex_pre(ipol,block_source(ib)).matrix;
        Ey = Ey_pre(ipol,block_source(ib)).matrix;
        Ez = Ez_pre(ipol,block_source(ib)).matrix;
        ex = zeros(Nx_now,Ny_now+1,Nz_now+1);
        ey = zeros(Nx_now+1,Ny_now,Nz_now+1);
        ez = zeros(Nx_now+1,Ny_now+1,Nz_now);
        
        Ex = reshape(Ex,Nedge_pre(1),1);
        Ey = reshape(Ey,Nedge_pre(2),1);
        Ez = reshape(Ez,Nedge_pre(3),1);
        ex = reshape(ex,Nedge_now(1),1);
        ey = reshape(ey,Nedge_now(2),1);
        ez = reshape(ez,Nedge_now(3),1);
        
        % replace
        [ex,ey,ez] = replace_overlap_edge_value(Ex,Ey,Ez,ex,ey,ez,sigma_avg_pre,sigma_avg_now,...
                                                                    Nedge_pre,Nedge_now,overlap_edge,ratio);
        % interpolate
        [ex,ey,ez] = replace_ghost_edge_value(Ex,Ey,Ez,ex,ey,ez,sigma_avg_pre,sigma_avg_now, ...
                                            Nedge_pre,Nedge_now,coor_edge_pre,coor_edge_now,ghost,ghost_control,ratio);
        ex = reshape(ex,Nx_now,Ny_now+1,Nz_now+1);
        ey = reshape(ey,Nx_now+1,Ny_now,Nz_now+1);
        ez = reshape(ez,Nx_now+1,Ny_now+1,Nz_now);
        if ib == 1
            Ex_now(ipol,iblo).matrix = ex;
            Ey_now(ipol,iblo).matrix = ey;
            Ez_now(ipol,iblo).matrix = ez;
        else
            ratio_now = (model(end).ratio) ./ (model(ilevel).ratio);
            block_area_pre = block(ilevel-1).area(block_source(ib),:);
            block_area_pre(1:2) = block_area_pre(1:2)/ratio_now(1);
            block_area_pre(3:4) = block_area_pre(3:4)/ratio_now(2);
            block_area_pre(5:6) = block_area_pre(5:6)/ratio_now(3);
            %
            x_range = (block_area_pre(1)+1) : (Nx_now - block_area_pre(2));
            y_range = (block_area_pre(3)+1) : (Ny_now+1 - block_area_pre(4));
            z_range = (block_area_pre(5)+1) : (Nz_now+1 - block_area_pre(6));
            Ex_now(ipol,iblo).matrix = ex(x_range,y_range,z_range);

            x_range = (block_area_pre(1)+1) : (Nx_now+1 - block_area_pre(2));
            y_range = (block_area_pre(3)+1) : (Ny_now - block_area_pre(4));
            z_range = (block_area_pre(5)+1) : (Nz_now+1 - block_area_pre(6));
            Ey_now(ipol,iblo).matrix = ey(x_range,y_range,z_range);
            
            x_range = (block_area_pre(1)+1) : (Nx_now+1 - block_area_pre(2));
            y_range = (block_area_pre(3)+1) : (Ny_now+1 - block_area_pre(4));
            z_range = (block_area_pre(5)+1) : (Nz_now - block_area_pre(6));
            Ez_now(ipol,iblo).matrix = ez(x_range,y_range,z_range);
        end
    end
end
    % fill in
    model(ilevel).Ex(1,iblo,iper) = Ex_now(1,iblo);
    model(ilevel).Ex(2,iblo,iper) = Ex_now(2,iblo);
    model(ilevel).Ey(1,iblo,iper) = Ey_now(1,iblo);
    model(ilevel).Ey(2,iblo,iper) = Ey_now(2,iblo);
    model(ilevel).Ez(1,iblo,iper) = Ez_now(1,iblo);
    model(ilevel).Ez(2,iblo,iper) = Ez_now(2,iblo);
    ratio_now = (model(end).ratio) ./ (model(ilevel).ratio);
    block_area_now = block(ilevel).area(iblo,:);
    block_area_now(1:2) = block_area_now(1:2)/ratio_now(1);
    block_area_now(3:4) = block_area_now(3:4)/ratio_now(2);
    block_area_now(5:6) = block_area_now(5:6)/ratio_now(3);
    x_range = (block_area_now(1)+1) : (Nx_now - block_area_now(2));
    y_range = (block_area_now(3)+1) : (Ny_now+1 - block_area_now(4));
    z_range = (block_area_now(5)+1) : (Nz_now+1 - block_area_now(6));
    block(ilevel).Ex(1,iblo,iper).matrix = Ex_now(1,iblo).matrix(x_range,y_range,z_range);
    block(ilevel).Ex(2,iblo,iper).matrix = Ex_now(2,iblo).matrix(x_range,y_range,z_range);
    x_range = (block_area_now(1)+1) : (Nx_now+1 - block_area_now(2));
    y_range = (block_area_now(3)+1) : (Ny_now - block_area_now(4));
    z_range = (block_area_now(5)+1) : (Nz_now+1 - block_area_now(6));
    block(ilevel).Ey(1,iblo,iper).matrix = Ey_now(1,iblo).matrix(x_range,y_range,z_range);
    block(ilevel).Ey(2,iblo,iper).matrix = Ey_now(2,iblo).matrix(x_range,y_range,z_range);
    x_range = (block_area_now(1)+1) : (Nx_now+1 - block_area_now(2));
    y_range = (block_area_now(3)+1) : (Ny_now+1 - block_area_now(4));
    z_range = (block_area_now(5)+1) : (Nz_now - block_area_now(6));
    block(ilevel).Ez(1,iblo,iper).matrix = Ez_now(1,iblo).matrix(x_range,y_range,z_range);
    block(ilevel).Ez(2,iblo,iper).matrix = Ez_now(2,iblo).matrix(x_range,y_range,z_range);
else
    return;
end
end