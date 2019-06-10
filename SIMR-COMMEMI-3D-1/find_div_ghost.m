function [ghost2,ghost4,ghost8,ghost_control] = find_div_ghost(Nx_big,Ny_big,Nz_big, Nx_small,Ny_small,Nz_small,ratio,mode)
    mode = upper(mode);
    switch mode
        case 'NODE'
            [Nedge,Nface,Nnode] = cal_fen(Nx_big,Ny_big,Nz_big);
            Nnode_big = Nnode;%
            [Nedge,Nface,Nnode] = cal_fen(Nx_small,Ny_small,Nz_small);
            Nnode_small = Nnode;%
%% ghost 8
            % 准备
            node_big = zeros(Nx_big+1,Ny_big+1,Nz_big+1);
            node_small = zeros(Nx_small+1,Ny_small+1,Nz_small+1);
            node_big(1:(end-1),1:(end-1),1:(end-1)) = 1;
            for iz = 1 : (ratio(3)-1)
                for iy = 1 : (ratio(2)-1)
                    for ix = 1 : (ratio(1)-1)
                        node_small((1+ix):ratio(1):(end-ratio(1)+ix),  (1+iy):ratio(2):(end-ratio(2)+iy),...
                                                    (1+iz):ratio(3):(end-ratio(3)+iz) ) =...
                                                            ix + (iy-1)*(ratio(1)-1) + ...
                                                                 (iz-1)*(ratio(2)-1)*(ratio(3)-1);
                    end
                end
            end
            % 矩阵变形
            node_big = reshape(node_big,Nnode_big,1);
            node_small = reshape(node_small,Nnode_small,1);
            % 查询对应特征的位置
            node_big_temp = find(node_big == 1);
            add = [0,1,(Nx_big+1),(Nx_big+1+1), ...
                (Nx_big+1)*(Ny_big+1),((Nx_big+1)*(Ny_big+1) + 1), ...
                ((Nx_big+1)*(Ny_big+1)+(Nx_big+1)),((Nx_big+1)*(Ny_big+1)+(Nx_big+1) + 1)];
            for in = 1:8
                if in == 1
                    node_big_vec = node_big_temp;
                else
                    node_big_vec = [node_big_vec,(node_big_temp + add(in)) ];
                end
            end
            % 构建索引 ghost8
            for in = 1:( (ratio(1)-1)*(ratio(2)-1)*(ratio(3)-1) )
                node_small_vec = find(node_small == in);
                if in == 1
                    ghost8 = [node_big_vec,node_small_vec];
                else
                    ghost8 = [ghost8;
                              node_big_vec,node_small_vec];
                end
            end
 %% ghost 4
 % face x 面上四个点
            % 准备
            node_big = zeros(Nx_big+1,Ny_big+1,Nz_big+1);
            node_small = zeros(Nx_small+1,Ny_small+1,Nz_small+1);
            node_big(:,1:(end-1),1:(end-1)) = 1;
            for iz = 1:(ratio(3)-1)
                for iy = 1:(ratio(2)-1)
                    node_small(1:ratio(1):end, (1+iy):ratio(2):(end-ratio(2)+iy), ...
                                                        (1+iz):ratio(3):(end-ratio(3)+iz)) = ...
                                                            iy + (iz-1)*(ratio(2)-1);
                end
            end
            % 矩阵变形
            node_big = reshape(node_big,Nnode_big,1);
            node_small = reshape(node_small,Nnode_small,1);
            % 查询对应特征的位置
            node_big_temp = find(node_big == 1);
            add = [0,(Nx_big+1),(Nx_big+1)*(Ny_big+1),((Nx_big+1)*(Ny_big+1)+(Nx_big+1))];
            for in = 1:4
                if in == 1
                    node_big_vec = node_big_temp;
                else
                    node_big_vec = [node_big_vec,(node_big_temp + add(in)) ];
                end
            end
            % 构建索引 ghost4
            for in = 1:( (ratio(2)-1)*(ratio(3)-1) )
                node_small_vec = find(node_small == in);
                if in == 1
                    ghost4 = [node_big_vec,node_small_vec];
                else
                    ghost4 = [ghost4;
                              node_big_vec,node_small_vec];
                end
            end
 % face y 面上四个点
 % face z 面上四个点
 % edge x 棱上二个点
 % edge y 棱上二个点
 % edge z 棱上二个点
        case 'EDGE'
            [Nedge,Nface,Nnode] = cal_fen(Nx_big,Ny_big,Nz_big);
            Nedge_big = Nedge(1,:);%
            
            [Nedge,Nface,Nnode] = cal_fen(Nx_small,Ny_small,Nz_small);
            Nedge_small = Nedge(1,:);%
            
            edge_x_big = zeros(Nx_big,Ny_big+1,Nz_big+1);
            edge_y_big = zeros(Nx_big+1,Ny_big,Nz_big+1);
            edge_z_big = zeros(Nx_big+1,Ny_big+1,Nz_big);
            
            edge_x_small = zeros(Nx_small,Ny_small+1,Nz_small+1);
            edge_y_small = zeros(Nx_small+1,Ny_small,Nz_small+1);
            edge_z_small = zeros(Nx_small+1,Ny_small+1,Nz_small);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % edge x ghost 4 in y and z direction
            edge_x_big(:,1:(end-1),1:(end-1)) = edge_x_big(:,1:(end-1),1:(end-1)) + 1;
            edge_x_big(:,2:end,1:(end-1)) = edge_x_big(:,2:end,1:(end-1)) + 10;
            edge_x_big(:,1:(end-1),2:end) = edge_x_big(:,1:(end-1),2:end) + 1000i;
            edge_x_big(:,2:end,2:end) = edge_x_big(:,2:end,2:end) + 10000i;
            
            edge_x_big = reshape(edge_x_big,Nedge_big(1),1);
            index_edge_x_ghostz1y1 = find(real(edge_x_big)==1 | real(edge_x_big)==11);
            index_edge_x_ghostz1y2 = find(real(edge_x_big)>=10);
            index_edge_x_ghostz2y1 = find(imag(edge_x_big)== 1000 | imag(edge_x_big)== 11000);
            index_edge_x_ghostz2y2 = find(imag(edge_x_big)>= 10000);
            
            for iz = 1:(ratio(3)-1)
                for iy = 1:(ratio(2)-1)
                    for ix = 1:ratio(1)
                        edge_x_small(ix:ratio(1):(end-ratio(1)+ix), ...
                                   (1+iy):ratio(2):(end-ratio(2)+iy), ...
                                   (1+iz):ratio(3):(end-ratio(3)+iz)) = ix + (iy-1)*ratio(1) + (iz-1)*(ratio(2)-1)*ratio(1);
                    end
                end
            end
            edge_x_small = reshape(edge_x_small,Nedge_small(1),1);
            for ie = 1:(ratio(1)*(ratio(2)-1)*(ratio(3)-1))
                edge_x_small_temp = find(edge_x_small == ie);
                temp = [index_edge_x_ghostz1y1,index_edge_x_ghostz1y2,...
                                   index_edge_x_ghostz2y1,index_edge_x_ghostz2y2,edge_x_small_temp];
                if ie == 1
                    edge_x_vec4 = temp;
                else
                    edge_x_vec4 = [edge_x_vec4;temp];
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % edge x ghost 2 in y direction
            edge_x_big = zeros(Nx_big,Ny_big+1,Nz_big+1);
            edge_x_small = zeros(Nx_small,Ny_small+1,Nz_small+1);
            edge_x_big(:,1:(end-1),:) = edge_x_big(:,1:(end-1),:) + 1;
            edge_x_big(:,2:end,:) = edge_x_big(:,2:end,:) + 10;
            
            edge_x_big = reshape(edge_x_big,Nedge_big(1),1);
            index_edge_x_ghosty1 = find(real(edge_x_big)==1 | real(edge_x_big)==11);
            index_edge_x_ghosty2 = find(real(edge_x_big)>=10);
            
            for iy = 1:(ratio(2)-1)
                  for ix = 1:ratio(1)
                      edge_x_small(ix:ratio(1):(end-ratio(1)+ix), ...
                                 (1+iy):ratio(2):(end-ratio(2)+iy), 1:ratio(3):end) = ...
                                        ix + (iy-1)*ratio(1);
                  end
            end
            edge_x_small = reshape(edge_x_small,Nedge_small(1),1);
            for ie = 1:(ratio(1)*(ratio(2)-1))
                edge_x_small_temp = find(edge_x_small == ie);
                temp = [index_edge_x_ghosty1,index_edge_x_ghosty2,edge_x_small_temp];
                if ie == 1
                    edge_x_vec2 = temp;
                else
                    edge_x_vec2 = [edge_x_vec2;temp];
                end
            end
            xY = size(edge_x_vec2,1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % edge x ghost 2 in z direction
            edge_x_big = zeros(Nx_big,Ny_big+1,Nz_big+1);
            edge_x_small = zeros(Nx_small,Ny_small+1,Nz_small+1);
            edge_x_big(:,:,1:(end-1)) = edge_x_big(:,:,1:(end-1)) + 1;
            edge_x_big(:,:,2:end) = edge_x_big(:,:,2:end) + 10;
            
            edge_x_big = reshape(edge_x_big,Nedge_big(1),1);
            index_edge_x_ghostz1 = find(real(edge_x_big)==1 | real(edge_x_big)==11);
            index_edge_x_ghostz2 = find(real(edge_x_big)>=10);
            
            for iz = 1:(ratio(3)-1)
                  for ix = 1:ratio(1)
                      edge_x_small(ix:ratio(1):(end-ratio(1)+ix), ...
                                 1:ratio(2):end,(1+iz):ratio(3):(end-ratio(3)+iz)) = ...
                                        ix + (iz-1)*ratio(1);
                  end
            end
            edge_x_small = reshape(edge_x_small,Nedge_small(1),1);
            for ie = 1:(ratio(1)*(ratio(3)-1))
                edge_x_small_temp = find(edge_x_small == ie);
                temp = [index_edge_x_ghostz1,index_edge_x_ghostz2,edge_x_small_temp];
                
                edge_x_vec2 = [edge_x_vec2;temp];
            end
            xZ = size(edge_x_vec2,1) - xY;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% edge y ghost 4 in x and z direction
            edge_y_big(1:(end-1),:,1:(end-1)) = edge_y_big(1:(end-1),:,1:(end-1)) + 1;
            edge_y_big(2:end,:,1:(end-1)) = edge_y_big(2:end,:,1:(end-1)) + 10;
            edge_y_big(1:(end-1),:,2:end) = edge_y_big(1:(end-1),:,2:end) + 1000i;
            edge_y_big(2:end,:,2:end) = edge_y_big(2:end,:,2:end) + 10000i;
            
            edge_y_big = reshape(edge_y_big,Nedge_big(2),1);
            index_edge_y_ghostz1y1 = find(real(edge_y_big)==1 | real(edge_y_big)==11);
            index_edge_y_ghostz1y2 = find(real(edge_y_big)>=10);
            index_edge_y_ghostz2y1 = find(imag(edge_y_big)== 1000 | imag(edge_y_big)== 11000);
            index_edge_y_ghostz2y2 = find(imag(edge_y_big)>= 10000);
            
            for iz = 1:(ratio(3)-1)
                for iy = 1:ratio(2)
                    for ix = 1:(ratio(1)-1)
                        edge_y_small( (ix+1):ratio(1):(end-ratio(1)+ix), ...
                                        iy:ratio(2):(end-ratio(2)+iy), ...
                                        (1+iz):ratio(3):(end-ratio(3)+iz)) = ix + (iy-1)*(ratio(1)-1) + (iz-1)*(ratio(1)-1)*ratio(2);
                    end
                end
            end
            edge_y_small = reshape(edge_y_small,Nedge_small(2),1);
            for ie = 1:(ratio(2)*(ratio(1)-1)*(ratio(3)-1))
                edge_y_small_temp = find(edge_y_small == ie);
                temp = [index_edge_y_ghostz1y1,index_edge_y_ghostz1y2,...
                                   index_edge_y_ghostz2y1,index_edge_y_ghostz2y2,edge_y_small_temp];
                if ie == 1
                    edge_y_vec4 = temp;
                else
                    edge_y_vec4 = [edge_y_vec4;temp];
                end
            end
            edge_y_vec4(:,1:4) = edge_y_vec4(:,1:4) + Nedge_big(1);
            edge_y_vec4(:,5) = edge_y_vec4(:,5) + Nedge_small(1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % edge y ghost 2 in x direction
            edge_y_big = zeros(Nx_big+1,Ny_big,Nz_big+1);
            edge_y_small = zeros(Nx_small+1,Ny_small,Nz_small+1);
            edge_y_big(1:(end-1),:,:) = edge_y_big(1:(end-1),:,:) + 1;
            edge_y_big(2:end,:,:) = edge_y_big(2:end,:,:) + 10;
            
            edge_y_big = reshape(edge_y_big,Nedge_big(2),1);
            index_edge_y_ghostx1 = find(real(edge_y_big)==1 | real(edge_y_big)==11);
            index_edge_y_ghostx2 = find(real(edge_y_big)>=10);
            
            for iy = 1:ratio(2)
                  for ix = 1:(ratio(1)-1)
                      edge_y_small((1+ix):ratio(1):(end-ratio(1)+ix), ...
                                 iy:ratio(2):(end-ratio(2)+iy), 1:ratio(3):end) = ...
                                        ix + (iy-1)*(ratio(1)-1);
                  end
            end
            edge_y_small = reshape(edge_y_small,Nedge_small(2),1);
            for ie = 1:(ratio(2)*(ratio(1)-1))
                edge_y_small_temp = find(edge_y_small == ie);
                temp = [index_edge_y_ghostx1,index_edge_y_ghostx2,edge_y_small_temp];
                if ie == 1
                    edge_y_vec2 = temp;
                else
                    edge_y_vec2 = [edge_y_vec2;temp];
                end
            end
            yX = size(edge_y_vec2,1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % edge y ghost 2 in z direction
            edge_y_big = zeros(Nx_big+1,Ny_big,Nz_big+1);
            edge_y_small = zeros(Nx_small+1,Ny_small,Nz_small+1);
            edge_y_big(:,:,1:(end-1)) = edge_y_big(:,:,1:(end-1)) + 1;
            edge_y_big(:,:,2:end) = edge_y_big(:,:,2:end) + 10;
            
            edge_y_big = reshape(edge_y_big,Nedge_big(2),1);
            index_edge_y_ghostz1 = find(real(edge_y_big)==1 | real(edge_y_big)==11);
            index_edge_y_ghostz2 = find(real(edge_y_big)>=10);
            
            for iz = 1:(ratio(3)-1)
                  for iy = 1:ratio(2)
                      edge_y_small(1:ratio(1):end, ...
                                 iy:ratio(2):(end-ratio(2)+iy),(1+iz):ratio(3):(end-ratio(3)+iz)) = ...
                                        iy + (iz-1)*ratio(2);
                  end
            end
            edge_y_small = reshape(edge_y_small,Nedge_small(2),1);
            for ie = 1:(ratio(2)*(ratio(3)-1))
                edge_y_small_temp = find(edge_y_small == ie);
                temp = [index_edge_y_ghostz1,index_edge_y_ghostz2,edge_y_small_temp];
                
                edge_y_vec2 = [edge_y_vec2;temp];
            end
            yZ = size(edge_y_vec2,1) - yX;
            edge_y_vec2(:,1:2) = edge_y_vec2(:,1:2) + Nedge_big(1);
            edge_y_vec2(:,3) = edge_y_vec2(:,3) + Nedge_small(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% edge z ghost 4 in x and y direction
            edge_z_big(1:(end-1),1:(end-1),:) = edge_z_big(1:(end-1),1:(end-1),:) + 1;
            edge_z_big(2:end,1:(end-1),:) = edge_z_big(2:end,1:(end-1),:) + 10;
            edge_z_big(1:(end-1),2:end,:) = edge_z_big(1:(end-1),2:end,:) + 1000i;
            edge_z_big(2:end,2:end,:) = edge_z_big(2:end,2:end,:) + 10000i;
            
            edge_z_big = reshape(edge_z_big,Nedge_big(3),1);
            index_edge_z_ghostx1y1 = find(real(edge_z_big)==1 | real(edge_z_big)==11);
            index_edge_z_ghostx2y1 = find(real(edge_z_big)>=10);
            index_edge_z_ghostx1y2 = find(imag(edge_z_big)== 1000 | imag(edge_z_big)== 11000);
            index_edge_z_ghostx2y2 = find(imag(edge_z_big)>= 10000);
            
            for iz = 1:ratio(3)
                for iy = 1:(ratio(2)-1)
                    for ix = 1:(ratio(1)-1)
                        edge_z_small( (ix+1):ratio(1):(end-ratio(1)+ix), ...
                            (iy+1):ratio(2):(end-ratio(2)+iy), ...
                            iz:ratio(3):(end-ratio(3)+iz)) = ix + (iy-1)*(ratio(1)-1) + (iz-1)*(ratio(1)-1)*(ratio(2)-1);
                    end
                end
            end
            edge_z_small = reshape(edge_z_small,Nedge_small(3),1);
            for ie = 1:(ratio(3)*(ratio(1)-1)*(ratio(2)-1))
                edge_z_small_temp = find(edge_z_small == ie);
                temp = [index_edge_z_ghostx1y1,index_edge_z_ghostx2y1,...
                                   index_edge_z_ghostx1y2,index_edge_z_ghostx2y2,edge_z_small_temp];
                if ie == 1
                    edge_z_vec4 = temp;
                else
                    edge_z_vec4 = [edge_z_vec4;temp];
                end
            end
            edge_z_vec4(:,1:4) = edge_z_vec4(:,1:4) + Nedge_big(1) + Nedge_big(2);
            edge_z_vec4(:,5) = edge_z_vec4(:,5) + Nedge_small(1) + Nedge_small(2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % edge z ghost 2 in x direction
            edge_z_big = zeros(Nx_big+1,Ny_big+1,Nz_big);
            edge_z_small = zeros(Nx_small+1,Ny_small+1,Nz_small);
            edge_z_big(1:(end-1),:,:) = edge_z_big(1:(end-1),:,:) + 1;
            edge_z_big(2:end,:,:) = edge_z_big(2:end,:,:) + 10;
            
            edge_z_big = reshape(edge_z_big,Nedge_big(3),1);
            index_edge_z_ghostx1 = find(real(edge_z_big)==1 | real(edge_z_big)==11);
            index_edge_z_ghostx2 = find(real(edge_z_big)>=10);
            
            for iz = 1:ratio(3)
                  for ix = 1:(ratio(1)-1)
                      edge_z_small((1+ix):ratio(1):(end-ratio(1)+ix), ...
                                 1:ratio(2):end, iz:ratio(3):(end-ratio(3)+iz)) = ...
                                        ix + (iz-1)*(ratio(1)-1);
                  end
            end
            edge_z_small = reshape(edge_z_small,Nedge_small(3),1);
            for ie = 1:(ratio(3)*(ratio(1)-1))
                edge_z_small_temp = find(edge_z_small == ie);
                temp = [index_edge_z_ghostx1,index_edge_z_ghostx2,edge_z_small_temp];
                if ie == 1
                    edge_z_vec2 = temp;
                else
                    edge_z_vec2 = [edge_z_vec2;temp];
                end
            end
            zX = size(edge_z_vec2,1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % edge z ghost 2 in y direction
            edge_z_big = zeros(Nx_big+1,Ny_big+1,Nz_big);
            edge_z_small = zeros(Nx_small+1,Ny_small+1,Nz_small);
            edge_z_big(:,1:(end-1),:) = edge_z_big(:,1:(end-1),:) + 1;
            edge_z_big(:,2:end,:) = edge_z_big(:,2:end,:) + 10;
            
            edge_z_big = reshape(edge_z_big,Nedge_big(3),1);
            index_edge_z_ghosty1 = find(real(edge_z_big)==1 | real(edge_z_big)==11);
            index_edge_z_ghosty2 = find(real(edge_z_big)>=10);
            
            for iz = 1:ratio(3)
                  for iy = 1:(ratio(2)-1)
                      edge_z_small(1:ratio(1):end, ...
                                 (1+iy):ratio(2):(end-ratio(2)+iy), iz:ratio(3):(end-ratio(3)+iz)) = ...
                                        iy + (iz-1)*(ratio(2)-1);
                  end
            end
            edge_z_small = reshape(edge_z_small,Nedge_small(3),1);
            for ie = 1:(ratio(3)*(ratio(2)-1))
                edge_z_small_temp = find(edge_z_small == ie);
                temp = [index_edge_z_ghosty1,index_edge_z_ghosty2,edge_z_small_temp];

                edge_z_vec2 = [edge_z_vec2;temp];
            end
            zY = size(edge_z_vec2,1) - zX;
            edge_z_vec2(:,1:2) = edge_z_vec2(:,1:2) + Nedge_big(1) + Nedge_big(2);
            edge_z_vec2(:,3) = edge_z_vec2(:,3) + Nedge_small(1) + Nedge_small(2);
            
            ghost2 = [edge_x_vec2;edge_y_vec2;edge_z_vec2];
            ghost4 = [edge_x_vec4;edge_y_vec4;edge_z_vec4];
            ghost8 = [];
            
            ghost_control.e2 = [xY,xZ,yX,yZ,zX,zY];
            ghost_control.e4 = [size(edge_x_vec4,1),size(edge_y_vec4,1),size(edge_z_vec4,1)];
            ghost_control.e8 = [0,0,0];
        otherwise
            warndlg('mode is wrong! (find_div_overlap.m)');
    end
return