function [overlap] = find_div_overlap(Nx_big,Ny_big,Nz_big, Nx_small,Ny_small,Nz_small,ratio,mode)
    mode = upper(mode);
    switch mode
        case 'NODE'
            [Nedge,Nface,Nnode] = cal_fen(Nx_big,Ny_big,Nz_big);
            Nnode_big = Nnode;%
            [Nedge,Nface,Nnode] = cal_fen(Nx_small,Ny_small,Nz_small);
            Nnode_small = Nnode;%
            
            node_big = zeros(Nx_big+1,Ny_big+1,Nz_big+1);
            node_small = zeros(Nx_small+1,Ny_small+1,Nz_small+1);
            node_big(:,:,:) = 1;
            node_small(1:ratio(1):end,1:ratio(2):end,1:ratio(3):end) = 1;
            
            node_big = reshape(node_big,Nnode_big,1);
            node_small = reshape(node_small,Nnode_small,1);
            
            node_big_vec = find(node_big == 1);
            node_small_vec = find(node_small == 1);

            overlap = [node_big_vec,node_small_vec];
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
            % edge x
            edge_x_big(:,:,:) = 1;
            for ix = 1:ratio(1)
                edge_x_small(ix:ratio(1):(end-ratio(1)+ix),1:ratio(2):end,1:ratio(3):end) = ix;
            end
            edge_x_big = reshape(edge_x_big,Nedge_big(1),1);
            edge_x_small = reshape(edge_x_small,Nedge_small(1),1);
            edge_x_big_vec = find(edge_x_big == 1);
            for ix = 1:ratio(1)
                edge_x_small_temp = find(edge_x_small == ix);
                if ix == 1
                    edge_x_vec = [edge_x_big_vec,edge_x_small_temp];
                else
                    edge_x_vec = [edge_x_vec,edge_x_small_temp];
                end
            end
            % edge y
            edge_y_big(:,:,:) = 1;
            for iy = 1:ratio(2)
                edge_y_small(1:ratio(1):end,iy:ratio(2):(end-ratio(2)+iy),1:ratio(3):end) = iy;
            end
            edge_y_big = reshape(edge_y_big,Nedge_big(2),1);
            edge_y_small = reshape(edge_y_small,Nedge_small(2),1);
            edge_y_big_vec = find(edge_y_big == 1);
            edge_y_big_vec = edge_y_big_vec + Nedge_big(1);
            for iy = 1:ratio(2)
                edge_y_small_temp = find(edge_y_small == iy);
                edge_y_small_temp = edge_y_small_temp + Nedge_small(1);
                if iy == 1
                    edge_y_vec = [edge_y_big_vec,edge_y_small_temp];
                else
                    edge_y_vec = [edge_y_vec,edge_y_small_temp];
                end
            end
            % edge z
            edge_z_big(:,:,:) = 1;
            for iz = 1:ratio(3)
                edge_z_small(1:ratio(1):end,1:ratio(2):end,iz:ratio(3):(end-ratio(3)+iz)) = iz;
            end
            edge_z_big = reshape(edge_z_big,Nedge_big(3),1);
            edge_z_small = reshape(edge_z_small,Nedge_small(3),1);
            edge_z_big_vec = find(edge_z_big == 1);
            edge_z_big_vec = edge_z_big_vec + Nedge_big(1) + Nedge_big(2);
            for iz = 1:ratio(3)
                edge_z_small_temp = find(edge_z_small == iz);
                edge_z_small_temp = edge_z_small_temp + Nedge_small(1) + Nedge_small(2);
                if iz == 1
                    edge_z_vec = [edge_z_big_vec,edge_z_small_temp];
                else
                    edge_z_vec = [edge_z_vec,edge_z_small_temp];
                end
            end
            overlap.edge_x = edge_x_vec;
            overlap.edge_y = edge_y_vec;
            overlap.edge_z = edge_z_vec;
            
        otherwise
            warndlg('mode is wrong! (find_div_overlap.m)');
    end
return