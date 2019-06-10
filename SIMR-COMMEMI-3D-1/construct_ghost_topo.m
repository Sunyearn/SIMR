function [ghost_topo] = construct_ghost_topo(sigma_avg_pre,sigma_avg_now,coor_edge_pre,coor_edge_now,...
                                                                                                        ghost,ghost_control)
%% edge2 GHOST
    if ~isempty(ghost.edge2)
        ghost_temp = ghost.edge2;
        topo_temp2 = zeros(size(ghost_temp));
        topo_temp2(:,3) = ghost_temp(:,3);
        
        control_temp = ghost_control.e2;
        control_sum = cumsum(control_temp);
        
        if control_temp(1)~=0
            % topo about length 根据棱边的 坐标值 调整散度算子的拓扑结构
            x_index_Y = ghost_temp(1:control_sum(1),:);
            Lall = abs( coor_edge_pre(x_index_Y(:,1),2) - coor_edge_pre(x_index_Y(:,2),2) );
            L1 = abs( coor_edge_pre(x_index_Y(:,1),2) - coor_edge_now(x_index_Y(:,3),2) );
            L2 = abs( coor_edge_pre(x_index_Y(:,2),2) - coor_edge_now(x_index_Y(:,3),2) );
            
            temp_1 = L2./Lall;
            temp_2 = L1./Lall;
            topo_temp2(1:control_sum(1),1) = temp_1;
            topo_temp2(1:control_sum(1),2) = temp_2;
            % topo about sigma 根据棱边的 sigma值 调整散度算子的拓扑结构
%             temp1 = sigma_avg_pre( x_index_Y(:,1) )./sigma_avg_now( x_index_Y(:,3) );
%             temp2 = sigma_avg_pre( x_index_Y(:,2) )./sigma_avg_now( x_index_Y(:,3) );
            temp1 = 1;
            temp2 = 1;
            topo_temp2(1:control_sum(1),1) = topo_temp2(1:control_sum(1),1).*temp1;
            topo_temp2(1:control_sum(1),2) = topo_temp2(1:control_sum(1),2).*temp2;
        else
            disp('There is no ghost edge x2 in Y direction');
        end
        
        if control_temp(2)~=0
            % topo about length 根据棱边的 坐标值 调整散度算子的拓扑结构
            x_index_Z = ghost_temp((control_sum(1)+1):control_sum(2),:);
            Lall = abs( coor_edge_pre(x_index_Z(:,1),3) - coor_edge_pre(x_index_Z(:,2),3) );
            L1 = abs( coor_edge_pre(x_index_Z(:,1),3) - coor_edge_now(x_index_Z(:,3),3) );
            L2 = abs( coor_edge_pre(x_index_Z(:,2),3) - coor_edge_now(x_index_Z(:,3),3) );
            
            temp_1 = L2./Lall;
            temp_2 = L1./Lall;
            topo_temp2((control_sum(1)+1):control_sum(2),1) = temp_1;
            topo_temp2((control_sum(1)+1):control_sum(2),2) = temp_2;
            % topo about sigma 根据棱边的 sigma值 调整散度算子的拓扑结构
%             temp1 = sigma_avg_pre( x_index_Z(:,1) )./sigma_avg_now( x_index_Z(:,3) );
%             temp2 = sigma_avg_pre( x_index_Z(:,2) )./sigma_avg_now( x_index_Z(:,3) );
            temp1 = 1;
            temp2 = 1;
            topo_temp2( (control_sum(1)+1):control_sum(2),1 ) = topo_temp2( (control_sum(1)+1):control_sum(2),1 ).*temp1;
            topo_temp2( (control_sum(1)+1):control_sum(2),2 ) = topo_temp2( (control_sum(1)+1):control_sum(2),2 ).*temp2;
        else
            disp('There is no ghost edge x2 in Z direction');
        end
        
        if control_temp(3)~=0
            % topo about length 根据棱边的 坐标值 调整散度算子的拓扑结构
            y_index_X = ghost_temp((control_sum(2)+1):control_sum(3),:);
            Lall = abs( coor_edge_pre(y_index_X(:,1),1) - coor_edge_pre(y_index_X(:,2),1) );
            L1 = abs( coor_edge_pre(y_index_X(:,1),1) - coor_edge_now(y_index_X(:,3),1) );
            L2 = abs( coor_edge_pre(y_index_X(:,2),1) - coor_edge_now(y_index_X(:,3),1) );
            
            temp_1 = L2./Lall;
            temp_2 = L1./Lall;
            topo_temp2((control_sum(2)+1):control_sum(3),1) = temp_1;
            topo_temp2((control_sum(2)+1):control_sum(3),2) = temp_2;
            % topo about sigma 根据棱边的 sigma值 调整散度算子的拓扑结构
%             temp1 = sigma_avg_pre( y_index_X(:,1) )./sigma_avg_now( y_index_X(:,3) );
%             temp2 = sigma_avg_pre( y_index_X(:,2) )./sigma_avg_now( y_index_X(:,3) );
            temp1 = 1;
            temp2 = 1;
            topo_temp2( (control_sum(2)+1):control_sum(3),1 ) = topo_temp2( (control_sum(2)+1):control_sum(3),1 ).*temp1;
            topo_temp2( (control_sum(2)+1):control_sum(3),2 ) = topo_temp2( (control_sum(2)+1):control_sum(3),2 ).*temp2;
        else
            disp('There is no ghost edge y2 in X direction');
        end
        
        if control_temp(4)~=0
            % topo about length 根据棱边的 坐标值 调整散度算子的拓扑结构
            y_index_Z = ghost_temp((control_sum(3)+1):control_sum(4),:);
            Lall = abs( coor_edge_pre(y_index_Z(:,1),3) - coor_edge_pre(y_index_Z(:,2),3) );
            L1 = abs( coor_edge_pre(y_index_Z(:,1),3) - coor_edge_now(y_index_Z(:,3),3) );
            L2 = abs( coor_edge_pre(y_index_Z(:,2),3) - coor_edge_now(y_index_Z(:,3),3) );
            
            temp_1 = L2./Lall;
            temp_2 = L1./Lall;
            topo_temp2((control_sum(3)+1):control_sum(4),1) = temp_1;
            topo_temp2((control_sum(3)+1):control_sum(4),2) = temp_2;
            % topo about sigma 根据棱边的 sigma值 调整散度算子的拓扑结构
%             temp1 = sigma_avg_pre( y_index_Z(:,1) )./sigma_avg_now( y_index_Z(:,3) );
%             temp2 = sigma_avg_pre( y_index_Z(:,2) )./sigma_avg_now( y_index_Z(:,3) );
            temp1 = 1;
            temp2 = 1;
            topo_temp2( (control_sum(3)+1):control_sum(4),1 ) = topo_temp2( (control_sum(3)+1):control_sum(4),1 ).*temp1;
            topo_temp2( (control_sum(3)+1):control_sum(4),2 ) = topo_temp2( (control_sum(3)+1):control_sum(4),2 ).*temp2;
        else
            disp('There is no ghost edge y2 in Z direction');
        end
        
        if control_temp(5)~=0
            % topo about length 根据棱边的 坐标值 调整散度算子的拓扑结构
            z_index_X = ghost_temp( (control_sum(4)+1):control_sum(5),:);
            Lall = abs( coor_edge_pre(z_index_X(:,1),1) - coor_edge_pre(z_index_X(:,2),1) );
            L1 = abs( coor_edge_pre(z_index_X(:,1),1) - coor_edge_now(z_index_X(:,3),1) );
            L2 = abs( coor_edge_pre(z_index_X(:,2),1) - coor_edge_now(z_index_X(:,3),1) );
            
            temp_1 = L2./Lall;
            temp_2 = L1./Lall;
            topo_temp2((control_sum(4)+1):control_sum(5),1) = temp_1;
            topo_temp2((control_sum(4)+1):control_sum(5),2) = temp_2;
            % topo about sigma 根据棱边的 sigma值 调整散度算子的拓扑结构
%             temp1 = sigma_avg_pre( z_index_X(:,1) )./sigma_avg_now( z_index_X(:,3) );
%             temp2 = sigma_avg_pre( z_index_X(:,2) )./sigma_avg_now( z_index_X(:,3) );
            temp1 = 1;
            temp2 = 1;
            topo_temp2( (control_sum(4)+1):control_sum(5),1 ) = topo_temp2( (control_sum(4)+1):control_sum(5),1 ).*temp1;
            topo_temp2( (control_sum(4)+1):control_sum(5),2 ) = topo_temp2( (control_sum(4)+1):control_sum(5),2 ).*temp2;
        else
            disp('There is no ghost edge z2 in X direction');
        end
        
        if control_temp(6)~=0
            % topo about length 根据棱边的 坐标值 调整散度算子的拓扑结构
            z_index_Y = ghost_temp((control_sum(5)+1):control_sum(6),:);
            Lall = abs( coor_edge_pre(z_index_Y(:,1),2) - coor_edge_pre(z_index_Y(:,2),2) );
            L1 = abs( coor_edge_pre(z_index_Y(:,1),2) - coor_edge_now(z_index_Y(:,3),2) );
            L2 = abs( coor_edge_pre(z_index_Y(:,2),2) - coor_edge_now(z_index_Y(:,3),2) );
            
            temp_1 = L2./Lall;
            temp_2 = L1./Lall;
            topo_temp2((control_sum(5)+1):control_sum(6),1) = temp_1;
            topo_temp2((control_sum(5)+1):control_sum(6),2) = temp_2;
            % topo about sigma 根据棱边的 sigma值 调整散度算子的拓扑结构
%             temp1 = sigma_avg_pre( z_index_Y(:,1) )./sigma_avg_now( z_index_Y(:,3) );
%             temp2 = sigma_avg_pre( z_index_Y(:,2) )./sigma_avg_now( z_index_Y(:,3) );
            temp1 = 1;
            temp2 = 1;
            topo_temp2( (control_sum(5)+1):control_sum(6),1 ) = topo_temp2( (control_sum(5)+1):control_sum(6),1 ).*temp1;
            topo_temp2( (control_sum(5)+1):control_sum(6),2 ) = topo_temp2( (control_sum(5)+1):control_sum(6),2 ).*temp2;
        else
            disp('There is no ghost edge z2 in Y direction');
        end
        ghost_topo.edge2 = topo_temp2;
    else
        ghost_topo.edge2 = [];
        warndlg('There is no ghost edge 2 !(construct_ghost_topo.m)');
    end
%% edge4 GHOST
    if ~isempty(ghost.edge4)
        ghost_temp = ghost.edge4;
        topo_temp4 = zeros(size(ghost_temp));
        topo_temp4(:,5) = ghost_temp(:,5);
        
        control_temp = ghost_control.e4;
        control_sum = cumsum(control_temp);
        
        if control_temp(1)~=0
            % topo about length 根据棱边的 坐标值 调整散度算子的拓扑结构
            x_index_YZ = ghost_temp(1:control_sum(1),:);
            Sall_y = abs( coor_edge_pre(x_index_YZ(:,1),2) - coor_edge_pre(x_index_YZ(:,2),2) );
            Sall_z = abs( coor_edge_pre(x_index_YZ(:,1),3) - coor_edge_pre(x_index_YZ(:,3),3) );
            Sall   = Sall_y.*Sall_z;
            L_y1 = abs( coor_edge_pre(x_index_YZ(:,1),2) - coor_edge_now(x_index_YZ(:,5),2) );
            L_z1 = abs( coor_edge_pre(x_index_YZ(:,1),3) - coor_edge_now(x_index_YZ(:,5),3) );
            S1   = L_y1.* L_z1;
            L_y2 = abs( coor_edge_pre(x_index_YZ(:,2),2) - coor_edge_now(x_index_YZ(:,5),2) );
            L_z2 = abs( coor_edge_pre(x_index_YZ(:,2),3) - coor_edge_now(x_index_YZ(:,5),3) );
            S2   = L_y2.* L_z2;
            L_y3 = abs( coor_edge_pre(x_index_YZ(:,3),2) - coor_edge_now(x_index_YZ(:,5),2) );
            L_z3 = abs( coor_edge_pre(x_index_YZ(:,3),3) - coor_edge_now(x_index_YZ(:,5),3) );
            S3   = L_y3.* L_z3;
            L_y4 = abs( coor_edge_pre(x_index_YZ(:,4),2) - coor_edge_now(x_index_YZ(:,5),2) );
            L_z4 = abs( coor_edge_pre(x_index_YZ(:,4),3) - coor_edge_now(x_index_YZ(:,5),3) );
            S4   = L_y4.* L_z4;
            
            temp_1 = S4./Sall;
            temp_2 = S3./Sall;
            temp_3 = S2./Sall;
            temp_4 = S1./Sall;
            topo_temp4(1:control_sum(1),1) = temp_1;
            topo_temp4(1:control_sum(1),2) = temp_2;
            topo_temp4(1:control_sum(1),3) = temp_3;
            topo_temp4(1:control_sum(1),4) = temp_4;
            % topo about sigma 根据棱边的 sigma值 调整散度算子的拓扑结构
%             temp1 = sigma_avg_pre( x_index_YZ(:,1) )./sigma_avg_now( x_index_YZ(:,5) );
%             temp2 = sigma_avg_pre( x_index_YZ(:,2) )./sigma_avg_now( x_index_YZ(:,5) );
%             temp3 = sigma_avg_pre( x_index_YZ(:,3) )./sigma_avg_now( x_index_YZ(:,5) );
%             temp4 = sigma_avg_pre( x_index_YZ(:,4) )./sigma_avg_now( x_index_YZ(:,5) );
            temp1 = 1;
            temp2 = 1;
            temp3 = 1;
            temp4 = 1;
            topo_temp4(1:control_sum(1),1) = topo_temp4(1:control_sum(1),1).*temp1;
            topo_temp4(1:control_sum(1),2) = topo_temp4(1:control_sum(1),2).*temp2;
            topo_temp4(1:control_sum(1),3) = topo_temp4(1:control_sum(1),3).*temp3;
            topo_temp4(1:control_sum(1),4) = topo_temp4(1:control_sum(1),4).*temp4;
        else
            disp('There is no ghost edge x4 in Y,Z direction');
        end
        
        if control_temp(2)~=0
            % topo about length 根据棱边的 坐标值 调整散度算子的拓扑结构
            y_index_XZ = ghost_temp((control_sum(1)+1):control_sum(2),:);
            Sall_x = abs( coor_edge_pre(y_index_XZ(:,1),1) - coor_edge_pre(y_index_XZ(:,2),1) );
            Sall_z = abs( coor_edge_pre(y_index_XZ(:,1),3) - coor_edge_pre(y_index_XZ(:,3),3) );
            Sall   = Sall_x.*Sall_z;
            L_x1 = abs( coor_edge_pre(y_index_XZ(:,1),1) - coor_edge_now(y_index_XZ(:,5),1) );
            L_z1 = abs( coor_edge_pre(y_index_XZ(:,1),3) - coor_edge_now(y_index_XZ(:,5),3) );
            S1   = L_x1.* L_z1;
            L_x2 = abs( coor_edge_pre(y_index_XZ(:,2),1) - coor_edge_now(y_index_XZ(:,5),1) );
            L_z2 = abs( coor_edge_pre(y_index_XZ(:,2),3) - coor_edge_now(y_index_XZ(:,5),3) );
            S2   = L_x2.* L_z2;
            L_x3 = abs( coor_edge_pre(y_index_XZ(:,3),1) - coor_edge_now(y_index_XZ(:,5),1) );
            L_z3 = abs( coor_edge_pre(y_index_XZ(:,3),3) - coor_edge_now(y_index_XZ(:,5),3) );
            S3   = L_x3.* L_z3;
            L_x4 = abs( coor_edge_pre(y_index_XZ(:,4),1) - coor_edge_now(y_index_XZ(:,5),1) );
            L_z4 = abs( coor_edge_pre(y_index_XZ(:,4),3) - coor_edge_now(y_index_XZ(:,5),3) );
            S4   = L_x4.* L_z4;
            
            temp_1 = S4./Sall;
            temp_2 = S3./Sall;
            temp_3 = S2./Sall;
            temp_4 = S1./Sall;
            topo_temp4((control_sum(1)+1):control_sum(2),1) = temp_1;
            topo_temp4((control_sum(1)+1):control_sum(2),2) = temp_2;
            topo_temp4((control_sum(1)+1):control_sum(2),3) = temp_3;
            topo_temp4((control_sum(1)+1):control_sum(2),4) = temp_4;
            % topo about sigma 根据棱边的 sigma值 调整散度算子的拓扑结构
%             temp1 = sigma_avg_pre( y_index_XZ(:,1) )./sigma_avg_now( y_index_XZ(:,5) );
%             temp2 = sigma_avg_pre( y_index_XZ(:,2) )./sigma_avg_now( y_index_XZ(:,5) );
%             temp3 = sigma_avg_pre( y_index_XZ(:,3) )./sigma_avg_now( y_index_XZ(:,5) );
%             temp4 = sigma_avg_pre( y_index_XZ(:,4) )./sigma_avg_now( y_index_XZ(:,5) );
            temp1 = 1;
            temp2 = 1;
            temp3 = 1;
            temp4 = 1;
            topo_temp4((control_sum(1)+1):control_sum(2),1) = topo_temp4((control_sum(1)+1):control_sum(2),1).*temp1;
            topo_temp4((control_sum(1)+1):control_sum(2),2) = topo_temp4((control_sum(1)+1):control_sum(2),2).*temp2;
            topo_temp4((control_sum(1)+1):control_sum(2),3) = topo_temp4((control_sum(1)+1):control_sum(2),3).*temp3;
            topo_temp4((control_sum(1)+1):control_sum(2),4) = topo_temp4((control_sum(1)+1):control_sum(2),4).*temp4;
        else
            disp('There is no ghost edge y4 in X,Z direction');
        end
        
        if control_temp(3)~=0
            % topo about length 根据棱边的 坐标值 调整散度算子的拓扑结构
            z_index_XY = ghost_temp((control_sum(2)+1):control_sum(3),:);
            Sall_x = abs( coor_edge_pre(z_index_XY(:,1),1) - coor_edge_pre(z_index_XY(:,2),1) );
            Sall_y = abs( coor_edge_pre(z_index_XY(:,1),2) - coor_edge_pre(z_index_XY(:,3),2) );
            Sall   = Sall_x.*Sall_y;
            L_x1 = abs( coor_edge_pre(z_index_XY(:,1),1) - coor_edge_now(z_index_XY(:,5),1) );
            L_y1 = abs( coor_edge_pre(z_index_XY(:,1),2) - coor_edge_now(z_index_XY(:,5),2) );
            S1   = L_x1.* L_y1;
            L_x2 = abs( coor_edge_pre(z_index_XY(:,2),1) - coor_edge_now(z_index_XY(:,5),1) );
            L_y2 = abs( coor_edge_pre(z_index_XY(:,2),2) - coor_edge_now(z_index_XY(:,5),2) );
            S2   = L_x2.* L_y2;
            L_x3 = abs( coor_edge_pre(z_index_XY(:,3),1) - coor_edge_now(z_index_XY(:,5),1) );
            L_y3 = abs( coor_edge_pre(z_index_XY(:,3),2) - coor_edge_now(z_index_XY(:,5),2) );
            S3   = L_x3.* L_y3;
            L_x4 = abs( coor_edge_pre(z_index_XY(:,4),1) - coor_edge_now(z_index_XY(:,5),1) );
            L_y4 = abs( coor_edge_pre(z_index_XY(:,4),2) - coor_edge_now(z_index_XY(:,5),2) );
            S4   = L_x4.* L_y4;
            
            temp_1 = S4./Sall;
            temp_2 = S3./Sall;
            temp_3 = S2./Sall;
            temp_4 = S1./Sall;
            topo_temp4((control_sum(2)+1):control_sum(3),1) = temp_1;
            topo_temp4((control_sum(2)+1):control_sum(3),2) = temp_2;
            topo_temp4((control_sum(2)+1):control_sum(3),3) = temp_3;
            topo_temp4((control_sum(2)+1):control_sum(3),4) = temp_4;
            % topo about sigma 根据棱边的 sigma值 调整散度算子的拓扑结构
%             temp1 = sigma_avg_pre( z_index_XY(:,1) )./sigma_avg_now( z_index_XY(:,5) );
%             temp2 = sigma_avg_pre( z_index_XY(:,2) )./sigma_avg_now( z_index_XY(:,5) );
%             temp3 = sigma_avg_pre( z_index_XY(:,3) )./sigma_avg_now( z_index_XY(:,5) );
%             temp4 = sigma_avg_pre( z_index_XY(:,4) )./sigma_avg_now( z_index_XY(:,5) );
            temp1 = 1;
            temp2 = 1;
            temp3 = 1;
            temp4 = 1;
            topo_temp4((control_sum(2)+1):control_sum(3),1) = topo_temp4((control_sum(2)+1):control_sum(3),1).*temp1;
            topo_temp4((control_sum(2)+1):control_sum(3),2) = topo_temp4((control_sum(2)+1):control_sum(3),2).*temp2;
            topo_temp4((control_sum(2)+1):control_sum(3),3) = topo_temp4((control_sum(2)+1):control_sum(3),3).*temp3;
            topo_temp4((control_sum(2)+1):control_sum(3),4) = topo_temp4((control_sum(2)+1):control_sum(3),4).*temp4;
        else
            disp('There is no ghost edge z4 in X,Y direction');
        end
        ghost_topo.edge4 = topo_temp4;
    else
        ghost_topo.edge4 = [];
%         warndlg('There is no ghost edge 4 !(construct_ghost_topo.m)');
    end
%% edge8 GHOST
    if ~isempty(ghost.edge8)
    else
        ghost_topo.edge8 = [];
%         warndlg('There is no ghost edge 8 !(construct_ghost_topo.m)');
    end
return