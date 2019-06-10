function [ex,ey,ez] = replace_ghost_edge_value(Ex,Ey,Ez,ex,ey,ez,sigma_avg_pre,sigma_avg_now, ...
                                            Nedge_pre,Nedge_now,coor_edge_pre,coor_edge_now,ghost,ghost_control,ratio)
    % 构造 ghost 的拓扑值
    [ghost_topo] = construct_ghost_topo(sigma_avg_pre,sigma_avg_now, ...
                                                                 coor_edge_pre,coor_edge_now,ghost,ghost_control);
%% edge 2
if ~isempty(ghost.edge2)
    control = ghost_control.e2;
    sum_control = cumsum(control);
    % prepare index and topo
    ghost_index_xY = ghost.edge2( 1:sum_control(1),:);
    ghost_index_xZ = ghost.edge2((sum_control(1)+1):sum_control(2),:);
    ghost_index_yX = ghost.edge2((sum_control(2)+1):sum_control(3),:);
    ghost_index_yZ = ghost.edge2((sum_control(3)+1):sum_control(4),:);
    ghost_index_zX = ghost.edge2((sum_control(4)+1):sum_control(5),:);
    ghost_index_zY = ghost.edge2((sum_control(5)+1):sum_control(6),:);
    
    ghost_topo_xY = ghost_topo.edge2( 1:sum_control(1),:);
    ghost_topo_xZ = ghost_topo.edge2((sum_control(1)+1):sum_control(2),:);
    ghost_topo_yX = ghost_topo.edge2((sum_control(2)+1):sum_control(3),:);
    ghost_topo_yZ = ghost_topo.edge2((sum_control(3)+1):sum_control(4),:);
    ghost_topo_zX = ghost_topo.edge2((sum_control(4)+1):sum_control(5),:);
    ghost_topo_zY = ghost_topo.edge2((sum_control(5)+1):sum_control(6),:);
    
    % 更新矩阵拓扑结构
    % ghost edge_x in Y direction
    for ie = 1:2 % 与大网格棱边有关的循环
        temp_ghost = ghost_index_xY(:,ie);
        temp_topo = ghost_topo_xY(:,ie);
        [value,index] = sort(temp_ghost);
        temp_topo = temp_topo(index);
        num = (ratio(2)-1)*ratio(1);
        reshape_value = reshape(value,num,length(value)/num);
        reshape_value = reshape_value';
        reshape_index = reshape(index,num,length(value)/num);
        reshape_index = reshape_index';
        reshape_topo = reshape(temp_topo,num,length(value)/num);
        reshape_topo = reshape_topo';
        for iy = 1:(ratio(2)-1)
            for ix = 1:ratio(1)
                reshape_small_value = ghost_index_xY(:,3);
                col_now = ix+(iy-1)*ratio(1);
                
                reshape_small_value = reshape_small_value( reshape_index(:,col_now) );
                ex(reshape_small_value ) = ex(reshape_small_value ) + ...
                                           Ex(reshape_value(:,col_now) ).* reshape_topo(:,col_now);
            end
        end
    end
    % ghost edge_x in Z direction
    for ie = 1:2 % 与大网格棱边有关的循环
        temp_ghost = ghost_index_xZ(:,ie);
        temp_topo = ghost_topo_xZ(:,ie);
        [value,index] = sort(temp_ghost);
        temp_topo = temp_topo(index);
        num = (ratio(3)-1)*ratio(1);
        reshape_value = reshape(value,num,length(value)/num);
        reshape_value = reshape_value';
        reshape_index = reshape(index,num,length(value)/num);
        reshape_index = reshape_index';
        reshape_topo = reshape(temp_topo,num,length(value)/num);
        reshape_topo = reshape_topo';
        for iz = 1:(ratio(3)-1)
            for ix = 1:ratio(1)
                reshape_small_value = ghost_index_xZ(:,3);
                col_now = ix+(iz-1)*ratio(1);
                
                reshape_small_value = reshape_small_value( reshape_index(:,col_now) );
                ex(reshape_small_value ) = ex(reshape_small_value ) + ...
                                           Ex(reshape_value(:,col_now) ).* reshape_topo(:,col_now);
            end
        end
    end
    % ghost edge_y in X direction
    for ie = 1:2 % 与大网格棱边有关的循环
        temp_ghost = ghost_index_yX(:,ie);
        temp_topo = ghost_topo_yX(:,ie);
        [value,index] = sort(temp_ghost);
        temp_topo = temp_topo(index);
        num = (ratio(1)-1)*ratio(2);
        reshape_value = reshape(value,num,length(value)/num);
        reshape_value = reshape_value';
        reshape_index = reshape(index,num,length(value)/num);
        reshape_index = reshape_index';
        reshape_topo = reshape(temp_topo,num,length(value)/num);
        reshape_topo = reshape_topo';
        minus_pre = Nedge_pre(1);
        minus_now = Nedge_now(1);
        for iy = 1:ratio(2)
            for ix = 1:(ratio(1)-1)
                reshape_small_value = ghost_index_yX(:,3);
                col_now = ix+(iy-1)*(ratio(1)-1);
                
                reshape_small_value = reshape_small_value( reshape_index(:,col_now) );
                ey(reshape_small_value-minus_now ) = ey(reshape_small_value-minus_now )+ ...
                                                     Ey(reshape_value(:,col_now)-minus_pre ).* reshape_topo(:,col_now);
            end
        end
    end
    % ghost edge_y in Z direction
    for ie = 1:2 % 与大网格棱边有关的循环
        temp_ghost = ghost_index_yZ(:,ie);
        temp_topo = ghost_topo_yZ(:,ie);
        [value,index] = sort(temp_ghost);
        temp_topo = temp_topo(index);
        num = (ratio(3)-1)*ratio(2);
        reshape_value = reshape(value,num,length(value)/num);
        reshape_value = reshape_value';
        reshape_index = reshape(index,num,length(value)/num);
        reshape_index = reshape_index';
        reshape_topo = reshape(temp_topo,num,length(value)/num);
        reshape_topo = reshape_topo';
        minus_pre = Nedge_pre(1);
        minus_now = Nedge_now(1);
        for iz = 1:(ratio(3)-1)
            for iy = 1:ratio(2)
                reshape_small_value = ghost_index_yZ(:,3);
                col_now = iy+(iz-1)*ratio(2);
                
                reshape_small_value = reshape_small_value( reshape_index(:,col_now) );
                ey(reshape_small_value-minus_now ) = ey(reshape_small_value-minus_now )+ ...
                                                     Ey(reshape_value(:,col_now)-minus_pre ).* reshape_topo(:,col_now);
            end
        end
    end
    % ghost edge_z in X direction
    for ie = 1:2 % 与大网格棱边有关的循环
        temp_ghost = ghost_index_zX(:,ie);
        temp_topo = ghost_topo_zX(:,ie);
        [value,index] = sort(temp_ghost);
        temp_topo = temp_topo(index);
        num = (ratio(1)-1)*ratio(3);
        reshape_value = reshape(value,num,length(value)/num);
        reshape_value = reshape_value';
        reshape_index = reshape(index,num,length(value)/num);
        reshape_index = reshape_index';
        reshape_topo = reshape(temp_topo,num,length(value)/num);
        reshape_topo = reshape_topo';
        minus_pre = Nedge_pre(1)+Nedge_pre(2);
        minus_now = Nedge_now(1)+Nedge_now(2);
        for iz = 1:ratio(3)
            for ix = 1:(ratio(1)-1)
                reshape_small_value = ghost_index_zX(:,3);
                col_now = ix+(iz-1)*(ratio(1)-1);
                
                reshape_small_value = reshape_small_value( reshape_index(:,col_now) );
                ez(reshape_small_value-minus_now ) = ez(reshape_small_value-minus_now ) + ...
                                                     Ez(reshape_value(:,col_now)-minus_pre ).* reshape_topo(:,col_now);
            end
        end
    end
    % ghost edge_z in Y direction
    for ie = 1:2 % 与大网格棱边有关的循环
        temp_ghost = ghost_index_zY(:,ie);
        temp_topo = ghost_topo_zY(:,ie);
        [value,index] = sort(temp_ghost);
        temp_topo = temp_topo(index);
        num = (ratio(2)-1)*ratio(3);
        reshape_value = reshape(value,num,length(value)/num);
        reshape_value = reshape_value';
        reshape_index = reshape(index,num,length(value)/num);
        reshape_index = reshape_index';
        reshape_topo = reshape(temp_topo,num,length(value)/num);
        reshape_topo = reshape_topo';
        minus_pre = Nedge_pre(1)+Nedge_pre(2);
        minus_now = Nedge_now(1)+Nedge_now(2);
        for iz = 1:ratio(3)
            for iy = 1:(ratio(2)-1)
                reshape_small_value = ghost_index_zY(:,3);
                col_now = iy+(iz-1)*(ratio(2)-1);
                
                reshape_small_value = reshape_small_value( reshape_index(:,col_now) );
                ez(reshape_small_value-minus_now ) = ez(reshape_small_value-minus_now ) + ...
                                                     Ez(reshape_value(:,col_now)-minus_pre ).* reshape_topo(:,col_now);
            end
        end
    end
else
end
%% edge 4
if ~isempty(ghost.edge4)
    control = ghost_control.e4;
    sum_control = cumsum(control);
    % prepare index and topo
    ghost_index_xYZ = ghost.edge4( 1:sum_control(1),:);
    ghost_index_yXZ = ghost.edge4((sum_control(1)+1):sum_control(2),:);
    ghost_index_zXY = ghost.edge4((sum_control(2)+1):sum_control(3),:);

    ghost_topo_xYZ = ghost_topo.edge4( 1:sum_control(1),:);
    ghost_topo_yXZ = ghost_topo.edge4((sum_control(1)+1):sum_control(2),:);
    ghost_topo_zXY = ghost_topo.edge4((sum_control(2)+1):sum_control(3),:);
    % ghost edge_x in YZ direction
    for ie = 1:4 % 与大网格棱边有关的循环
        temp_ghost = ghost_index_xYZ(:,ie);
        temp_topo = ghost_topo_xYZ(:,ie);
        [value,index] = sort(temp_ghost);
        temp_topo = temp_topo(index);
        num = (ratio(3)-1)*(ratio(2)-1)*ratio(1);
        reshape_value = reshape(value,num,length(value)/num);
        reshape_value = reshape_value';
        reshape_index = reshape(index,num,length(value)/num);
        reshape_index = reshape_index';
        reshape_topo = reshape(temp_topo,num,length(value)/num);
        reshape_topo = reshape_topo';
        for iz = 1:(ratio(3)-1)
            for iy = 1:(ratio(2)-1)
                for ix = 1:ratio(1)
                reshape_small_value =ghost_index_xYZ(:,5);
                col_now = ix+(iy-1)*ratio(1)+(iz-1)*ratio(1)*(ratio(2)-1);
                reshape_small_value =reshape_small_value(reshape_index(:,col_now) );
                
                ex(reshape_small_value ) = ex(reshape_small_value ) + ...
                                           Ex(reshape_value(:,col_now) ).* reshape_topo(:,col_now);
                end
            end
        end
    end
    % ghost edge_y in XZ direction
    for ie = 1:4 % 与大网格棱边有关的循环
        temp_ghost = ghost_index_yXZ(:,ie);
        temp_topo = ghost_topo_yXZ(:,ie);
        [value,index] = sort(temp_ghost);
        temp_topo = temp_topo(index);
        num = (ratio(3)-1)*(ratio(1)-1)*ratio(2);
        reshape_value = reshape(value,num,length(value)/num);
        reshape_value = reshape_value';
        reshape_index = reshape(index,num,length(value)/num);
        reshape_index = reshape_index';
        reshape_topo = reshape(temp_topo,num,length(value)/num);
        reshape_topo = reshape_topo';
        minus_pre = Nedge_pre(1);
        minus_now = Nedge_now(1);
        for iz = 1:(ratio(3)-1)
            for iy = 1:ratio(2)
                for ix = 1:(ratio(1)-1)
                reshape_small_value =ghost_index_yXZ(:,5);
                col_now = ix+(iy-1)*(ratio(1)-1)+(iz-1)*ratio(2)*(ratio(1)-1);
                reshape_small_value =reshape_small_value(reshape_index(:,col_now) );
                
                ey(reshape_small_value-minus_now ) = ey(reshape_small_value-minus_now ) + ...
                                                     Ey(reshape_value(:,col_now)-minus_pre ).* reshape_topo(:,col_now);
                end
            end
        end
    end
    % ghost edge_z in XY direction
    for ie = 1:4 % 与大网格棱边有关的循环
        temp_ghost = ghost_index_zXY(:,ie);
        temp_topo = ghost_topo_zXY(:,ie);
        [value,index] = sort(temp_ghost);
        temp_topo = temp_topo(index);
        num = (ratio(2)-1)*(ratio(1)-1)*ratio(3);
        reshape_value = reshape(value,num,length(value)/num);
        reshape_value = reshape_value';
        reshape_index = reshape(index,num,length(value)/num);
        reshape_index = reshape_index';
        reshape_topo = reshape(temp_topo,num,length(value)/num);
        reshape_topo = reshape_topo';
        minus_pre = Nedge_pre(1)+Nedge_pre(2);
        minus_now = Nedge_now(1)+Nedge_now(2);
        for iz = 1:ratio(3)
            for iy = 1:(ratio(2)-1)
                for ix = 1:(ratio(1)-1)
                reshape_small_value =ghost_index_zXY(:,5);
                col_now = ix+(iy-1)*(ratio(1)-1)+(iz-1)*(ratio(2)-1)*(ratio(1)-1);
                reshape_small_value =reshape_small_value(reshape_index(:,col_now) );
                
                ez(reshape_small_value-minus_now ) = ez(reshape_small_value-minus_now ) + ...
                                                     Ez(reshape_value(:,col_now)-minus_pre ).* reshape_topo(:,col_now);
                end
            end
        end
    end
else
end
%% edge 8
if ~isempty(ghost.edge8)
    control = ghost_control.e8;
    sum_control = cumsum(control);
    % prepare index and topo
    ghost_index_xYZ = ghost.edge8( 1:sum_control(1),:);
    ghost_index_yXZ = ghost.edge8((sum_control(1)+1):sum_control(2),:);
    ghost_index_zXY = ghost.edge8((sum_control(2)+1):sum_control(3),:);

    ghost_topo_xYZ = ghost_topo.edge8( 1:sum_control(1),:);
    ghost_topo_yXZ = ghost_topo.edge8((sum_control(1)+1):sum_control(2),:);
    ghost_topo_zXY = ghost_topo.edge8((sum_control(2)+1):sum_control(3),:);
else
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sparse函数中，如果