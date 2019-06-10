function [ex,ey,ez] = replace_overlap_edge_value(Ex,Ey,Ez,ex,ey,ez,sigma_avg_pre,sigma_avg_now, ...
                                                                        Nedge_pre,Nedge_now,overlap_edge,ratio)
    
    overlap_edge_x = overlap_edge.edge_x;
    overlap_edge_y = overlap_edge.edge_y;
    overlap_edge_z = overlap_edge.edge_z;
    
    pick_sigma_avg_pre = sigma_avg_pre(overlap_edge_x(:,1));
    for ix = 1:ratio(1)
        pick_sigma_avg_now = sigma_avg_now(overlap_edge_x(:,1+ix));
%         ex(overlap_edge_x(:,1+ix) ) = Ex(overlap_edge_x(:,1) ).* (pick_sigma_avg_pre./pick_sigma_avg_now);
        ex(overlap_edge_x(:,1+ix) ) = Ex(overlap_edge_x(:,1) );
    end
    
    pick_sigma_avg_pre = sigma_avg_pre(overlap_edge_y(:,1));
    minus_pre = Nedge_pre(1);
    minus_now = Nedge_now(1);
    for iy = 1:ratio(2)
        pick_sigma_avg_now = sigma_avg_now(overlap_edge_y(:,1+iy));
%         ey(overlap_edge_y(:,1+iy)-minus_now ) = Ey(overlap_edge_y(:,1)-minus_pre ).* ...
%                                                                             (pick_sigma_avg_pre./pick_sigma_avg_now);
        ey(overlap_edge_y(:,1+iy)-minus_now ) = Ey(overlap_edge_y(:,1)-minus_pre );
    end
    
    pick_sigma_avg_pre = sigma_avg_pre(overlap_edge_z(:,1));
    minus_pre = Nedge_pre(1)+Nedge_pre(2);
    minus_now = Nedge_now(1)+Nedge_now(2);
    for iz = 1:ratio(3)
        pick_sigma_avg_now = sigma_avg_now(overlap_edge_z(:,1+iz));
%         ez(overlap_edge_z(:,1+iz)-minus_now ) = Ez(overlap_edge_z(:,1)-minus_pre ).* ...
%                                                                             (pick_sigma_avg_pre./pick_sigma_avg_now);
        ez(overlap_edge_z(:,1+iz)-minus_now ) = Ez(overlap_edge_z(:,1)-minus_pre );
    end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sparse函数中，如果