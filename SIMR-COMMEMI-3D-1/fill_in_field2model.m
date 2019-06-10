function [model,block] = fill_in_field2model(model,block,ilevel,iblo,iper,field,skip)
    Nx = length(model(ilevel).Dx);
    Ny = length(model(ilevel).Dy);
    Nz = length(model(ilevel).Dz);
    
    ratio = (model((model(1).level(2))).ratio) ./ (model(ilevel).ratio);
    ratio_change = [ratio(1),ratio(1),ratio(2),ratio(2),ratio(3),ratio(3)];
    area = block(ilevel).area(iblo,:)./ratio_change; % the grids to skip, at current scale of grid
    
    Ex = zeros(Nx,Ny+1,Nz+1);
    Ey = zeros(Nx+1,Ny,Nz+1);
    Ez = zeros(Nx+1,Ny+1,Nz);
    
    skip_x = skip(1);
    skip_y = skip(2);
    skip_z = skip(3);
    
    for ipol = 1:2 % 2 polarization mode
        if ilevel == 1
            model(ilevel).Ex(ipol,iblo,iper).matrix = field(ipol).Ex;
            model(ilevel).Ey(ipol,iblo,iper).matrix = field(ipol).Ey;
            model(ilevel).Ez(ipol,iblo,iper).matrix = field(ipol).Ez;
        block(ilevel).Ex(ipol,iblo,iper).matrix = field(ipol).Ex;
        block(ilevel).Ey(ipol,iblo,iper).matrix = field(ipol).Ey;
        block(ilevel).Ez(ipol,iblo,iper).matrix = field(ipol).Ez;
        else
            model(ilevel).Ex(ipol,iblo,iper).matrix( ((area(1)+skip_x+1):(end-skip_x-area(2))),...
                                        ((area(3)+skip_y+1):(end-skip_y-area(4))),...
                                        ((area(5)+skip_z+1):(end-skip_z-area(6))) )=field(ipol).Ex((skip_x+1):(end-skip_x),...
                                                                                                   (skip_y+1):(end-skip_y),...
                                                                                                   (skip_z+1):(end-skip_z));
            model(ilevel).Ey(ipol,iblo,iper).matrix( ((area(1)+skip_x+1):(end-skip_x-area(2))),...
                                        ((area(3)+skip_y+1):(end-skip_y-area(4))),...
                                        ((area(5)+skip_z+1):(end-skip_z-area(6))) )=field(ipol).Ey((skip_x+1):(end-skip_x),...
                                                                                                   (skip_y+1):(end-skip_y),...
                                                                                                   (skip_z+1):(end-skip_z));
            model(ilevel).Ez(ipol,iblo,iper).matrix( ((area(1)+skip_x+1):(end-skip_x-area(2))),...
                                        ((area(3)+skip_y+1):(end-skip_y-area(4))),...
                                        ((area(5)+skip_z+1):(end-skip_z-area(6))) )=field(ipol).Ez((skip_x+1):(end-skip_x),...
                                                                                                   (skip_y+1):(end-skip_y),...
                                                                                                   (skip_z+1):(end-skip_z));
            block(ilevel).Ex(ipol,iblo,iper).matrix = field(ipol).Ex;
            block(ilevel).Ey(ipol,iblo,iper).matrix = field(ipol).Ey;
            block(ilevel).Ez(ipol,iblo,iper).matrix = field(ipol).Ez;
        end
        if ilevel == 1
            model(ilevel).Hx(ipol,iblo,iper).matrix = field(ipol).Hx;
            model(ilevel).Hy(ipol,iblo,iper).matrix = field(ipol).Hy;
            model(ilevel).Hz(ipol,iblo,iper).matrix = field(ipol).Hz;
            block(ilevel).Hx(ipol,iblo,iper).matrix = field(ipol).Hx;
            block(ilevel).Hy(ipol,iblo,iper).matrix = field(ipol).Hy;
            block(ilevel).Hz(ipol,iblo,iper).matrix = field(ipol).Hz;
        else
            block(ilevel).Hx(ipol,iblo,iper).matrix = field(ipol).Hx;
            block(ilevel).Hy(ipol,iblo,iper).matrix = field(ipol).Hy;
            block(ilevel).Hz(ipol,iblo,iper).matrix = field(ipol).Hz;
        end
        
    end
return