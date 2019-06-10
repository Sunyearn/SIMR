function [model,block,per] = creat_model_commemi_3d_1()
% The struct "model" contains the parameters of coarse-grid model and
% fine-grid models. 
% The struct "block" contains the parameters of coarse-grid block and
% fine-grid blocks.
% vetror "per" is the vector of period

% Each model, coarse-grid model or fine-grid models, covers the whole domain
% of computational space.
% The coarse-grid block covers the whole domain. The two fine-grid blocks respectively covers two 
% areas that could be smaller than the coarse-grid block.
% In this test, the coarse-grid block is equal to the coarse-grid model,
% the two fine-grid blocks  may only take up part of domain of fine-grid model.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% Preparation %%%%%%%%%%%%%%%%%%%%%%%%%%
    air = 50 * 2.^(1:12);
    Nza = length(air);
    air = fliplr(air);
    air_max_depth = sum(air);
    
    Nright = 30; % the number of extended grid in positive direction
    Nmiddle = 60; % the number of centre grid
    right = 100 * 1.2.^(1:Nright); % steps of extended grid in positve direction
    left = fliplr(right);
    middle = 100*ones(1,Nmiddle); % centre grid steps 
% to reasonably divide the grid, there are four different parts of grid in
% the earth. The grid number of four parts are N1, N2, N3 and N4.
    N1 = 4;
    temp_1 = 50 * 1.1.^(1:N1);
    N2 = 8;
    temp_2 = (log10(2000)-log10(250))/N2;
    N3 = 10;
    temp_3 = (log10(20000)-log10(2000))/N3;
    N4 = 10;
    temp_4 = (log10(100000)-log10(20000))/N4;
    z1 = cumsum(temp_1);
    z1(end) = 250;
    z2 = 10.^(log10(250):temp_2:log10(2000));
    z3 = 10.^(log10(2000):temp_3:log10(20000));
    z4 = 10.^(log10(20000):temp_4:log10(100000));
    medium = [z1,z2(2:end),z3(2:end),z4(2:end)];
    medium = medium - [0,medium(1:(end-1))]; % grid steps in depth
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% Setup the grid steps and coordinates of node in x,y,z direction %%%%%%%%%%%%%%%%%%%%%%%%%%
    Dx = [left,middle,right]'; % grid steps in X direction
    Dy = [left,middle,right]'; % grid steps in Y direction
    Dz = [air,medium]'; % % grid steps in Z direction
    
    Cx = cumsum([0;Dx]); % coordinates of grid node in X direction
    % Set the centre node as origin of coordinate in X direction
    NCx = length(Cx);
    temp_n = floor(NCx/2);
    Cx = Cx - Cx(temp_n+1);
    Cy = cumsum([0;Dy]); % coordinates of grid node in Y direction
    % Set the centre node as origin of coordinate in Y direction
    NCy = length(Cy);
    temp_n = floor(NCy/2);
    Cy = Cy - Cy(temp_n+1);
    Cz = cumsum([air_max_depth;-Dz]); % coordinates of grid node in Z direction
    
    Nx=size(Dx,1);
    Ny=size(Dy,1);
    Nz=size(Dz,1);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    per = [10]; % Setup a frequency sequence
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Setup a frequency sequence
%     tt1 = 1;
%     tt2 = 1000;
%     dt = (log10(tt2) - log10(tt1))/12;
%     per = 10.^(log10(tt1):dt:log10(tt2));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%  Construct the basic resistivity model  %%%%%%%%%%%%%%%%%%%%%%%
% res is the resistivity model of inital fine-grid which covers the whole domain
    Ra = 1e10; % resistivity value of air layers
    res = Ra * ones(length(Dx),length(Dy),length(Dz));
    res(:,:,(Nza+1):end) = 100;
    res((30+20+1):(30+20+20),(30+25+1):(30+25+10),(Nza+N1+1):(Nza+N1+N2)) = 0.5;
    
    [Ly,Lx,Lz] = meshgrid(Dy,Dx,Dz);
    V = Ly.*Lx.*Lz;
    V_res = V.* res;
%% %%%%%%%%%%%%%%%%%%   Construct coarse-grid parameters  %%%%%%%%%%%%%%%%%%%%%%
    Dx2 = ( Dx(1:2:end) + Dx(2:2:end) ); % Coarse-grid steps in X direction
    Dy2 = ( Dy(1:2:end) + Dy(2:2:end) ); % Coarse-grid steps in Y direction
    Dz2 = ( Dz(1:2:end) + Dz(2:2:end) ); % Coarse-grid steps in Z direction
    
    % would rather use the coordinates from sequence of Cx than get Cx2 from calculation
    Cx2 = Cx(1:2:end); % Coarse-grid coordinates in X direction
    Cy2 = Cy(1:2:end); % Coarse-grid coordinates in Y direction
    Cz2 = Cz(1:2:end); % Coarse-grid coordinates in Z direction
% res2 is the resistivity model of coarse-grid that covers the whole domain
    res2 = V_res(1:2:(end-1),1:2:(end-1),1:2:(end-1)) + V_res(2:2:end,2:2:end,2:2:end)+...
      V_res(2:2:end,1:2:(end-1),1:2:(end-1)) + V_res(1:2:(end-1),2:2:end,1:2:(end-1))+V_res(1:2:(end-1),1:2:(end-1),2:2:end)+...
      V_res(2:2:end,2:2:end,1:2:(end-1)) + V_res(1:2:(end-1),2:2:end,2:2:end)+V_res(2:2:end,1:2:(end-1),2:2:end);
    temp_V = V(1:2:(end-1),1:2:(end-1),1:2:(end-1)) + V(2:2:end,2:2:end,2:2:end)+...
      V(2:2:end,1:2:(end-1),1:2:(end-1)) + V(1:2:(end-1),2:2:end,1:2:(end-1))+V(1:2:(end-1),1:2:(end-1),2:2:end)+...
      V(2:2:end,2:2:end,1:2:(end-1)) + V(1:2:(end-1),2:2:end,2:2:end)+V(2:2:end,1:2:(end-1),2:2:end);
    res2 = res2./temp_V;
%% 
model(1).Dx = Dx2; % The grid-step vector of coarse-grid model in X direction
model(1).Dy = Dy2; % The grid-step vector of coarse-grid model in Y direction
model(1).Dz = Dz2; % The grid-step vector of coarse-grid model in Z direction

model(1).Cx = Cx2; % The node-coordinate vector of coarse-grid model in X direction
model(1).Cy = Cy2; % The node-coordinate vector of coarse-grid model in Y direction
model(1).Cz = Cz2; % The node-coordinate vector of coarse-grid model in Z direction

model(1).level = [1,2]; % model(A1).level=[A1,B1]: A1 is the current level of grid, B1 is the maximum level of grid.
model(1).ratio = [1,1,1]; % model(A1).ratio = [a,b,c]: Current level of grid is A1. The [a,b,c] means:
                          % [a] *[current step in x direction] = [the coarsest grid step in x direction]
                          % [b] *[current step in y direction] = [the coarsest grid step in y direction]
                          % [c] *[current step in z direction] = [the coarsest grid step in z direction]
model(1).res = res2; % the resistivity model of 1st level of grid is res2. The 1st level of grid is the coarsest grid
model(1).Nza=find(model(1).res(1,1,:)<Ra,1)-1; % assume the suface is horizontal, then find the number of air layers in 1st lever of grid.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = length(Dx2); % the grid number of first-level grid in x direction
Ny = length(Dy2); % the grid number of first-level grid in y direction
Nz = length(Dz2); % the grid number of first-level grid in z direction

Nedge_x = Nx*(Ny+1)*(Nz+1); % the edge number of first-level grid in x direction
Nedge_y = (Nx+1)*Ny*(Nz+1); % the edge number of first-level grid in y direction
Nedge_z = (Nx+1)*(Ny+1)*Nz; % the edge number of first-level grid in z direction

Nface_x = (Nx+1)*Ny*Nz; % the face number of first-level grid in x direction
Nface_y = Nx*(Ny+1)*Nz; % the face number of first-level grid in y direction
Nface_z = Nx*Ny*(Nz+1); % the face number of first-level grid in z direction

Nedge = [Nedge_x,Nedge_y,Nedge_z;...
         Nedge_x,(Nedge_x+Nedge_y),(Nedge_x+Nedge_y+Nedge_z)];
Nface = [Nface_x,Nface_y,Nface_z;...
         Nface_x,(Nface_x+Nface_y),(Nface_x+Nface_y+Nface_z)];
Nnode = (Nx+1)*(Ny+1)*(Nz+1); % the node number of first-level grid

model(1).fen.edge = Nedge; % fen means face, edge and node
model(1).fen.face = Nface;
model(1).fen.node = Nnode;

%%
model(2).Dx = Dx;
model(2).Dy = Dy;
model(2).Dz = Dz;

model(2).Cx = Cx;
model(2).Cy = Cy;
model(2).Cz = Cz;

model(2).level = [2,2]; % model(A1).level=[A1,B1]: A1 is the current level of grid, B1 is the maximum level of grid.
model(2).ratio = [2,2,2]; % model(A1).ratio = [a,b,c]: Current level of grid is A1. The [a,b,c] means:
                          % [a] *[current step in x direction] = [the coarsest grid step in x direction]
                          % [b] *[current step in y direction] = [the coarsest grid step in y direction]
                          % [c] *[current step in z direction] = [the coarsest grid step in z direction]
model(2).res = res; %This is the fine-grid resistivity model which covers the whole computational domain
model(2).Nza=find(model(2).res(1,1,:)<Ra,1)-1;  % assume the suface is horizontal, then find the number of air layers in 2nd lever of grid.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = length(Dx);
Ny = length(Dy);
Nz = length(Dz);

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

model(2).fen.edge = Nedge; % fen: face, edge, node
model(2).fen.face = Nface;
model(2).fen.node = Nnode;
%% the block of coarse-grid
block.N = 3; % the number of subgrid
block(1).flag = [1,1];% block(2).flag = [ A1, B1; A2, B2]; A1 is the level of grid, B1 is the number of fine-grid areas at A1 level
block(1).area = [0,0,0,0,0,0]; % "block(A1).area(b)": the amount of rows is equal to the amount of fine-grid areas at A1 level
                               % "block(A1).area(b) = [x1,x2,y1,y2,z1,z2]":
                               % x1 and x2 are the boundary-control coefficient in X direction at the finest-grid scale;
                               %    there are x1 grids at finest-grid scale not included in negative direction of x axis in A1 level b block
                               %    there are x2 grids at finest-grid scale not included in positive direction of x axis in A1 level b block
                               % y1 and y2 are the boundary-control coefficient in Y direction at the finest-grid scale;
                               %    there are y1 grids at finest-grid scale not included in negative direction of y axis in A1 level b block
                               %    there are y2 grids at finest-grid scale not included in positive direction of y axis in A1 level b block
                               % z1 and z2 are the boundary-control coefficient in Z direction at the finest-grid scale;
                               %    there are z1 grids at finest-grid scale not included in negative direction of z axis in A1 level b block
                               %    there are z2 grids at finest-grid scale not included in positive direction of z axis in A1 level b block
                               % to undersdant this, you need to see the following 5-11 rows.
block(1).source = [0,0]; % "block(A1).source": the amount of rows is equal to the amount of subgrid areas at A1 level
                         % "block(A1).source = [a,b,c,...,m];" : The current block is A1. The amount of 
                         % upper-level blocks used to calculate the current block is m.
                         % The index of these upper-level blocks are [a,b,c,...](all before m).
block(1).Dx(1).vec = Dx2( (block(1).area(1,1)+1) : (end-block(1).area(1,2)) ); % block(A1).Dx(B1).vec: the grid-step sequence of (B1)th block of (A1)th level grid in x direction
block(1).Dy(1).vec = Dy2( (block(1).area(1,3)+1) : (end-block(1).area(1,4)) );
block(1).Dz(1).vec = Dz2( (block(1).area(1,5)+1) : (end-block(1).area(1,6)) );

block(1).Cx(1).vec = Cx2( (block(1).area(1,1)+1) : (end-block(1).area(1,2)) ); % block(A1).Cx(B1).vec: the node-coordinate sequence of (B1)th block of (A1)th level grid in x direction
block(1).Cy(1).vec = Cy2( (block(1).area(1,3)+1) : (end-block(1).area(1,4)) );
block(1).Cz(1).vec = Cz2( (block(1).area(1,5)+1) : (end-block(1).area(1,6)) );

block(1).res(1).matrix = res2( (block(1).area(1,1)+1) : (end-block(1).area(1,2)),...
                              (block(1).area(1,3)+1) : (end-block(1).area(1,4)),...
                              (block(1).area(1,5)+1) : (end-block(1).area(1,6)) ); % block(A1).res(B1).matrix: the resistivity matrix of (B1)th block of (A1)th level grid
block(1).Nza(1)=find(abs((block(1).res(1).matrix(1,1,:)-Ra)./Ra)>1e-9,1)-1;  % block(A1).Nza(B1): the air-layer number of (B1)th block of (A1)th level grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = length(block(1).Dx(1).vec);
Ny = length(block(1).Dy(1).vec);
Nz = length(block(1).Dz(1).vec);

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

block(1).fen(1).edge = Nedge; % fen means face, edge and node
block(1).fen(1).face = Nface; % fen means face, edge and node
block(1).fen(1).node = Nnode; % fen means face, edge and node
block(1).Nlevel = model(1).level(2);
%% the first block of fine-grid
block(2).flag = [1,1; % block(2).flag = [ A1, B1; A2, B2]; A1 is the level of grid, B1 is the number of fine-grid areas at A1 level
                 2,1]; % B2 is less than A2, and it means that the second block fine-grid is inactivated ( second row means second block of fine-grid)
% if mod(length(Dy)/2,2) == 0
%     y_right = (length(Dy)/2-4);
%     y_left = y_right;
% else
%     y_right = ((length(Dy)/2-1)-4);
%     y_left = y_right;
% end
% block(2).area = [(5*2),(5*2),(5*2),y_right,0,0;
%                  (5*2),(5*2), y_left,(5*2),0,0]; 
block(2).area = [30,30,30,30,0,0;
                 30,30,30,30,0,0];
block(2).source = [1,1;
                   1,1];
block(2).Dx(1).vec = Dx( (block(2).area(1,1)+1) : (end-block(2).area(1,2)) );
block(2).Dy(1).vec = Dy( (block(2).area(1,3)+1) : (end-block(2).area(1,4)) );
block(2).Dz(1).vec = Dz( (block(2).area(1,5)+1) : (end-block(2).area(1,6)) );

block(2).Cx(1).vec = Cx( (block(2).area(1,1)+1) : (end-block(2).area(1,2)) );
block(2).Cy(1).vec = Cy( (block(2).area(1,3)+1) : (end-block(2).area(1,4)) );
block(2).Cz(1).vec = Cz( (block(2).area(1,5)+1) : (end-block(2).area(1,6)) );

block(2).res(1).matrix = res( (block(2).area(1,1)+1) : (end-block(2).area(1,2)),...
                              (block(2).area(1,3)+1) : (end-block(2).area(1,4)),...
                              (block(2).area(1,5)+1) : (end-block(2).area(1,6)) );
block(2).Nza(1)=find(abs((block(2).res(1).matrix(1,1,:)-Ra)./Ra)>1e-9,1)-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = length(block(2).Dx(1).vec);
Ny = length(block(2).Dy(1).vec);
Nz = length(block(2).Dz(1).vec);

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

block(2).fen(1).edge = Nedge;
block(2).fen(1).face = Nface;
block(2).fen(1).node = Nnode;
%% the second block of fine-grid
block(2).Dx(2).vec = Dx( (block(2).area(2,1)+1) : (end-block(2).area(2,2)) );
block(2).Dy(2).vec = Dy( (block(2).area(2,3)+1) : (end-block(2).area(2,4)) );
block(2).Dz(2).vec = Dz( (block(2).area(2,5)+1) : (end-block(2).area(2,6)) );

block(2).Cx(2).vec = Cx( (block(2).area(2,1)+1) : (end-block(2).area(2,2)) );
block(2).Cy(2).vec = Cy( (block(2).area(2,3)+1) : (end-block(2).area(2,4)) );
block(2).Cz(2).vec = Cz( (block(2).area(2,5)+1) : (end-block(2).area(2,6)) );

block(2).res(2).matrix = res( (block(2).area(2,1)+1) : (end-block(2).area(2,2)),...
                              (block(2).area(2,3)+1) : (end-block(2).area(2,4)),...
                              (block(2).area(2,5)+1) : (end-block(2).area(2,6)) );
% number of air layers
block(2).Nza(2)=find(abs((block(2).res(2).matrix(1,1,:)-Ra)./Ra)>1e-9,1)-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = length(block(2).Dx(2).vec);
Ny = length(block(2).Dy(2).vec);
Nz = length(block(2).Dz(2).vec);

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

block(2).fen(2).edge = Nedge;
block(2).fen(2).face = Nface;
block(2).fen(2).node = Nnode;
block(2).Nlevel = model(1).level(2);
end