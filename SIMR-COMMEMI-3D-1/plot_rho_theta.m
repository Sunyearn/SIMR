function plot_rho_theta(field3D,ilevel,iblock,iper,per,mode,direction)
%   a simple script to plot the field3D data structure (which contains) 
if strcmp(mode,'LINE')
% figure('Units', 'normalized', 'Position', [0.3 0.05 0.3 0.8])
Cx=field3D(ilevel,iblock).x;
Cy=field3D(ilevel,iblock).y;
Cz=field3D(ilevel,iblock).z;
% Cz=field3D(ilevel,iblock).z;
Cz_rho_xy=field3D(ilevel,iblock).rho_xy(iper).matrix;
Cz_rho_xx=field3D(ilevel,iblock).rho_xx(iper).matrix;
Cz_theta_xy=field3D(ilevel,iblock).theta_xy(iper).matrix;
Cz_theta_xx=field3D(ilevel,iblock).theta_xx(iper).matrix;
Cz_rho_yx=field3D(ilevel,iblock).rho_yx(iper).matrix;
Cz_rho_yy=field3D(ilevel,iblock).rho_yy(iper).matrix;
Cz_theta_yx=field3D(ilevel,iblock).theta_yx(iper).matrix;
Cz_theta_yy=field3D(ilevel,iblock).theta_yy(iper).matrix;

figure;
h1 = subplot(2,1,1);
h2 = subplot(2,1,2);
period = per(iper);
%    case 'rho'  
if strcmp(direction,'X')
    semilogy(h1,Cx,Cz_rho_xy,'r+'); 
    hold(h1,'on')
    semilogy(h1,Cx,Cz_rho_yx,'bx');
    xlabel(h1,'Northing(km)');
    title(h1,['Line X Resistivity ' num2str(period) ' s'])
elseif strcmp(direction,'Y')
    semilogy(h1,Cy,Cz_rho_xy,'r+'); 
    hold(h1,'on')
    semilogy(h1,Cy,Cz_rho_yx,'bx');
    xlabel(h1,'Easting(km)');
    title(h1,['Line Y Resistivity ' num2str(period) ' s'])
else
    errordlg('wrong direction mode');
end
legend(h1,'Rho-XY','Rho-YX');
ylabel(h1,'Res(Ohm*m)');
%%
%    case 'theta'
if min([min(Cz_theta_xy),min(Cz_theta_yx)]) >=0
floor_y = 0;
elseif min([min(Cz_theta_xy),min(Cz_theta_yx)])<0 && min([min(Cz_theta_xy),min(Cz_theta_yx)])>= -90
floor_y = -90;
elseif min([min(Cz_theta_xy),min(Cz_theta_yx)])<-90 && min([min(Cz_theta_xy),min(Cz_theta_yx)])>= -180
floor_y = -180;
end

if max([max(Cz_theta_xy),max(Cz_theta_yx)]) >= 90 
ceil_y = 180;
elseif max([max(Cz_theta_xy),max(Cz_theta_yx)])<90 && max([max(Cz_theta_xy),max(Cz_theta_yx)])>= 0
ceil_y = 90;
elseif max([max(Cz_theta_xy),max(Cz_theta_yx)])<0
ceil_y = 0;
end
step_y = (ceil_y-floor_y)/6;

floor_y = -180;
ceil_y = 180;
step_y = (ceil_y-floor_y)/8;
        
if strcmp(direction,'X')
    plot(h2,Cx,Cz_theta_xy,'r+'); 
    hold(h2,'on')
    plot(h2,Cx,Cz_theta_yx,'bx');
    xlabel(h2,'Northing(km)');
    title(h2,['Line X Theta ' num2str(period) ' s'])
elseif strcmp(direction,'Y')
    plot(h2,Cy,Cz_theta_xy,'r+'); 
    hold(h2,'on')
    plot(h2,Cy,Cz_theta_yx,'bx');
    xlabel(h2,'Easting(km)');
    title(h2,['Line Y Theta ' num2str(period) ' s'])
else
    errordlg('wrong direction mode');
end
legend(h2,'Theta-XY','Theta-YX');
ylabel(h2,'Theta');
set(h2,'YLim',[floor_y ceil_y]);
set(h2,'YTick',[floor_y:step_y:ceil_y]);
set(h2,'YTickLabel',[floor_y:step_y:ceil_y]);
elseif strcmp(mode,'POINT')
end
return















