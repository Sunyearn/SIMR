function [points] = cal_rho(points,ilevel,iblock,per)
% calculate the apparent resistivity
    mju0=4.0*pi*1E-7; % unit is Wb/(A·m) in SI. 
    for iper = 1:length(per)
        
        omega=2*pi./per(iper); %per是频率总表
        
        Zxx = points(ilevel,iblock).Zxx(iper).matrix;
        Zxy = points(ilevel,iblock).Zxy(iper).matrix;
        Zyx = points(ilevel,iblock).Zyx(iper).matrix;
        Zyy = points(ilevel,iblock).Zyy(iper).matrix;
        
        rho_xx = Zxx.*conj(Zxx)/(mju0*omega);
        rho_xy = Zxy.*conj(Zxy)/(mju0*omega);
        rho_yx = Zyx.*conj(Zyx)/(mju0*omega);
        rho_yy = Zyy.*conj(Zyy)/(mju0*omega);
%%
        theta_xx = atan2(imag(Zxx),real(Zxx))*180/pi;
        theta_xy = atan2(imag(Zxy),real(Zxy))*180/pi;
        theta_yx = atan2(imag(Zyx),real(Zyx))*180/pi;
        theta_yy = atan2(imag(Zyy),real(Zyy))*180/pi;
        
        points(ilevel,iblock).rho_xx(iper).matrix = rho_xx;
        points(ilevel,iblock).rho_xy(iper).matrix = rho_xy;
        points(ilevel,iblock).rho_yx(iper).matrix = rho_yx;
        points(ilevel,iblock).rho_yy(iper).matrix = rho_yy;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        points(ilevel,iblock).theta_xx(iper).matrix = theta_xx;
        points(ilevel,iblock).theta_xy(iper).matrix = theta_xy;
        points(ilevel,iblock).theta_yx(iper).matrix = theta_yx;
        points(ilevel,iblock).theta_yy(iper).matrix = theta_yy;
    end
return