function [points] = cal_impedance(block,points,ilevel,iblock,per)
% calculate the impedence
Nlevel = block(ilevel).Nlevel;
Nblock = block(ilevel).flag(1,2);
            
Cx = block(ilevel).Cx(iblock).vec;
Cy = block(ilevel).Cy(iblock).vec;
% Cz = block(ilevel).Cz(iblock).vec((Nza+1):end);
Cz = block(ilevel).Cz(iblock).vec;
if length(Cx) ~= 1
    Cxm = ( Cx(1:(end-1)) + Cx(2:end) )/2;
else
    Cxm = Cx;
end
if length(Cy) ~= 1
    Cym = ( Cy(1:(end-1)) + Cy(2:end) )/2;
else
    Cym = Cy;
end
if length(Cz) ~= 1
    Czm = ( Cz(1:(end-1)) + Cz(2:end) )/2;
else
	Czm = Cz;
end
            
x = points(ilevel,iblock).x;
y = points(ilevel,iblock).y;
z = points(ilevel,iblock).z;
for iper = 1:length(per)
        for ipol = 1:2
            % Ex
            [My,Mx,Mz] = meshgrid(Cy,Cxm,Cz);
%             V = block(ilevel).Ex(ipol,iblock,iper).matrix(:,:,((Nza+1):end));
            V = block(ilevel).Ex(ipol,iblock,iper).matrix;
            v = interp3(My,Mx,Mz,V,y,x,z);
            points(ilevel,iblock).Ex(ipol).matrix = v;
            % Ey
            [My,Mx,Mz] = meshgrid(Cym,Cx,Cz);
%             V = block(ilevel).Ey(ipol,iblock,iper).matrix(:,:,((Nza+1):end));
            V = block(ilevel).Ey(ipol,iblock,iper).matrix;
            v = interp3(My,Mx,Mz,V,y,x,z);
            points(ilevel,iblock).Ey(ipol).matrix = v;
            % Ez
            [My,Mx,Mz] = meshgrid(Cy,Cx,Czm);
%             V = block(ilevel).Ez(ipol,iblock,iper).matrix(:,:,((Nza+1):end));
            V = block(ilevel).Ez(ipol,iblock,iper).matrix;
            v = interp3(My,Mx,Mz,V,y,x,z);
            points(ilevel,iblock).Ez(ipol).matrix = v;
            % Hx
            [My,Mx,Mz] = meshgrid(Cym,Cx,Czm);
%             V = block(ilevel).Hx(ipol,iblock,iper).matrix(:,:,((Nza+1):end));
            V = block(ilevel).Hx(ipol,iblock,iper).matrix;
            v = interp3(My,Mx,Mz,V,y,x,z);
            points(ilevel,iblock).Hx(ipol).matrix = v;
            % Hy
            [My,Mx,Mz] = meshgrid(Cy,Cxm,Czm);
%             V = block(ilevel).Hy(ipol,iblock,iper).matrix(:,:,((Nza+1):end));
            V = block(ilevel).Hy(ipol,iblock,iper).matrix;
            v = interp3(My,Mx,Mz,V,y,x,z);
            points(ilevel,iblock).Hy(ipol).matrix = v;
            % Hz
            [My,Mx,Mz] = meshgrid(Cym,Cxm,Cz);
%             V = block(ilevel).Hz(ipol,iblock,iper).matrix(:,:,((Nza+1):end));
            V = block(ilevel).Hz(ipol,iblock,iper).matrix;
            v = interp3(My,Mx,Mz,V,y,x,z);
            points(ilevel,iblock).Hz(ipol).matrix = v;
        end
        Ex1 = points(ilevel,iblock).Ex(1).matrix;
%         Ex1 = reshape(Ex1,numel(Ex1),1);
        Ex2 = points(ilevel,iblock).Ex(2).matrix;
%         Ex2 = reshape(Ex2,numel(Ex2),1);
        Ey1 = points(ilevel,iblock).Ey(1).matrix;
%         Ey1 = reshape(Ey1,numel(Ey1),1);
        Ey2 = points(ilevel,iblock).Ey(2).matrix;
%         Ey2 = reshape(Ey2,numel(Ey2),1);
        Ez1 = points(ilevel,iblock).Ez(1).matrix;
%         Ez1 = reshape(Ez1,numel(Ez1),1);
        Ez2 = points(ilevel,iblock).Ez(2).matrix;
%         Ez2 = reshape(Ez2,numel(Ez2),1);
        
        Hx1 = points(ilevel,iblock).Hx(1).matrix;
%         Hx1 = reshape(Hx1,numel(Hx1),1);
        Hx2 = points(ilevel,iblock).Hx(2).matrix;
%         Hx2 = reshape(Hx2,numel(Hx2),1);
        Hy1 = points(ilevel,iblock).Hy(1).matrix;
%         Hy1 = reshape(Hy1,numel(Hy1),1);
        Hy2 = points(ilevel,iblock).Hy(2).matrix;
%         Hy2 = reshape(Hy2,numel(Hy2),1);
        Hz1 = points(ilevel,iblock).Hz(1).matrix;
%         Hy1 = reshape(Hy1,numel(Hy1),1);
        Hz2 = points(ilevel,iblock).Hz(2).matrix;
%         Hy2 = reshape(Hy2,numel(Hy2),1);
        Zxx = ( Ex1.*Hy2-Ex2.*Hy1 )...
                ./( Hx1.*Hy2-Hx2.*Hy1 );
        Zxy = ( Ex2.*Hx1-Ex1.*Hx2 )...
                ./( Hx1.*Hy2-Hx2.*Hy1 );
        Zyx = ( Ey1.*Hy2-Ey2.*Hy1 )...
                ./( Hx1.*Hy2-Hx2.*Hy1 );
        Zyy = ( Ey2.*Hx1-Ey1.*Hx2 )...
                ./( Hx1.*Hy2-Hx2.*Hy1 );

        points(ilevel,iblock).Zxx(iper).matrix = Zxx;
        points(ilevel,iblock).Zxy(iper).matrix = Zxy;
        points(ilevel,iblock).Zyx(iper).matrix = Zyx;
        points(ilevel,iblock).Zyy(iper).matrix = Zyy;
end
points(ilevel,iblock).Nlevel = Nlevel;
points(ilevel,iblock).Nblock = Nblock;
return