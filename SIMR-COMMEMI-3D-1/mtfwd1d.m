function F = mtfwd1d(per,res,l,pol)
%=========================checking parameters=============================%
switch nargin
    case 0
        error('not enough input arguments, 2 at least')
    case 1
        error('not enough input arguments, 2 at least')       
    case 2
        l = 100; % treat the earth as half space
        pol = 'E';
    case 3
        pol = 'E';
end
if (size(l)~=size(res))
    if length(l)==length(res)
        res=res';
    else
	    disp('please check the input parametres ')
        error('res and l should have same size');
    end
end
%=========================================================================%
mju0=4.0*pi*1E-7; % unit is Wb/(A��m) in SI. 
% Top Boundary condition: E=1+0j
omega=2*pi./per;
% now adding air layers
sigma=1./res;
Dl=(l(1:end-1)+l(2:end))/2; 
% distances of mid points of each layer, used for FD method
Nl=length(l); % number of layers
NP=length(per); % number of freqs
F=zeros(Nl+1,NP); % Nl grids with Nl+1 interfaces
%   now we set up the differential (gradient) operators
G = diag(ones(Nl,1),1)-diag(ones(Nl+1,1),0);
G(Nl+1,:) = [];
%   which is something like: 
%  -1    1   0   0
%   0   -1   1   0
%   0    0  -1   1
%   so diag(1/l)*G will become the first order differencial operator
D = diag(ones(Nl,1),0)-diag(ones(Nl-1,1),1);
D(Nl,:) = [];
%   which is something like: 
%   1   -1    0    0
%   0    1   -1    0
%   0    0    1   -1
%   so diag(1/l)*D will become the (second half of) second order 
%   differencial operator
for iper=1:NP
    switch upper(pol)
        % essentially sets up (second-order) differencial operators for 1D 
        % problems
        case 'E'
            %   calculate conductance 
            cond = l.*sigma; 
            %   and average sigma - centered differences
            sigma_avg = (cond(1:end-1)+cond(2:end))./...
                (l(1:end-1)+l(2:end));
            %   now setting up the 2nd order difference operator
            %   (which will become curl-curl in higher dimensions)
            A = diag(1./Dl)*D*diag(1./l)*G;
            %   padding zeros to the top and bottom rows for boundary 
            %   condition
            A = [zeros(1,Nl+1); A; zeros(1,Nl+1)]; 
            %   now setup the operator for the last line -- for lower BC
            %   with one-sided (backward) FD 
            ln=l(end);
            lnm1=l(end-1);
            Dln=Dl(end);
            a1=(4*Dln-lnm1)/(2*Dln*ln);
            a2=-2*Dln/(lnm1*ln);
            a3=ln/(2*Dln*lnm1);
            % c = 1/(1j*k)
            c = 1/(1i*sqrt(1j*omega(iper)*mju0*sigma(end)));
            A(end,end) = c*a1-2;
            A(end,end-1) = c*a2;
            A(end,end-2) = c*a3;      
            %   i*omega*mu*sigma
%             iwus = [ 1; (isign*1j*omega(iper)*mju0*sigma_avg); 1];
            iwus = [ 1; (1j*omega(iper)*mju0*sigma_avg); 1];
%             iwus = [ 1; (isign*1j*omega(iper)*mju0*sigma_avg + isign * mju0 * omega * omega*1); 1];
            A = A + diag(iwus,0);
        case 'B'
            %   now setting up the 2nd order difference operator
            %   (which will become curl-curl in higher dimensions)
            A = diag(1./Dl)*D*diag(1./sigma./l)*G;
            %   padding zeros to the top and bottom rows for boundary 
            %   condition
            A = [zeros(1,Nl+1); A; zeros(1,Nl+1)]; 
            %   now setup the operator for the last line -- for lower BC
            %   with one-sided (backward) FD 
            ln=l(end);
            lnm1=l(end-1);
            Dln=Dl(end);
            a1=(4*Dln-lnm1)/(2*Dln*ln);
            a2=-2*Dln/(lnm1*ln);
            a3=ln/(2*Dln*lnm1);
            % c = 1/(1j*k)
            c = 1/(1i*sqrt(1j*omega(iper)*mju0*sigma(end)));
            A(end,end) = c*a1-2;
            A(end,end-1) = c*a2;
            A(end,end-2) = c*a3;                  
            %   i*omega*mu*sigma
%             iwus = [ 1; isign*1j*omega(iper)*mju0*ones(Nl-1,1); 1];
            iwus = [ 1; 1j*omega(iper)*mju0*ones(Nl-1,1); 1];
            A = A + diag(iwus,0);
        otherwise
            error('polarization option not recognized, must be "E" or "B"')
    end
    %   direct solve with a simple (jacobian) pre-conditioner ...
    T = diag(1./diag(A),0); % T=P^(-1)
    A = T*A;
    %   Boundary - assume 1 on the top and 0 on the bottom
    b = [1; zeros(Nl,1)];
    b = T*b;
    F(:,iper) =A\b; % now direct solve the problem...
end 
return