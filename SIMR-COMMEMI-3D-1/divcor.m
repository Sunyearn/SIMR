function [e,isfailed, residual, iter, phi] = divcor(e,Dts,DtstG,G,...
    idxni,idxnb,phi0)
% a simple (iterative) function to make the divergence correction for 3D
% fwd responses. 
% 
% essentially we assume that there's no electric source in the model domain
% except the 'source' part, hence 
% 
%   div(E) => 0 % where there is no source 
%
% however, when a E field is calculated numerically, it is almost always
% not the case, especially for the longer periods:
%
%   div(E) = rho/epsilon
%
% hence what we want to do is to substract this part of the "non-free
% divergence electric field" for the calculated field. Apparently, the
% electric potential related to this part of field can be solved by:
%
%   div(grad(phi)) = rho/epsilon (Poission's equation)
% 
% so once we get the phi, we can simply correct the field as:
% 
%   Ec = E - grad(phi);
%   

% the real calculation (solving the Poission's equation) still utilize the 
% 'solver' functional for now, using PCG with modified IC preconditioner.
b = Dts*e;
% b2 = -Dts *e;
if nargin ==7
    b = b-phi0;
end
% calls the 'solver' function again 
phi=zeros(size(idxni));%%%%%%
[phi,isfailed,residual,iter] = solver(DtstG,b,phi,'pcg',100,1e-8);%%%%%
phiAll = zeros(length(idxni)+length(idxnb),1);
phiAll(idxni) = phi;
e = e-G*phiAll;
phi = Dts*e;
return