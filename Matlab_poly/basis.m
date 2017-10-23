function [phi, dphi] = basis(bb, E2P, blist, Fk, nod3)
% function [phi, dphi] = basis(bb, E2P, blist, Fk, nod3)
% This function evaluates the basis functions at the 3D quadrature nodes
% for every polyhedral element.
%
% Author: Andrea Vescovini

nq3 = size(nod3,2);
Np = size(blist,1);
Ntet = size(Fk,3);

%three derivatives, nq3 quadrature points and Np functions
dphi = zeros(3, Np, nq3, Ntet);
phi = zeros(Np, nq3, Ntet);

for ie = 1:Ntet
    pt = Fk(:,:,ie)*nod3;
    for f = 1:Np

        [valx, dvalx] = LegendreP(pt(1,:), blist(f,1), bb(1,:,E2P(ie)));
        [valy, dvaly] = LegendreP(pt(2,:), blist(f,2), bb(2,:,E2P(ie)));
        [valz, dvalz] = LegendreP(pt(3,:), blist(f,3), bb(3,:,E2P(ie)));

        phi(f,:,ie) = valx.*valy.*valz;
        dphi(1,f,:,ie) = dvalx.*valy.*valz;
        dphi(2,f,:,ie) = valx.*dvaly.*valz;
        dphi(3,f,:,ie) = valx.*valy.*dvalz;
    end

end

end
