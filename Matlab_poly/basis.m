function [phi, dphi] = basis(bb, blist, Fk, nod3)
% function [phi, dphi] = basis(bb, blist, Fk, nod3)
% For every element evaluate the basis functions at the 3D quadrature nodes.

nq3 = size(nod3,2);
Np = size(blist,1);
K = size(Fk,3);

%three derivatives, nq3 quadrature points and Np functions
dphi = zeros(3, Np, nq3, K);
phi = zeros(Np, nq3, K);

for ie = 1:K
    pt = Fk(:,:,ie)*nod3;
    for f = 1:Np
    
        [valx, dvalx] = LegendreP(pt(1,:), blist(f,1), bb(1,:,ie));
        [valy, dvaly] = LegendreP(pt(2,:), blist(f,2), bb(2,:,ie));
        [valz, dvalz] = LegendreP(pt(3,:), blist(f,3), bb(3,:,ie));
       
        phi(f,:,ie) = valx.*valy.*valz;
        dphi(1,f,:,ie) = dvalx.*valy.*valz;
        dphi(2,f,:,ie) = valx.*dvaly.*valz;
        dphi(3,f,:,ie) = valx.*valy.*dvalz;
    end

end

end

