function [phi, dphi, phi_bordo, grad_bordo] = basis(bb, blist, Fk, nod3, node_maps, nod2)
% function [phi, dphi, phi_bordo, grad_bordo] = basis(bb, blist, Fk, nod3, node_maps, nod2)
% For every element evaluate the basis functions at the 3D and 2D quadrature nodes.

nq3 = size(nod3,2);
nq2 = size(nod2,2);
Np = size(blist,1);
K = size(Fk,3);

%three derivatives, nq3 quadrature points and Np functions
dphi = zeros(3, nq3, Np, K);
phi = zeros(nq3, Np, K);

%evaluation on 3 derivatives, nq2 quadrature boundary points, 4 faces, Np functions.
grad_bordo = zeros(3, nq2, 4, Np, K);
phi_bordo = zeros(nq2, 4, Np, K);

for ie = 1:K
    
    pt = Fk(:,:,ie)*nod3;
    for f = 1:Np
    
        [valx, dvalx] = LegendreP(pt(1,:), blist(f,1), bb(1,:,ie));
        [valy, dvaly] = LegendreP(pt(2,:), blist(f,2), bb(2,:,ie));
        [valz, dvalz] = LegendreP(pt(3,:), blist(f,3), bb(3,:,ie));
%         dvalx = LegendrePder(pt(1,:), blist(f,1), bb(1,:,ie));
%         dvaly = LegendrePder(pt(2,:), blist(f,2), bb(2,:,ie));
%         dvalz = LegendrePder(pt(3,:), blist(f,3), bb(3,:,ie));
       
        phi(:,f,ie) = (valx.*valy.*valz)';
        dphi(1,:,f,ie) = dvalx.*valy.*valz;
        dphi(2,:,f,ie) = valx.*dvaly.*valz;
        dphi(3,:,f,ie) = valx.*valy.*dvalz;
    end

    for e = 1:4
        pt = Fk(:,:,ie)*node_maps(:,:,e)*nod2;
        for f = 1:Np
        [valx, dvalx] = LegendreP(pt(1,:), blist(f,1), bb(1,:,ie));
        [valy, dvaly] = LegendreP(pt(2,:), blist(f,2), bb(2,:,ie));
        [valz, dvalz] = LegendreP(pt(3,:), blist(f,3), bb(3,:,ie));
%         dvalx = LegendrePder(pt(1,:), blist(f,1), bb(1,:,ie));
%         dvaly = LegendrePder(pt(2,:), blist(f,2), bb(2,:,ie));
%         dvalz = LegendrePder(pt(3,:), blist(f,3), bb(3,:,ie));
    
            phi_bordo(:,e,f,ie) = (valx.*valy.*valz)';
            grad_bordo(1,:,e,f,ie) = dvalx.*valy.*valz;
            grad_bordo(2,:,e,f,ie) = valx.*dvaly.*valz;
            grad_bordo(3,:,e,f,ie) = valx.*valy.*dvalz;
        end
    end

end

end

