function [phi_bordo, grad_bordo] = basis_boundary(bb, blist, Fk, node_maps, nod2)
% function [phi_bordo, grad_bordo] = basis(bb, blist, Fk, node_maps, nod2)
% For every element evaluate the basis functions at the 2D quadrature nodes.

nq2 = size(nod2,2);
Np = size(blist,1);
K = size(Fk,3);

%evaluation on 3 derivatives, nq2 quadrature boundary points, 4 faces, Np functions.
grad_bordo = zeros(3, Np, nq2, 4, K);
phi_bordo = zeros(Np, nq2, 4, K);

for ie = 1:K
    
    for e = 1:4
        pt = Fk(:,:,ie)*node_maps(:,:,e)*nod2;
        for f = 1:Np
        [valx, dvalx] = LegendreP(pt(1,:), blist(f,1), bb(1,:,ie));
        [valy, dvaly] = LegendreP(pt(2,:), blist(f,2), bb(2,:,ie));
        [valz, dvalz] = LegendreP(pt(3,:), blist(f,3), bb(3,:,ie));
    
            phi_bordo(f,:,e,ie) = valx.*valy.*valz;
            grad_bordo(1,f,:,e,ie) = dvalx.*valy.*valz;
            grad_bordo(2,f,:,e,ie) = valx.*dvaly.*valz;
            grad_bordo(3,f,:,e,ie) = valx.*valy.*dvalz;
        end
    end

end

end

