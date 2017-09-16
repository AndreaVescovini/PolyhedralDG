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

for q = 1:nq3
    pt = Fk(:,:,ie)*nod3(:,q);
    for f = 1:Np
        valx = LegendreP(pt(1), blist(f,1), bb(1,:,ie));
        valy = LegendreP(pt(2), blist(f,2), bb(2,:,ie));
        valz = LegendreP(pt(3), blist(f,3), bb(3,:,ie));
        
        phi(q,f,ie) = valx*valy*valz;
        dphi(1,q,f,ie) = LegendrePder(pt(1), blist(f,1), bb(1,:,ie))*valy*valz;
        dphi(2,q,f,ie) = valx*LegendrePder(pt(2), blist(f,2), bb(2,:,ie))*valz;
        dphi(3,q,f,ie) = valx*valy*LegendrePder(pt(3), blist(f,3), bb(3,:,ie));
    end
end

for q = 1:nq2
    for e = 1:4
        pt = Fk(:,:,ie)*node_maps(:,:,e)*nod2(:,q);
        for f = 1:Np
            valx = LegendreP(pt(1), blist(f,1), bb(1,:,ie));
            valy = LegendreP(pt(2), blist(f,2), bb(2,:,ie));
            valz = LegendreP(pt(3), blist(f,3), bb(3,:,ie));
        
            phi_bordo(q,e,f,ie) = valx*valy*valz;
            grad_bordo(1,q,e,f,ie) = LegendrePder(pt(1), blist(f,1), bb(1,:,ie))*valy*valz;
            grad_bordo(2,q,e,f,ie) = valx*LegendrePder(pt(2), blist(f,2), bb(2,:,ie))*valz;
            grad_bordo(3,q,e,f,ie) = valx*valy*LegendrePder(pt(3), blist(f,3), bb(3,:,ie));
            
            
        end
    end
end

end

end

