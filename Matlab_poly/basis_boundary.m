function [phi_bordo, grad_bordo] = basis_boundary(bb, E2P, faces_neig, blist, Fk, node_maps, nod2)
% function [phi_bordo, grad_bordo] = basis(bb, E2P, faces_neig, blist, Fk, node_maps, nod2)
% For face evaluate the basis functions both of E1 and E2 at the 2D quadrature nodes.

nq2 = size(nod2,2);
Np = size(blist,1);
Nfaces = size(faces_neig,1);

%evaluation on 3 derivatives, nq2 quadrature boundary points, Nf faces, 2
%neihbors
grad_bordo = zeros(3, Np, nq2, Nfaces, 2);
phi_bordo = zeros(Np, nq2, Nfaces, 2);

for e = 1:Nfaces
     E1 = faces_neig(e,1);
     e_E1 = faces_neig(e,2);
     pt = Fk(:,:,E1)*node_maps(:,:,e_E1)*nod2;
     for f = 1:Np
         [valx, dvalx] = LegendreP(pt(1,:), blist(f,1), bb(1,:,E2P(E1)));
         [valy, dvaly] = LegendreP(pt(2,:), blist(f,2), bb(2,:,E2P(E1)));
         [valz, dvalz] = LegendreP(pt(3,:), blist(f,3), bb(3,:,E2P(E1)));
    
         phi_bordo(f,:,e,1) = valx.*valy.*valz;
         grad_bordo(1,f,:,e,1) = dvalx.*valy.*valz;
         grad_bordo(2,f,:,e,1) = valx.*dvaly.*valz;
         grad_bordo(3,f,:,e,1) = valx.*valy.*dvalz;
     end
     E2 = faces_neig(e,3);
     if E2 ~= 0
         e_E2 = faces_neig(e,4);
         pt = Fk(:,:,E2)*node_maps(:,:,e_E2)*nod2;
         for f = 1:Np
            [valx, dvalx] = LegendreP(pt(1,:), blist(f,1), bb(1,:,E2P(E2)));
            [valy, dvaly] = LegendreP(pt(2,:), blist(f,2), bb(2,:,E2P(E2)));
            [valz, dvalz] = LegendreP(pt(3,:), blist(f,3), bb(3,:,E2P(E2)));
    
            phi_bordo(f,:,e,2) = valx.*valy.*valz;
            grad_bordo(1,f,:,e,2) = dvalx.*valy.*valz;
            grad_bordo(2,f,:,e,2) = valx.*dvaly.*valz;
            grad_bordo(3,f,:,e,2) = valx.*valy.*dvalz;
         end
     end
end

end

