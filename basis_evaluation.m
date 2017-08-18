function [phi, dphi, phi_bordo, grad_bordo] = basis_evaluation(nod3, nod2, maps2D)
%basis_evaluation. Gives the evaluation of the P1 basis functions and their gradient
%at the 3D quadrature nodes.

% phi_func = {@(p) (-0.5*(p(1,:)+p(2,:)+p(3,:)+1));
%             @(p) (0.5*(p(1,:)+1));
%             @(p) (0.5*(p(2,:)+1));
%             @(p) (0.5*(p(3,:)+1))};

phi_func = {@(p) (-p(1,:)-p(2,:)-p(3,:)+1);
            @(p) (p(1,:));
            @(p) (p(2,:));
            @(p) (p(3,:))};
        
grad_phi_func = [-1 1 0 0;
                 -1 0 1 0;
                 -1 0 0 1];
             
nq3 = size(nod3, 2);
nq2 = size(nod2, 2);

%three derivatives, nq3 quadrature points and four functions
dphi = zeros(3,nq3,4);
phi = zeros(nq3,4);

%evaluation on 3 derivatives, nq2 quadrature boundary points, 4 faces, 4 functions.
phi_bordo = zeros(nq2,4,4);
grad_bordo = zeros(3,nq2,4,4);

for f = 1:4
    phi(:,f) = phi_func{f}(nod3)';
    dphi(:,:,f) = repmat(grad_phi_func(:,f), [1 nq3]);
    for e = 1:4
        grad_bordo(:,:,e,f) = repmat(grad_phi_func(:,f), [1 nq2]);
        phi_bordo(:,e,f) = phi_func{f}(maps2D(:,:,e)*nod2);
    end
end

end
