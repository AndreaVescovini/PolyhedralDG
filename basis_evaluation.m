function [phi, dphi, phi_bordo] = basis_evaluation(nod3, nod2, maps2D)
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

%three derivatives, five quadrature points and four functions
dphi = zeros(3,5,4);
phi = zeros(5,4);

%evaluation on 4 faces, 4 quadrature boundary points, 4 functions.
phi_bordo = zeros(4,4,4);

for f = 1:4
    phi(:,f) = phi_func{f}(nod3)';
    dphi(:,:,f) = repmat(grad_phi_func(:,f), [1 5]);
    for e = 1:4
        phi_bordo(e,:,f) = phi_func{f}(maps2D(:,:,e)*nod2);
    end
end

end
