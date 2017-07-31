function [phi, dphi] = basis_evaluation(nod3)
%basis_evaluation. Gives the evaluation of the PI basis functions and their gradient
%at the 3D quadrature nodes.

phi_func = {@(p) (-0.5*(p(1,:)+p(2,:)+p(3,:)+1));
            @(p) (0.5*(p(1,:)+1));
            @(p) (0.5*(p(2,:)+1));
            @(p) (0.5*(p(3,:)+1))};
        
grad_phi_func = [-0.5 0.5  0   0;
                 -0.5  0  0.5  0;
                 -0.5  0   0  0.5];

%three derivatives, five quadrature points and four functions
dphi = zeros(3,5,4);
phi = zeros(5,4);

for q = 1:4
    phi(:,q) = phi_func{q}(nod3)';
    dphi(:,:,q) = repmat(grad_phi_func(:,q), [1 5]);
end

end
