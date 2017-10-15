function [err_L2, err_H10] = errors(uex, uex_grad, u, x, y, z, N)
%function [err_L2, err_H10] = errors(uex, uex_grad, u, x, y, z, N)
%   Compute L2 norm and H1 seminorm of the error u-uex

Np = (N+1)*(N+2)*(N+3)/6; %number of nodes for every element
Nfaces = 4; % number of faces for every element
[r, s, t] = set_dof(N); % dof on the reference tetrahedron

[nod2, ~, nod3, wei3, node_maps, ~] = quadrature(Nfaces);
[phi, dphi, ~, ~, ~] = basis_evaluation(nod3, nod2, node_maps, N, Np);

err_L2 = 0;
err_H10 = 0;

[trasl, J, Jcof, Jdet] = jacobians(x(1:4,:),y(1:4,:),z(1:4,:),r(1:4),s(1:4),t(1:4));

for ie = 1:size(x,2)
    
    for q = 1:length(wei3)
        % evaluation of the dg-fem solution and its gradient at the quadrature node
        u_h = phi(:,q)'*u(:,ie);
        grad_uh = zeros(3,1);
        for ii = 1:3
            grad_uh(ii) = dphi(ii,:,q)*u(:,ie);
        end
        
        point_diff = uex_grad(trasl(:,ie)+J(:,:,ie)*nod3(:,q)) - Jcof(:,:,ie)*grad_uh./Jdet(ie);
        err_H10 = err_H10 + wei3(q)*abs(Jdet(ie))*(point_diff'*point_diff);
        err_L2 = err_L2 + wei3(q)*abs(Jdet(ie))*(u_h - uex(trasl(:,ie)+J(:,:,ie)*nod3(:,q)))^2;
        
    end
end

err_L2 = sqrt(err_L2);
err_H10 = sqrt(err_H10);

end

