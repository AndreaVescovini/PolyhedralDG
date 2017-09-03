function [err_L2, err_H10] = errors(uex, uex_grad, u, x, y, z, N)
%function [err_L2, err_H10] = errors(uex, uex_grad, u, x, y, z, N)
%   Compute L2 norm and H1 seminorm of the error u-uex

Np = (N+1)*(N+2)*(N+3)/6; %number of nodes for every element
Nfaces = 4; % number of faces for every element
[r, s, t] = set_dof_lin; % dof on the reference tetrahedron

[nod2, ~, nod3, wei3, node_maps] = quadrature(Nfaces);
[phi, dphi, ~, ~] = basis_evaluation(nod3, nod2, node_maps);

err_L2 = 0;
err_H10 = 0;

for ie = 1:size(x,2)
    [J, Jcof, Jdet, trasl] = jacobians(x(:,ie),y(:,ie),z(:,ie),r,s,t);
    for q = 1:length(wei3)
        % evaluation of the dg-fem solution and its gradient at the quadrature node
        u_h = phi(q,:)*u(:,ie);
        grad_uh = 0;
        for i = 1:Np
            grad_uh = grad_uh + dphi(:,q,i)*u(i,ie);
        end
        
        point_diff = uex_grad(trasl+J*nod3(:,q)) - 1/Jdet*Jcof*grad_uh;
        err_H10 = err_H10 + wei3(q)*abs(Jdet)*(point_diff'*point_diff);
        err_L2 = err_L2 + wei3(q)*abs(Jdet)*(u_h - uex(trasl+J*nod3(:,q)))^2;
        
    end
end

err_L2 = sqrt(err_L2);
err_H10 = sqrt(err_H10);

end

