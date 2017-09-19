function [err_L2, err_H10] = errors(uex, uex_grad, u, x, y, z, N)
%function [err_L2, err_H10] = errors(uex, uex_grad, u, x, y, z, N)
%   Compute L2 norm and H1 seminorm of the error u-uex

Np = (N+1)*(N+2)*(N+3)/6; %number of nodes for every element

[~, ~, nod3, wei3, ~, ~] = quadrature(N);
blist = basis_list(N, Np);
bb = box(x,y,z);
[Fk, ~, Jdet] = jacobians(x,y,z);
[phi, dphi] = basis(bb, blist, Fk, nod3);

err_L2 = 0;
err_H10 = 0;

for ie = 1:size(x,2)  
    for q = 1:length(wei3)
        
        % evaluation of the dg-fem solution and its gradient at the quadrature node
        u_h = phi(:,q,ie)'*u(:,ie);
        grad_uh = 0;
        for ii = 1:Np
            grad_uh = grad_uh + dphi(:,ii,q,ie)*u(ii,ie);
        end
        
        point_diff = uex_grad(Fk(:,:,ie)*nod3(:,q)) - grad_uh;
        err_H10 = err_H10 + wei3(q)*abs(Jdet(ie))*(point_diff'*point_diff);
        err_L2 = err_L2 + wei3(q)*abs(Jdet(ie))*(u_h - uex(Fk(:,:,ie)*nod3(:,q)))^2;
        
    end
end

err_L2 = sqrt(err_L2);
err_H10 = sqrt(err_H10);

end