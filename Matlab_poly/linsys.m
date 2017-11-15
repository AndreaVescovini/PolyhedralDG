function [u] = linsys(geom, f, gd, N, sigma, epsilon)
%[u] = linsys(geom, f, gd, N, sigma, epsilon).
%Function that computes the matrix A and the vector b of the linear
%system that solves the problem.
%
% Author: Andrea Vescovini

fprintf('Computing the basis...');

% get quadrature nodes for the required order
[nod2, wei2, nod3, wei3, node_maps, node_maps_inv] = quadrature(N);

% Computed values of basis functions.
Np = (N+1)*(N+2)*(N+3)/6; %number of dof for every element
blist = basis_list(N, Np);
[phi, dphi] = basis(geom.bb, geom.E2P, blist, geom.Fk, nod3);
[phi_bordo, grad_bordo] = basis_boundary(geom.bb, geom.E2P, geom.faces_neig, blist, geom.Fk, node_maps, nod2);

fprintf('done!\nAssembling the linear system...');

% Allocate matrix and vector
A = spalloc(geom.K*Np, geom.K*Np, geom.K*Np*Np);
b = zeros(geom.K*Np,1);

for ie = 1:geom.Ntet %loop ever tetrahedra

    % Recover the indexes relative to the element inside which ie is.
    index = (geom.E2P(ie)-1)*Np+1:geom.E2P(ie)*Np;

    for q = 1:length(wei3) %loop on 3D quadrature points
        A(index, index) = A(index, index) + (wei3(q)*abs(geom.Jdet(ie)))*(dphi(:,:,q,ie)'*dphi(:,:,q,ie));
        b(index) = b(index) + abs(geom.Jdet(ie))*wei3(q)*f(geom.Fk(:,:,ie)*nod3(:,q))*phi(:,q,ie);
    end
end

for e = 1:geom.Nfaces % loop over faces

    % Recover the neighboring tetrahedra the the local faces number of e
    E1 = geom.faces_neig(e,1);
    E2 = geom.faces_neig(e,3);
    e_E1 = geom.faces_neig(e,2);
    e_E2 = geom.faces_neig(e,4);

    if E2 == 0
        sig = sigma*N^2/geom.hk(geom.E2P(E1));
    else
        sig = sigma*N^2/min(geom.hk(geom.E2P([E1, E2])));
    end

    % Recover the indexes relative to the element inside which E1 is.
    index1 = (geom.E2P(E1)-1)*Np+1:geom.E2P(E1)*Np;

    for q = 1:length(wei2) %loop on 2D quadrature nodes
        if E2 == 0 % if true then e is a boundary face

            A(index1, index1) = A(index1, index1) + sig*wei2(q)*phi_bordo(:,q,e,1)*phi_bordo(:,q,e,1)'*geom.area(e)...; % S
                        + epsilon*wei2(q)*(grad_bordo(:,:,q,e,1)'*geom.normal(:,e))*phi_bordo(:,q,e,1)'*geom.area(e)... % I
                        - wei2(q)*phi_bordo(:,q,e,1)*(geom.normal(:,e)'*grad_bordo(:,:,q,e,1))*geom.area(e); % I'

            % evaluate gd on the physical quadrature nodes on the boundary
                  b(index1) = b(index1) + sig*wei2(q)*gd(geom.Fk(:,:,E1)*node_maps(:,:,e_E1)*nod2(:,q))*phi_bordo(:,q,e,1)*geom.area(e)... % bs
                                  + epsilon*wei2(q)*gd(geom.Fk(:,:,E1)*node_maps(:,:,e_E1)*nod2(:,q))*(grad_bordo(:,:,q,e,1)'*geom.normal(:,e))*geom.area(e); % bi

        else
            % Recover the indexes relative to the element inside which E2 is.
            index2 = (geom.E2P(E2)-1)*Np+1:geom.E2P(E2)*Np;
            q_star = node_maps_inv(:,:,e_E2)*geom.Jinv(:,:,E2)*(-geom.Fk(:,4,E2)+geom.Fk(:,:,E1)*node_maps(:,:,e_E1)*nod2(:,q));
            q_2 = find((q_star(1)+0.00001 > nod2(1,:)) & (nod2(1,:) > q_star(1)-0.00001) & (q_star(2)+0.00001 > nod2(2,:)) & (nod2(2,:) > q_star(2)-0.00001));

            A(index1, index1) = A(index1, index1) + sig*wei2(q)*phi_bordo(:,q,e,1)*phi_bordo(:,q,e,1)'*geom.area(e)... % S1
                              + epsilon*0.5*wei2(q)*(grad_bordo(:,:,q,e,1)'*geom.normal(:,e))*phi_bordo(:,q,e,1)'*geom.area(e)... % I1
                              - 0.5*wei2(q)*phi_bordo(:,q,e,1)*(geom.normal(:,e)'*grad_bordo(:,:,q,e,1))*geom.area(e); % I'1

            A(index2, index2) = A(index2, index2) + sig*wei2(q)*phi_bordo(:,q_2,e,2)*phi_bordo(:,q_2,e,2)'*geom.area(e)... % S2
                              + epsilon*0.5*wei2(q)*(grad_bordo(:,:,q_2,e,2)'*(-geom.normal(:,e)))*phi_bordo(:,q_2,e,2)'*geom.area(e)... % I2
                              - 0.5*wei2(q)*phi_bordo(:,q_2,e,2)*(-geom.normal(:,e)'*grad_bordo(:,:,q_2,e,2))*geom.area(e); % I'2

            A(index1, index2) = A(index1, index2) - sig*wei2(q)*phi_bordo(:,q,e,1)*phi_bordo(:,q_2,e,2)'*geom.area(e)... % S1
                              - epsilon*0.5*wei2(q)*(grad_bordo(:,:,q,e,1)'*geom.normal(:,e))*phi_bordo(:,q_2,e,2)'*geom.area(e)... % I1
                              + 0.5*wei2(q)*phi_bordo(:,q,e,1)*(-geom.normal(:,e)'*grad_bordo(:,:,q_2,e,2))*geom.area(e); % I'2

            A(index2, index1) = A(index2, index1) - sig*wei2(q)*phi_bordo(:,q_2,e,2)*phi_bordo(:,q,e,1)'*geom.area(e)... % S2
                              - epsilon*0.5*wei2(q)*(grad_bordo(:,:,q_2,e,2)'*(-geom.normal(:,e)))*phi_bordo(:,q,e,1)'*geom.area(e)... % I2
                              + 0.5*wei2(q)*phi_bordo(:,q_2,e,2)*(geom.normal(:,e)'*grad_bordo(:,:,q,e,1))*geom.area(e); % I'1

        end
    end
end

% spy(A, 'k');

fprintf('done!\nSolving the linear system...');

% solve the linear system
u = A\b;
u = reshape(u, [Np, geom.K]);

fprintf('done!\n');
end
