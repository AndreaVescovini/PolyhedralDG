function [u, x, y, z] = linsys(mesh, f, gd, N, sigma, epsilon)
%[u, x, y, z] = linsys(mesh, f, gd, N, sigma, epsilon).
%Function that computes the matrix A and the vector b of the linear
%system that solves the problem.

% set vertices of tetrahedra, faces to be integrated on and bounding boxes
% of polyhedra.
[faces, faces_neig] = read_faces(mesh);
[x, y, z] = set_vertices(mesh.E2V, mesh.VX, mesh.VY, mesh.VZ);
[bb, hk] = bbox(x,y,z, mesh.E2P, mesh.K);

[Fk, Jinv, Jdet] = jacobians(x,y,z); % del determinante dovro' poi prenderne il valore assoluto

[nod2, wei2, nod3, wei3, node_maps, node_maps_inv] = quadrature(N);

% Computed values of basis functions.
Np = (N+1)*(N+2)*(N+3)/6; %number of dof for every element
blist = basis_list(N, Np);
[phi, dphi] = basis(bb, mesh.E2P, blist, Fk, nod3);
[phi_bordo, grad_bordo] = basis_boundary(bb, mesh.E2P, faces_neig, blist, Fk, node_maps, nod2);

A = spalloc(mesh.K*Np, mesh.K*Np, mesh.K*Np*Np);
b = zeros(mesh.K*Np,1);

for ie = 1:mesh.Ntet %loop ever tetrahedra
    
    % Recover the indexes relative to the element inside which ie is.
    index = (mesh.E2P(ie)-1)*Np+1:mesh.E2P(ie)*Np;
    
    for q = 1:length(wei3) %loop on 3D quadrature points
        A(index, index) = A(index, index) + (wei3(q)*abs(Jdet(ie)))*(dphi(:,:,q,ie)'*dphi(:,:,q,ie));
        b(index) = b(index) + abs(Jdet(ie))*wei3(q)*f(Fk(:,:,ie)*nod3(:,q))*phi(:,q,ie);
    end
end

for e = 1:size(faces,1) % loop over faces
    
    % Recover the neighboring tetrahedra the the local faces number of e
    E1 = faces_neig(e,1);
    E2 = faces_neig(e,3);    
    e_E1 = faces_neig(e,2);
    e_E2 = faces_neig(e,4);
    
    % Calulate the area of e and the normal in the direction E1 -> E2
    [area, normal] = metric2D(mesh.VX(faces(e,:)), mesh.VY(faces(e,:)), mesh.VZ(faces(e,:)), e_E1);
    if E2 == 0
        sig = sigma*N^2/hk(mesh.E2P(E1));
    else
        sig = sigma*N^2/min(mesh.E2P([E1, E2]));
    end
    
    % Recover the indexes relative to the element inside which E1 is.
    index1 = (mesh.E2P(E1)-1)*Np+1:mesh.E2P(E1)*Np;
       
    for q = 1:length(wei2) %loop on 2D quadrature nodes         
        if E2 == 0 % if true then e is a boundary face
                        
            A(index1, index1) = A(index1, index1) + sig*wei2(q)*phi_bordo(:,q,e,1)*phi_bordo(:,q,e,1)'*area...; % S
                        + epsilon*wei2(q)*(grad_bordo(:,:,q,e,1)'*normal)*phi_bordo(:,q,e,1)'*area... % I
                        - wei2(q)*phi_bordo(:,q,e,1)*(normal'*grad_bordo(:,:,q,e,1))*area; % I'
                    
            % evaluate gd on the physical quadrature nodes on the boundary
                  b(index1) = b(index1) + sig*wei2(q)*gd(Fk(:,:,E1)*node_maps(:,:,e_E1)*nod2(:,q))*phi_bordo(:,q,e,1)*area... % bs
                                  + epsilon*wei2(q)*gd(Fk(:,:,E1)*node_maps(:,:,e_E1)*nod2(:,q))*(grad_bordo(:,:,q,e,1)'*normal)*area; % bi
                        
        else
            % Recover the indexes relative to the element inside which E2 is.
            index2 = (mesh.E2P(E2)-1)*Np+1:mesh.E2P(E2)*Np;
            q_star = node_maps_inv(:,:,e_E2)*Jinv(:,:,E2)*(-Fk(:,4,E2)+Fk(:,:,E1)*node_maps(:,:,e_E1)*nod2(:,q));
            q_2 = find((q_star(1)+0.00001 > nod2(1,:)) & (nod2(1,:) > q_star(1)-0.00001) & (q_star(2)+0.00001 > nod2(2,:)) & (nod2(2,:) > q_star(2)-0.00001));
                       
            A(index1, index1) = A(index1, index1) + sig*wei2(q)*phi_bordo(:,q,e,1)*phi_bordo(:,q,e,1)'*area... % S1
                              + epsilon*0.5*wei2(q)*(grad_bordo(:,:,q,e,1)'*normal)*phi_bordo(:,q,e,1)'*area... % I1
                              - 0.5*wei2(q)*phi_bordo(:,q,e,1)*(normal'*grad_bordo(:,:,q,e,1))*area; % I'1
                    
            A(index2, index2) = A(index2, index2) + sig*wei2(q)*phi_bordo(:,q_2,e,2)*phi_bordo(:,q_2,e,2)'*area... % S2
                              + epsilon*0.5*wei2(q)*(grad_bordo(:,:,q_2,e,2)'*(-normal))*phi_bordo(:,q_2,e,2)'*area... % I2
                              - 0.5*wei2(q)*phi_bordo(:,q_2,e,2)*(-normal'*grad_bordo(:,:,q_2,e,2))*area; % I'2

            A(index1, index2) = A(index1, index2) - sig*wei2(q)*phi_bordo(:,q,e,1)*phi_bordo(:,q_2,e,2)'*area... % S1
                              - epsilon*0.5*wei2(q)*(grad_bordo(:,:,q,e,1)'*normal)*phi_bordo(:,q_2,e,2)'*area... % I1
                              + 0.5*wei2(q)*phi_bordo(:,q,e,1)*(-normal'*grad_bordo(:,:,q_2,e,2))*area; % I'2
                       
            A(index2, index1) = A(index2, index1) - sig*wei2(q)*phi_bordo(:,q_2,e,2)*phi_bordo(:,q,e,1)'*area... % S2
                              - epsilon*0.5*wei2(q)*(grad_bordo(:,:,q_2,e,2)'*(-normal))*phi_bordo(:,q,e,1)'*area... % I2
                              + 0.5*wei2(q)*phi_bordo(:,q_2,e,2)*(normal'*grad_bordo(:,:,q,e,1))*area; % I'1

        end     
    end
end

% sum up the contributions
% A = A + S + epsilon.*I - transpose(I);
% b = b + bs + epsilon.*bi;

% solve the linear system
u = A\b;
u = reshape(u, [Np, mesh.K]);
end