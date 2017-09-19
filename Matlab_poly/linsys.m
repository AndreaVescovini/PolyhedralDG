function [u, x, y, z] = linsys(mesh, f, gd, N, sigma, epsilon)
%[u, x, y, z] = linsys(mesh, f, gd, N, sigma, epsilon).
%Function that computes the matrix A and the vector b of the linear
%system that solves the problem.

Np = (N+1)*(N+2)*(N+3)/6; %number of dof for every element
blist = basis_list(N, Np);

% set vertices and faces
[x, y, z] = set_vertices(mesh.EToV, mesh.VX, mesh.VY, mesh.VZ);
[faces, faces_neig, elem_faces] = read_faces(mesh);

% A = zeros(mesh.K*Np);
A = spalloc(mesh.K*Np, mesh.K*Np, mesh.K*Np*Np);
b = zeros(mesh.K*Np,1);

[nod2, wei2, nod3, wei3, node_maps, node_maps_inv] = quadrature(N);
% compute the jacobians
[Fk, Jinv, Jdet] = jacobians(x,y,z); % del determinante dovro' poi prenderne il valore assoluto
bb = box(x,y,z);
[phi, dphi] = basis(bb, blist, Fk, nod3);
[phi_bordo, grad_bordo] = basis_boundary(bb, blist, Fk, node_maps, nod2);

for ie = 1:mesh.K %loop on the elements
    
    index = (ie-1)*Np+1:ie*Np;
    
    for q = 1:length(wei3) %loop on 3D quadrature points
        A(index, index) = A(index, index) + (wei3(q)*abs(Jdet(ie)))*(dphi(:,:,q,ie)'*dphi(:,:,q,ie));
        b(index) = b(index) + abs(Jdet(ie))*wei3(q)*f(Fk(:,:,ie)*nod3(:,q))*phi(:,q,ie);
    end
end

for e = 1:size(faces,1)
    
    E1 = faces_neig(e,1);
    e_E1 = faces_neig(e,2);
    E2 = faces_neig(e,3);
    e_E2 = faces_neig(e,4);
    [area, normal] = metric2D(mesh.VX(faces(e,:)), mesh.VY(faces(e,:)), mesh.VZ(faces(e,:)), e_E1);
    sig = sigma/(area/2)^0.5;
    
    index1 = (E1-1)*Np+1:E1*Np;
       
    for q = 1:length(wei2) %loop on 2D quadrature nodes         
        if E2 == 0 % if true then e is a boundary face
                        
            A(index1, index1) = A(index1, index1) + sig*wei2(q)*phi_bordo(:,q,e_E1,E1)*phi_bordo(:,q,e_E1,E1)'*area...; % S
                        + epsilon*wei2(q)*(grad_bordo(:,:,q,e_E1,E1)'*normal)*phi_bordo(:,q,e_E1,E1)'*area... % I
                        - wei2(q)*phi_bordo(:,q,e_E1,E1)*(normal'*grad_bordo(:,:,q,e_E1,E1))*area; % I'
                    
            % evaluate gd on the physical quadrature nodes on the boundary
                  b(index1) = b(index1) + sig*wei2(q)*gd(Fk(:,:,E1)*node_maps(:,:,e_E1)*nod2(:,q))*phi_bordo(:,q,e_E1,E1)*area... % bs
                                  + epsilon*wei2(q)*gd(Fk(:,:,E1)*node_maps(:,:,e_E1)*nod2(:,q))*(grad_bordo(:,:,q,e_E1,E1)'*normal)*area; % bi
                        
        else
            index2 = (E2-1)*Np+1:E2*Np;
            q_star = node_maps_inv(:,:,e_E2)*Jinv(:,:,E2)*(-Fk(:,4,E2)+Fk(:,:,E1)*node_maps(:,:,e_E1)*nod2(:,q));
            q_2 = find((q_star(1)+0.0001 > nod2(1,:)) & (nod2(1,:) > q_star(1)-0.0001) & (q_star(2)+0.0001 > nod2(2,:)) & (nod2(2,:) > q_star(2)-0.0001));
                       
            A(index1, index1) = A(index1, index1) + sig*wei2(q)*phi_bordo(:,q,e_E1,E1)*phi_bordo(:,q,e_E1,E1)'*area... % S1
                              + epsilon*0.5*wei2(q)*(grad_bordo(:,:,q,e_E1,E1)'*normal)*phi_bordo(:,q,e_E1,E1)'*area... % I1
                              - 0.5*wei2(q)*phi_bordo(:,q,e_E1,E1)*(normal'*grad_bordo(:,:,q,e_E1,E1))*area; % I'1
                    
            A(index2, index2) = A(index2, index2) + sig*wei2(q)*phi_bordo(:,q_2,e_E2,E2)*phi_bordo(:,q_2,e_E2,E2)'*area... % S2
                              + epsilon*0.5*wei2(q)*(grad_bordo(:,:,q_2,e_E2,E2)'*(-normal))*phi_bordo(:,q_2,e_E2,E2)'*area... % I2
                              - 0.5*wei2(q)*phi_bordo(:,q_2,e_E2,E2)*(-normal'*grad_bordo(:,:,q_2,e_E2,E2))*area; % I'2

            A(index1, index2) = A(index1, index2) - sig*wei2(q)*phi_bordo(:,q,e_E1,E1)*phi_bordo(:,q_2,e_E2,E2)'*area... % S1
                              - epsilon*0.5*wei2(q)*(grad_bordo(:,:,q,e_E1,E1)'*normal)*phi_bordo(:,q_2,e_E2,E2)'*area... % I1
                              + 0.5*wei2(q)*phi_bordo(:,q,e_E1,E1)*(-normal'*grad_bordo(:,:,q_2,e_E2,E2))*area; % I'2
                       
            A(index2, index1) = A(index2, index1) - sig*wei2(q)*phi_bordo(:,q_2,e_E2,E2)*phi_bordo(:,q,e_E1,E1)'*area... % S2
                              - epsilon*0.5*wei2(q)*(grad_bordo(:,:,q_2,e_E2,E2)'*(-normal))*phi_bordo(:,q,e_E1,E1)'*area... % I2
                              + 0.5*wei2(q)*phi_bordo(:,q_2,e_E2,E2)*(normal'*grad_bordo(:,:,q,e_E1,E1))*area; % I'1

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