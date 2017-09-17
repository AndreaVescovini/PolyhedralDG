function [u, x, y, z] = linsys(mesh, f, gd, N, sigma, epsilon)
%[u, x, y, z] = linsys(mesh, f, gd, N, sigma, epsilon).
%Function that computes the matrix A and the vector b of the linear
%system that solves the problem.

Np = (N+1)*(N+2)*(N+3)/6; %number of dof for every element
Nfaces = 4; % number of faces for every element

blist = basis_list(N, Np);

% set vertices
[x, y, z] = set_vertices(mesh.EToV, mesh.VX, mesh.VY, mesh.VZ);

A = zeros(mesh.K*Np);
% A = spalloc(mesh.K*Np, mesh.K*Np, mesh.K*Np*Np);
b = zeros(mesh.K*Np,1);

[nod2, wei2, nod3, wei3, node_maps, node_maps_inv] = quadrature(Nfaces);
% compute the jacobians
[Fk, Jcof, Jdet] = jacobians(x,y,z); % del determinante dovro' poi prenderne il valore assoluto
bb = box(x,y,z);
[phi, dphi, phi_bordo, grad_bordo] = basis(bb, blist, Fk, nod3, node_maps, nod2);

for ie = 1:mesh.K %loop on the elements
    
    [areas, normals] = metric2D(x(:,ie), y(:,ie), z(:,ie), Nfaces);
    
    for q = 1:length(wei3) %loop on 3D quadrature points
        for i = 1:Np %loop on the basis functions
            igl = (ie-1)*Np + i; %global node number
            for j = 1:Np
                jgl = (ie-1)*Np + j;
                A(igl,jgl) = A(igl,jgl) + (wei3(q)*abs(Jdet(ie)))*(dphi(:,q,j,ie)'*dphi(:,q,i,ie));                                
            end
            b(igl) = b(igl) + abs(Jdet(ie))*wei3(q)*f(Fk(:,:,ie)*nod3(:,q))*phi(q,i,ie);
        end
    end
    
    for e = 1:Nfaces
        sig = sigma/(areas(e)/2)^0.5;
        
        for q = 1:length(wei2) %loop on 2D quadrature nodes
            nei_q = node_maps_inv(:,:,mesh.EToF(ie,e))*Jcof(:,:,mesh.EToE(ie,e))'/Jdet(mesh.EToE(ie,e))*(-Fk(:,4,mesh.EToE(ie,e))+Fk(:,:,ie)*node_maps(:,:,e)*nod2(:,q));
%             nei_q = Jcof(:,:,mesh.EToE(ie,e))'/Jdet(mesh.EToE(ie,e))*(-trasl(:,mesh.EToE(ie,e))+trasl(:,ie)+J(:,:,ie)*node_maps(:,:,e)*nod2(:,q));
            q_star = find((nei_q(1)+0.0001 > nod2(1,:)) & (nod2(1,:) > nei_q(1)-0.0001) & (nei_q(2)+0.0001 > nod2(2,:)) & (nod2(2,:) > nei_q(2)-0.0001));
            for i = 1:Np
                igl = (ie-1)*Np + i;
                for j = 1:Np
                    jgl = (ie-1)*Np + j;
                      
                    if mesh.EToE(ie, e) == ie % if true then e is a boundary face
                        
                          A(igl, jgl) = A(igl, jgl) + sig*wei2(q)*phi_bordo(q,e,i,ie)*phi_bordo(q,e,j,ie)*areas(e)...; % S
                                                    + epsilon*wei2(q)*(normals(:,e)'*grad_bordo(:,q,e,i,ie))*phi_bordo(q,e,j,ie)*areas(e); % I
                          A(jgl, igl) = A(jgl, igl) - wei2(q)*(normals(:,e)'*grad_bordo(:,q,e,i,ie))*phi_bordo(q,e,j,ie)*areas(e); % I'
                        
                    else
                        jneigh = (mesh.EToE(ie, e)-1)*Np + j; %column index corresponding to the neighbooring element
                       
                          A(igl, jgl) = A(igl, jgl) + sig*wei2(q)*phi_bordo(q, e, i,ie)*phi_bordo(q, e, j,ie)*areas(e)... % S
                                                    + epsilon*0.5*wei2(q)*(normals(:,e)'*grad_bordo(:,q,e,i,ie))*phi_bordo(q,e,j,ie)*areas(e); % I
                          A(jgl, igl) = A(jgl, igl) - 0.5*wei2(q)*(normals(:,e)'*grad_bordo(:,q,e,i,ie))*phi_bordo(q,e,j,ie)*areas(e); % I'

                          A(igl, jneigh) = A(igl, jneigh) - sig*wei2(q)*phi_bordo(q,e,i,ie)*phi_bordo(q_star,mesh.EToF(ie,e),j,mesh.EToE(ie,e))*areas(e)... % S
                                                          - epsilon*0.5*wei2(q)*(normals(:,e)'*grad_bordo(:,q,e,i,ie))*phi_bordo(q_star, mesh.EToF(ie,e),j,mesh.EToE(ie,e))*areas(e); % I                         
                          A(jneigh, igl) = A(jneigh, igl) + 0.5*wei2(q)*(normals(:,e)'*grad_bordo(:,q,e,i,ie))*phi_bordo(q_star,mesh.EToF(ie,e),j,mesh.EToE(ie,e))*areas(e); % I'

                    end
                end
                if mesh.EToE(ie, e) == ie % boundary face
                    % evaluate gd on the physical quadrature nodes on the boundary
                      b(igl) = b(igl) + sig*wei2(q)*gd(Fk(:,:,ie)*node_maps(:,:,e)*nod2(:,q))*phi_bordo(q,e,i,ie)*areas(e)... % bs
                                      + epsilon*wei2(q)*gd(Fk(:,:,ie)*node_maps(:,:,e)*nod2(:,q))*(normals(:,e)'*grad_bordo(:,q,e,i,ie))*areas(e); % bi
                   
                end              
            end
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