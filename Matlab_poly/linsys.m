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

[nod2, wei2, nod3, wei3, node_maps, node_maps_inv] = quadrature(N,Nfaces);
% compute the jacobians
[Fk, Jinv, Jdet] = jacobians(x,y,z); % del determinante dovro' poi prenderne il valore assoluto
bb = box(x,y,z);
[phi, dphi] = basis(bb, blist, Fk, nod3);
[phi_bordo, grad_bordo] = basis_boundary(bb, blist, Fk, node_maps, nod2);

for ie = 1:mesh.K %loop on the elements
    
    [areas, normals] = metric2D(x(:,ie), y(:,ie), z(:,ie), Nfaces);
    index = (ie-1)*Np+1:ie*Np;
    
    for q = 1:length(wei3) %loop on 3D quadrature points
        A(index, index) = A(index, index) + (wei3(q)*abs(Jdet(ie)))*(dphi(:,:,q,ie)'*dphi(:,:,q,ie));
        b(index) = b(index) + abs(Jdet(ie))*wei3(q)*f(Fk(:,:,ie)*nod3(:,q))*phi(:,q,ie);
    end
    
    for e = 1:Nfaces
        sig = sigma/(areas(e)/2)^0.5;
        
        for q = 1:length(wei2) %loop on 2D quadrature nodes
            nei_q = node_maps_inv(:,:,mesh.EToF(ie,e))*Jinv(:,:,mesh.EToE(ie,e))*(-Fk(:,4,mesh.EToE(ie,e))+Fk(:,:,ie)*node_maps(:,:,e)*nod2(:,q));
%             nei_q = Jcof(:,:,mesh.EToE(ie,e))'/Jdet(mesh.EToE(ie,e))*(-trasl(:,mesh.EToE(ie,e))+trasl(:,ie)+J(:,:,ie)*node_maps(:,:,e)*nod2(:,q));
            q_star = find((nei_q(1)+0.0001 > nod2(1,:)) & (nod2(1,:) > nei_q(1)-0.0001) & (nei_q(2)+0.0001 > nod2(2,:)) & (nod2(2,:) > nei_q(2)-0.0001));
            
            if mesh.EToE(ie, e) == ie % if true then e is a boundary face
                        
                A(index, index) = A(index, index) + sig*wei2(q)*phi_bordo(:,q,e,ie)*phi_bordo(:,q,e,ie)'*areas(e)...; % S
                            + epsilon*wei2(q)*(grad_bordo(:,:,q,e,ie)'*normals(:,e))*phi_bordo(:,q,e,ie)'*areas(e)... % I
                            - wei2(q)*phi_bordo(:,q,e,ie)*(normals(:,e)'*grad_bordo(:,:,q,e,ie))*areas(e); % I'
                        
            else
                index_neig = (mesh.EToE(ie, e)-1)*Np+1:mesh.EToE(ie, e)*Np; %column index corresponding to the neighbooring element
                       
                A(index, index) = A(index, index) + sig*wei2(q)*phi_bordo(:,q,e,ie)*phi_bordo(:,q,e,ie)'*areas(e)... % S
                            + epsilon*0.5*wei2(q)*(grad_bordo(:,:,q,e,ie)'*normals(:,e))*phi_bordo(:,q,e,ie)'*areas(e)... % I
                            - 0.5*wei2(q)*phi_bordo(:,q,e,ie)*(normals(:,e)'*grad_bordo(:,:,q,e,ie))*areas(e); % I'

                A(index, index_neig) = A(index, index_neig) - sig*wei2(q)*phi_bordo(:,q,e,ie)*phi_bordo(:,q_star,mesh.EToF(ie,e),mesh.EToE(ie,e))'*areas(e)... % S
                               - epsilon*0.5*wei2(q)*(grad_bordo(:,:,q,e,ie)'*normals(:,e))*phi_bordo(:,q_star,mesh.EToF(ie,e),mesh.EToE(ie,e))'*areas(e); % I
                A(index_neig, index) = A(index_neig, index) + 0.5*wei2(q)*phi_bordo(:,q_star,mesh.EToF(ie,e),mesh.EToE(ie,e))*(normals(:,e)'*grad_bordo(:,:,q,e,ie))*areas(e); % I'

            end
 
            if mesh.EToE(ie, e) == ie % boundary face
                    % evaluate gd on the physical quadrature nodes on the boundary
                      b(index) = b(index) + sig*wei2(q)*gd(Fk(:,:,ie)*node_maps(:,:,e)*nod2(:,q))*phi_bordo(:,q,e,ie)*areas(e)... % bs
                                      + epsilon*wei2(q)*gd(Fk(:,:,ie)*node_maps(:,:,e)*nod2(:,q))*(grad_bordo(:,:,q,e,ie)'*normals(:,e))*areas(e); % bi
                   
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