function [u, x, y, z] = linsys(mesh, f, gd, N, sigma, epsilon)
%[u, x, y, z] = linsys(mesh, f, gd, N, sigma, epsilon).
%Function that computes the matrix A and the vector b of the linear
%system that solves the problem.

Np = (N+1)*(N+2)*(N+3)/6; %number of nodes for every element
Nfaces = 4; % number of faces for every element

% dof in the reference simplex
[r, s, t] = set_dof(N);

% dof in the phisical elements
va = mesh.EToV(:,1)'; vb = mesh.EToV(:,2)'; vc = mesh.EToV(:,3)'; vd = mesh.EToV(:,4)';
x = (1-r-s-t)*mesh.VX(va)+r*mesh.VX(vb)+s*mesh.VX(vc)+t*mesh.VX(vd);
y = (1-r-s-t)*mesh.VY(va)+r*mesh.VY(vb)+s*mesh.VY(vc)+t*mesh.VY(vd);
z = (1-r-s-t)*mesh.VZ(va)+r*mesh.VZ(vb)+s*mesh.VZ(vc)+t*mesh.VZ(vd);

A = zeros(mesh.K*Np);
% I = zeros(mesh.K*Np);
% S = zeros(mesh.K*Np);
b = zeros(mesh.K*Np,1);
% bi = zeros(mesh.K*Np,1);
% bs = zeros(mesh.K*Np,1);

% A = spalloc(mesh.K*Np, mesh.K*Np, mesh.K*Np*Np);
% I = spalloc(mesh.K*Np, mesh.K*Np, mesh.K*Np*Np*(Nfaces+1));
% S = spalloc(mesh.K*Np, mesh.K*Np, mesh.K*Np*Np*(Nfaces+1));

[nod2, wei2, nod3, wei3, node_maps, node_maps_inv] = quadrature(Nfaces);
[phi, dphi, phi_bordo, grad_bordo, ~] = basis_evaluation(nod3, nod2, node_maps, N, Np);
% compute the jacobians
[trasl, J, Jcof, Jdet] = jacobians(x(1:4,:),y(1:4,:),z(1:4,:),r(1:4),s(1:4),t(1:4)); % del determinante dovro' poi prenderne il valore assoluto

for ie = 1:mesh.K %loop on the elements
    
    [areas, normals] = metric2D(x(1:4,ie), y(1:4,ie), z(1:4,ie), Nfaces);
    
    for q = 1:length(wei3) %loop on 3D quadrature points
        for i = 1:Np %loop on the basis functions
            igl = (ie-1)*Np + i; %global node number
            for j = 1:Np
                jgl = (ie-1)*Np + j;
                A(igl,jgl) = A(igl,jgl) + (wei3(q)/abs(Jdet(ie)))*(dphi(:,i,q)'*(Jcof(:,:,ie)'*Jcof(:,:,ie))*dphi(:,j,q));                                
            end
            b(igl) = b(igl) + abs(Jdet(ie))*wei3(q)*f(trasl(:,ie)+J(:,:,ie)*nod3(:,q))*phi(i,q);
        end
    end
    
    for e = 1:Nfaces
        sig = sigma/(areas(e)/2)^0.5;
        
        for q = 1:length(wei2) %loop on 2D quadrature nodes
            nei_q = node_maps_inv(:,:,mesh.EToF(ie,e))*Jcof(:,:,mesh.EToE(ie,e))'/Jdet(mesh.EToE(ie,e))*(-trasl(:,mesh.EToE(ie,e))+trasl(:,ie)+J(:,:,ie)*node_maps(:,:,e)*nod2(:,q));
%             nei_q = Jcof(:,:,mesh.EToE(ie,e))'/Jdet(mesh.EToE(ie,e))*(-trasl(:,mesh.EToE(ie,e))+trasl(:,ie)+J(:,:,ie)*node_maps(:,:,e)*nod2(:,q));
            q_star = find((nei_q(1)+0.0001 > nod2(1,:)) & (nod2(1,:) > nei_q(1)-0.0001) & (nei_q(2)+0.0001 > nod2(2,:)) & (nod2(2,:) > nei_q(2)-0.0001));
            for i = 1:Np
                igl = (ie-1)*Np + i;
                for j = 1:Np
                    jgl = (ie-1)*Np + j;
                      
                    if mesh.EToE(ie, e) == ie % if true then e is a boundary face
                        
                          A(igl, jgl) = A(igl, jgl) + sig*wei2(q)*phi_bordo(i,q,e)*phi_bordo(j,q,e)*areas(e)...; % S
                                                    + epsilon*wei2(q)*(normals(:,e)'*Jcof(:,:,ie)*grad_bordo(:,i,q,e)/Jdet(ie))*phi_bordo(j,q,e)*areas(e); % I
                          A(jgl, igl) = A(jgl, igl) - wei2(q)*(normals(:,e)'*Jcof(:,:,ie)*grad_bordo(:,i,q,e)/Jdet(ie))*phi_bordo(j,q,e)*areas(e); % I'
                        
                    else
                        jneigh = (mesh.EToE(ie, e)-1)*Np + j; %column index corresponding to the neighbooring element
                       
                          A(igl, jgl) = A(igl, jgl) + sig*wei2(q)*phi_bordo(i,q,e)*phi_bordo(j,q,e)*areas(e)... % S
                                                    + epsilon*0.5*wei2(q)*(normals(:,e)'*Jcof(:,:,ie)*grad_bordo(:,i,q,e)/Jdet(ie))*phi_bordo(j,q,e)*areas(e); % I
                          A(jgl, igl) = A(jgl, igl) - 0.5*wei2(q)*(normals(:,e)'*Jcof(:,:,ie)*grad_bordo(:,i,q,e)/Jdet(ie))*phi_bordo(j,q,e)*areas(e); % I'

%                           A(igl, jneigh) = A(igl, jneigh) - sig*wei2(q)*phi_bordo(q, e, i)*phi_func{j}(nei_q)*areas(e)... % S
%                                                           - epsilon*0.5*wei2(q)*(normals(:,e)'*Jcof(:,:,ie)*grad_bordo(:,q,e,i)/Jdet(ie))*phi_func{j}(nei_q)*areas(e); % I                         
%                           A(jneigh, igl) = A(jneigh, igl) + 0.5*wei2(q)*(normals(:,e)'*Jcof(:,:,ie)*grad_bordo(:,q,e,i)/Jdet(ie))*phi_func{j}(nei_q)*areas(e); % I'
                          A(igl, jneigh) = A(igl, jneigh) - sig*wei2(q)*phi_bordo(i,q,e)*phi_bordo(j,q_star,mesh.EToF(ie,e))*areas(e)... % S
                                                          - epsilon*0.5*wei2(q)*(normals(:,e)'*Jcof(:,:,ie)*grad_bordo(:,i,q,e)/Jdet(ie))*phi_bordo(j,q_star,mesh.EToF(ie,e))*areas(e); % I                         
                          A(jneigh, igl) = A(jneigh, igl) + 0.5*wei2(q)*(normals(:,e)'*Jcof(:,:,ie)*grad_bordo(:,i,q,e)/Jdet(ie))*phi_bordo(j,q_star,mesh.EToF(ie,e))*areas(e); % I'

                    end
                end
                if mesh.EToE(ie, e) == ie % boundary face
                    % evaluate gd on the physical quadrature nodes on the boundary
                      b(igl) = b(igl) + sig*wei2(q)*gd(trasl(:,ie)+J(:,:,ie)*node_maps(:,:,e)*nod2(:,q))*phi_bordo(i,q,e)*areas(e)... % bs
                                      + epsilon*wei2(q)*gd(trasl(:,ie)+J(:,:,ie)*node_maps(:,:,e)*nod2(:,q))*(normals(:,e)'*Jcof(:,:,ie)*grad_bordo(:,i,q,e)/Jdet(ie))*areas(e); % bi
                   
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