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
% b = zeros(mesh.K*Np,1);
% bi = zeros(mesh.K*Np,1);
% bs = zeros(mesh.K*Np,1);

% A = spalloc(mesh.K*Np, mesh.K*Np, mesh.K*Np*Np);
% I = spalloc(mesh.K*Np, mesh.K*Np, mesh.K*Np*Np*(Nfaces+1));
% S = spalloc(mesh.K*Np, mesh.K*Np, mesh.K*Np*Np*(Nfaces+1));
b = zeros(mesh.K*Np,1);
% bi = zeros(mesh.K*Np,1);
% bs = zeros(mesh.K*Np,1);


[nod2, wei2, nod3, wei3, node_maps] = quadrature(Nfaces);
[phi, dphi, phi_bordo, grad_bordo] = basis_evaluation(nod3, nod2, node_maps, N, Np);

for ie = 1:mesh.K %loop on the elements
    
    % compute the jacobians
    [J, Jcof, Jdet, trasl] = jacobians(x(1:4,ie),y(1:4,ie),z(1:4,ie),r(1:4),s(1:4),t(1:4)); % del determinante dovro' poi prenderne il valore assoluto
    [areas, normals] = metric2D(x(1:4,ie), y(1:4,ie), z(1:4,ie), Nfaces);
    
    for q = 1:length(wei3) %loop on 3D quadrature points
        for i = 1:Np %loop on the basis functions
            igl = (ie-1)*Np + i; %global node number
            for j = 1:Np
                jgl = (ie-1)*Np + j;
                A(igl,jgl) = A(igl,jgl) + (wei3(q)/abs(Jdet))*(dphi(:,q,j)'*(Jcof'*Jcof)*dphi(:,q,i));                                
            end
            b(igl) = b(igl) + abs(Jdet)*wei3(q)*f(trasl+J*nod3(:,q))*phi(q,i);
        end
    end
    
    for e = 1:Nfaces
        sig = sigma/(areas(e)/2)^0.5;
        for q = 1:length(wei2) %loop on 2D quadrature nodes
            for i = 1:Np
                igl = (ie-1)*Np + i;
                for j = 1:Np
                    jgl = (ie-1)*Np + j;
                      
                    if mesh.EToE(ie, e) == ie % if true then e is a boundary face
%                         S(igl, jgl) = S(igl, jgl) + sig*wei2(q)*phi_bordo(q, e, i)*phi_bordo(q, e, j)*areas(e);
%                         I(igl, jgl) = I(igl, jgl) + wei2(q)*(normals(:,e)'*Jcof*grad_bordo(:,q,e,i)/Jdet)*phi_bordo(q,e,j)*areas(e);
                          A(igl, jgl) = A(igl, jgl) + sig*wei2(q)*phi_bordo(q, e, i)*phi_bordo(q, e, j)*areas(e)...; % S
                                                    + epsilon*wei2(q)*(normals(:,e)'*Jcof*grad_bordo(:,q,e,i)/Jdet)*phi_bordo(q,e,j)*areas(e); % I
                          A(jgl, igl) = A(jgl, igl) - wei2(q)*(normals(:,e)'*Jcof*grad_bordo(:,q,e,i)/Jdet)*phi_bordo(q,e,j)*areas(e); % I'
                        
                    else
                        jneigh = (mesh.EToE(ie, e)-1)*Np + j; %column index corresponding to the neighbooring element
                        %j_star = j_map(mesh.EToV(ie,:) == mesh.EToV(mesh.EToE(ie,e),j), e);
                        j_star = j_map(x,y,z,ie,j,mesh.EToE(ie,e),e);
                        %S(igl, jneigh) = S(igl, jneigh) - sig*wei2(q)*phi_bordo(q, e, i)*phi_bordo(q,EToF(ie, e), j)*areas(e, ie); %old and wrong
                        
%                         S(igl, jgl) = S(igl, jgl) + sig*wei2(q)*phi_bordo(q, e, i)*phi_bordo(q, e, j)*areas(e);
%                         I(igl, jgl) = I(igl, jgl) + 0.5*wei2(q)*(normals(:,e)'*Jcof*grad_bordo(:,q,e,i)/Jdet)*phi_bordo(q,e,j)*areas(e);
                          A(igl, jgl) = A(igl, jgl) + sig*wei2(q)*phi_bordo(q, e, i)*phi_bordo(q, e, j)*areas(e)... % S
                                                    + epsilon*0.5*wei2(q)*(normals(:,e)'*Jcof*grad_bordo(:,q,e,i)/Jdet)*phi_bordo(q,e,j)*areas(e); % I
                          A(jgl, igl) = A(jgl, igl) - 0.5*wei2(q)*(normals(:,e)'*Jcof*grad_bordo(:,q,e,i)/Jdet)*phi_bordo(q,e,j)*areas(e); % I'
%                         S(igl, jneigh) = S(igl, jneigh) - sig*wei2(q)*phi_bordo(q, e, i)*phi_bordo(q, e, j_star)*areas(e);                    
%                         I(igl, jneigh) = I(igl, jneigh) - 0.5*wei2(q)*(normals(:,e)'*Jcof*grad_bordo(:,q,e,i)/Jdet)*phi_bordo(q, e, j_star)*areas(e);
                          A(igl, jneigh) = A(igl, jneigh) - sig*wei2(q)*phi_bordo(q, e, i)*phi_bordo(q, e, j_star)*areas(e)... % S
                                                          - epsilon*0.5*wei2(q)*(normals(:,e)'*Jcof*grad_bordo(:,q,e,i)/Jdet)*phi_bordo(q, e, j_star)*areas(e); % I                         
                          A(jneigh, igl) = A(jneigh, igl) + 0.5*wei2(q)*(normals(:,e)'*Jcof*grad_bordo(:,q,e,i)/Jdet)*phi_bordo(q, e, j_star)*areas(e); % I'

                    end
                end
                if mesh.EToE(ie, e) == ie % boundary face
                    % evaluate gd on the physical quadrature nodes on the boundary
%                     bs(igl) = bs(igl) + sig*wei2(q)*gd(trasl+J*node_maps(:,:,e)*nod2(:,q))*phi_bordo(q,e,i)*areas(e);
%                     bi(igl) = bi(igl) + wei2(q)*gd(trasl+J*node_maps(:,:,e)*nod2(:,q))*(normals(:,e)'*Jcof*grad_bordo(:,q,e,i)/Jdet)*areas(e);
                      b(igl) = b(igl) + sig*wei2(q)*gd(trasl+J*node_maps(:,:,e)*nod2(:,q))*phi_bordo(q,e,i)*areas(e)... % bs
                                      + epsilon*wei2(q)*gd(trasl+J*node_maps(:,:,e)*nod2(:,q))*(normals(:,e)'*Jcof*grad_bordo(:,q,e,i)/Jdet)*areas(e); % bi
                   
                end              
            end
        end
    end

end

% sum up the contributions
% A = A + S + epsilon.*I - transpose(I);
%b = b + bs + epsilon.*bi;

% solve the linear system
u = A\b;
u = reshape(u, [Np, mesh.K]);
end