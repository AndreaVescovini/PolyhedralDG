function [A, b] = linsys(f, gd, sigma, K, Np, Nfaces, x, y, z, r, s, t, EToE, EToF)
%linsys function that computes the matrix A and the vector b of the linear
%system that solves the problem

A = zeros(K*Np); %poi da mettere sparsi
I = zeros(K*Np);
S = zeros(K*Np);
b = zeros(K*Np,1);
bi = zeros(K*Np,1);
bs = zeros(K*Np,1);

% compute the jacobians
%%%%%% BISOGNA CONTROLLARE CHE SE N>1 A ME INTERESSANO SOLO I VERTICI,
%%%%%% DIPENDE DA COME CREO PRIMA r s t, CHE TRAMITE EToV MAPPANO x y z,
%%%%%% POTREI PER ESEMPIO METTERE I VERTICI SUBITO NEI PRIMI NODI.
[J, Jcof, Jdet, trasl] = jacobians(x,y,z,r,s,t); % del determinante dovro' poi prenderne il valore assoluto
[areas, normals, maps2D] = metric2D(x, y, z, Nfaces);
[nod2, wei2, nod3, wei3] = quadrature;
[phi, dphi, phi_bordo, grad_bordo] = basis_evaluation(nod3, nod2, maps2D);

for ie = 1:K %loop on the elements
    for q = 1:5 %loop on 3D quadrature points
        for i = 1:Np %loop on the basis functions
            igl = (ie-1)*Np + i; %global node number
            for j = 1:Np
                jgl = (ie-1)*Np + j;
                A(igl,jgl) = A(igl,jgl) + (wei3(q)/abs(Jdet(ie)))*(dphi(:,q,j)'*Jcof(:,:,ie)'*Jcof(:,:,ie)*dphi(:,q,i));                                
            end
            b(igl) = b(igl) + abs(Jdet(ie))*wei3(q)*f(trasl(:,ie)+J(:,:,ie)*nod3(:,q))*phi(q,i);
        end
    end
    
    for e = 1:Nfaces
        for q = 1:4 %loop on 2D quadrature nodes
            for i = 1:Np
                igl = (ie-1)*Np + i;
                for j = 1:Np
                    jgl = (ie-1)*Np + j;
                    if EToE(ie, e) == ie % if true then e is a boundary face
                        S(igl, jgl) = S(igl, jgl) + sigma*wei2(q)*phi_bordo(q, e, i)*phi_bordo(q, e, j)*areas(e, ie);
                        I(igl, jgl) = I(igl, jgl) + wei2(q)*(normals(:,e,ie)'*Jcof(:,:,ie)*grad_bordo(:,q,e,i)/Jdet(ie))*phi_bordo(q,e,j)*areas(e,ie);      
                    else
                        jneigh = (EToE(ie, e)-1)*Np + j; %column index corresponding to the neighbooring element
                        S(igl, jgl) = S(igl, jgl) + sigma*wei2(q)*phi_bordo(q, e, i)*phi_bordo(q, e, j)*areas(e, ie);
                        S(igl, jneigh) = S(igl, jneigh) - sigma*wei2(q)*phi_bordo(q, e, i)*phi_bordo(q, EToF(ie, e), j)*areas(e, ie);
                        I(igl, jgl) = I(igl, jgl) + 0.5*wei2(q)*(normals(:,e,ie)'*Jcof(:,:,ie)*grad_bordo(:,q,e,i)/Jdet(ie))*phi_bordo(q,e,j)*areas(e,ie);
                        I(igl, jneigh) = I(igl, jneigh) - 0.5*wei2(q)*(normals(:,e,ie)'*Jcof(:,:,ie)*grad_bordo(:,q,e,i)/Jdet(ie))*phi_bordo(q,EToF(ie,e),j)*areas(e,ie);  
                    end
                end
                if EToE(ie, e) == ie % boundary face
                    % valuto il dato al bordo gd nei nodi fisici
                    bs(igl) = bs(igl) + sigma*wei2(q)*gd(trasl(:,ie)+J(:,:,ie)*maps2D(:,:,e)*nod2(:,q))*phi_bordo(q,e,i)*areas(e,ie);
                    bi(igl) = bi(igl) + wei2(q)*gd(trasl(:,ie)+J(:,:,ie)*maps2D(:,:,e)*nod2(:,q))*(normals(:,e,ie)'*Jcof(:,:,ie)*grad_bordo(:,q,e,i)/Jdet(ie))*areas(e,ie);
                end              
            end       
        end
    end

end

% sum up the contributions
A = A + S - I - transpose(I);
b = b + bs - bi;

end

