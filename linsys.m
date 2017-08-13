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
[areas, maps2D] = metric2D(x, y, z, Nfaces);
[nod2, wei2, nod3, wei3] = quadrature;
[phi, dphi, phi_bordo] = basis_evaluation(nod3, nod2, maps2D);


% gd = @(p) p(1,:).*p(2,:).*p(3,:);
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
%                 for j = 1:Np
%                     jgl = (ie-1)*Np + j;
%                     if EToE(ie, e) == ie % if true then e is a boundary face
%                 
%                     else
%                     
%                     end
%                 end
                if EToE(ie, e) == ie % boundary face
                    %bi(igl) = bi(igl) + 
                    bs(igl) = bs(igl) + sigma*wei2(q)*gd(trasl(:,ie)+J(:,:,ie)*maps2D(:,:,e)*nod2(:,q))*phi_bordo(e,q,i)*areas(e,ie);
                end              
            end       
        end
    end
        
end

bs


% sum up the contributions
A = A + S - I - transpose(I);
b = b + bs - bi;

end

