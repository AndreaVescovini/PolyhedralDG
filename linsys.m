function [A, b] = linsys(f, gd, sigma, K, Np, Nfaces, J, Jcof, Jdet, trasl, EToE, EToF)
%linsys function that computes the matrix A and the vector b of the linear
%system that solves the problem

A = zeros(K*Np); %poi da mettere sparsi
I = zeros(K*Np);
S = zeros(K*Np);
b = zeros(K*Np,1);
bi = zeros(K*Np,1);
bs = zeros(K*Np,1);

[nod2, wei2, nod3, wei3] = quadrature;
[phi, dphi] = basis_evaluation(nod3);

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
                
                    else
                    
                    end
                end
                if EToE(ie, e) == ie % boundary face
                    bi(igl) = bi(igl) + 
                    bs(igl) = bs(igl) + sigma*wei2(q)*
                end              
            end       
        end
    end
        
end


% sum up the contributions
A = A + S - I - transpose(I);
b = b + bs - bi;

end

