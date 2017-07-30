function [A, b] = linsys(K, Np, Jdet, Jcof)
%linsys function that computes the matrix A and the vector b of the linear
%system that solves the problem

A = zeros(K*Np); %poi da mettere sparsi
I = zeros(K*Np);
S = zeros(K*Np);
b = zeros(K*Np,1);
bi = zeros(K*Np,1);
bs = zeros(K*Np,1);

[nod2, wei2, nod3, wei3] = quadrature;
[dphi] = basis_grad(nod3);

for ie = 1:K %loop on the elements
    
    for i = (ie-1)*Np+1:ie*Np %loop on the basis functions
        for j = (ie-1)*Np+1:ie*Np
            for q = 1:5 %loop on quadrature points
                A(i,j) = A(i,j) + (wei3(q)/abs(Jdet(ie)))*(dphi(:,q,j-(ie-1)*Np)'*Jcof(:,:,ie)'*Jcof(:,:,ie)*dphi(:,q,i-(ie-1)*Np));
            end
        end
    end
        
end


% sum up the contributions
A = A + S - I - transpose(I);
b = b + bs - bi;

end

