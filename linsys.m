function [A, b] = linsys( input_args )
%linsys function that computes the matrix A and the vector b of the linear
%system that solves the problem

A = zeros(K*Np); %poi da mettere sparsi
I = zeros(K*Np);
S = zeros(K*Np);
b = zeros(K*Np,1);
bi = zeros(K*Np, 1);
bs = zeros(K*Np,1);


for ie = 1:K %loop on the elements
    
    for i = (ie-1)*Np+1:ie*Np %loop on the basis functions
        for j = i:ie*Np
            A(i,j) = A(i,j) + (..)/Jdet(k);
        end
    end
        
end


% sum up the contributions
A = A + S - I - transpose(I);
b = b + bs - bi;

end

