function [Fk, Jcof, Jdet] = jacobians(x, y, z, r, s, t)
%[Fk, Jcof, Jdet] = jacobians(x, y, z, r, s, t, k)
%Functions that computes for every element the map Fk, the determinant of
%the jacobian and the matrix jcof wrt the simplex (0,0,0),(1,0,0),(0,1,0),(0,0,1)
%(pag 176 of Modellistica numerica per problemi differenziali, Quarteroni).
% Fk = [J trasl]

K = size(x,2);
Fk = zeros(3,4,K);
Jcof = zeros(3,3,K);
Jdet = zeros(1,K);

for ie = 1:K

    Fk(1,:,ie) = [r s t ones(4,1)]\x(:,ie);
    Fk(2,:,ie) = [r s t ones(4,1)]\y(:,ie);
    Fk(3,:,ie) = [r s t ones(4,1)]\z(:,ie);
    
    Jdet(ie) = det(Fk(:,1:3,ie));
    Jcof(:,:,ie) = Jdet(ie)*inv(Fk(:,1:3,ie)');
end
   
end

