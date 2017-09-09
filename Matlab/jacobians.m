function [trasl, J, Jcof, Jdet] = jacobians(x, y, z, r, s, t)
%[trasl, J, Jcof, Jdet] = jacobians(x, y, z, r, s, t, k)
%Functions that computes for every element the determinant of
%the jacobian and the matrix jcof wrt the simplex (0,0,0),(1,0,0),(0,1,0),(0,0,1)
%(pag 176 of Modellistica numerica per problemi differenziali, Quarteroni).

K = size(x,2);
trasl = zeros(3, K);
J = zeros(3,3,K);
Jcof = zeros(3,3,K);
Jdet = zeros(1,K);

for ie = 1:K

    Jtemp = zeros(3,4);
    Jtemp(1,:) = [r s t ones(4,1)]\x(:,ie);
    Jtemp(2,:) = [r s t ones(4,1)]\y(:,ie);
    Jtemp(3,:) = [r s t ones(4,1)]\z(:,ie);
    
    trasl(:,ie) = Jtemp(:,4);
    J(:,:,ie) = Jtemp(:,1:3);
    Jdet(ie) = det(J(:,:,ie));
    Jcof(:,:,ie) = Jdet(ie)*inv(J(:,:,ie)');
end
   
end

