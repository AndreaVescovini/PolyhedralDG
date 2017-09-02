function [J, Jcof, Jdet, trasl] = jacobians(x, y, z, r, s, t)
%[J, Jcof, Jdet, trasl] = jacobians(x, y, z, r, s, t, k)
%Functions that computes for the current element the determinant of
%the jacobian and the matrix jcof wrt the simplex (0,0,0),(1,0,0),(0,1,0),(0,0,1)
%(pag 176 of Modellistica numerica per problemi differenziali, Quarteroni).
    
Jtemp = zeros(3,4);
Jtemp(1,:) = [r s t ones(4,1)]\x;
Jtemp(2,:) = [r s t ones(4,1)]\y;
Jtemp(3,:) = [r s t ones(4,1)]\z;
    
trasl = Jtemp(:,4);
J = Jtemp(:,1:3);
Jdet = det(J);
Jcof = Jdet*inv(J');
   
end

