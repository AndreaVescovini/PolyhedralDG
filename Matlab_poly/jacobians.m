function [Fk, Jinv, Jdet] = jacobians(x, y, z)
%[Fk, Jcof, Jdet] = jacobians(x, y, z)
%Functions that computes for every tetrahedral element the map Fk from the
%reference simplex (0,0,0),(1,0,0),(0,1,0),(0,0,1); the determinant of
%the jacobian Jdet and the inverse of the jacobian Jinv.
%(pag 176 of Modellistica numerica per problemi differenziali, Quarteroni).
% Fk = [J trasl]
%
% Author: Andrea Vescovini

K = size(x,2);
Fk = zeros(3,4,K);
Jinv = zeros(3,3,K);
Jdet = zeros(1,K);

% vertices in the reference tetrahedron plus ones(4,1)
Mat = [0 0 0 1;
       1 0 0 1;
       0 1 0 1;
       0 0 1 1];

for ie = 1:K %loop over tetrahedra

    Fk(1,:,ie) = Mat\x(:,ie);
    Fk(2,:,ie) = Mat\y(:,ie);
    Fk(3,:,ie) = Mat\z(:,ie);

    Jdet(ie) = det(Fk(:,1:3,ie));
    Jinv(:,:,ie) = inv(Fk(:,1:3,ie));
end

end
