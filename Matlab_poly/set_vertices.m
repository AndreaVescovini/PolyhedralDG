function [x, y, z] = set_vertices(E2V, VX, VY, VZ)
%[x, y, z] = set_vertices(E2V, VX, VY, VZ)
%Set vertices for every tetrahedron.

% vertices
r = [0; 1; 0; 0];
s = [0; 0; 1; 0];
t = [0; 0; 0; 1];

va = E2V(:,1)'; vb = E2V(:,2)'; vc = E2V(:,3)'; vd = E2V(:,4)';
x = (1-r-s-t)*VX(va)+r*VX(vb)+s*VX(vc)+t*VX(vd);
y = (1-r-s-t)*VY(va)+r*VY(vb)+s*VY(vc)+t*VY(vd);
z = (1-r-s-t)*VZ(va)+r*VZ(vb)+s*VZ(vc)+t*VZ(vd);

end