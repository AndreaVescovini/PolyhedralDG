function [x, y, z] = set_vertices(EToV, VX, VY, VZ)
%[x, y, z] = set_vertices(EToV, VX, VY, VZ)
%Set vertices for every element.

% vertices
r = [0; 1; 0; 0];
s = [0; 0; 1; 0];
t = [0; 0; 0; 1];

va = EToV(:,1)'; vb = EToV(:,2)'; vc = EToV(:,3)'; vd = EToV(:,4)';
x = (1-r-s-t)*VX(va)+r*VX(vb)+s*VX(vc)+t*VX(vd);
y = (1-r-s-t)*VY(va)+r*VY(vb)+s*VY(vc)+t*VY(vd);
z = (1-r-s-t)*VZ(va)+r*VZ(vb)+s*VZ(vc)+t*VZ(vd);

end