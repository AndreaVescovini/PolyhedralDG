function [area, normal] = metric2D(x, y, z, Nface)
%[areas, normals] = metric2D(x, y, z, Nfaces).
%Function that computes the area doubled and the normal vector of the face
%described by the vertex in x, y, z. Nface indicate the local number of the
%face and it affects the orientation of the vertices.

if (Nface == 1 || Nface == 3) % clockwise numeration
    normal = cross([x(1)-x(2); y(1)-y(2); z(1)-z(2)],[x(3)-x(2); y(3)-y(2); z(3)-z(2)]);
else % counter-clockwise numeration
    normal = cross([x(3)-x(2); y(3)-y(2); z(3)-z(2)],[x(1)-x(2); y(1)-y(2); z(1)-z(2)]);  
end

area = norm(normal);
normal = normal./area;

end