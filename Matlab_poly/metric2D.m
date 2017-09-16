function [areas, normals] = metric2D(x, y, z, Nfaces)
%[areas, normals] = metric2D(x, y, z, Nfaces).
%Function that computes for every face of every element the area doubled
%and the normal vector.

areas = zeros(Nfaces); % sfruttando EToE e EToF potrei non calcolare le aree che così calcolo due volte
normals = zeros(3,Nfaces);


%face 1, clockwise numeration
areas(1) = norm(cross([x(1)-x(2) y(1)-y(2) z(1)-z(2)],[x(3)-x(2) y(3)-y(2) z(3)-z(2)]));
normals(:,1) = cross([x(1)-x(2) y(1)-y(2) z(1)-z(2)],[x(3)-x(2) y(3)-y(2) z(3)-z(2)])./areas(1);
%face 2, counter-clockwise numeration
areas(2) = norm(cross([x(2)-x(1) y(2)-y(1) z(2)-z(1)],[x(4)-x(1) y(4)-y(1) z(4)-z(1)]));
normals(:,2) = cross([x(2)-x(1) y(2)-y(1) z(2)-z(1)],[x(4)-x(1) y(4)-y(1) z(4)-z(1)])./areas(2);
%face 3, clockwise numeration
areas(3) = norm(cross([x(4)-x(1) y(4)-y(1) z(4)-z(1)],[x(3)-x(1) y(3)-y(1) z(3)-z(1)]));
normals(:,3) = cross([x(4)-x(1) y(4)-y(1) z(4)-z(1)],[x(3)-x(1) y(3)-y(1) z(3)-z(1)])./areas(3);
%face 4, counter-clockwise numeration
areas(4) = norm(cross([x(4)-x(3) y(4)-y(3) z(4)-z(3)],[x(2)-x(3) y(2)-y(3) z(2)-z(3)]));
normals(:,4) = cross([x(4)-x(3) y(4)-y(3) z(4)-z(3)],[x(2)-x(3) y(2)-y(3) z(2)-z(3)])./areas(4);

end

