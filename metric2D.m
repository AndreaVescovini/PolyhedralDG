function [areas, maps2D] = metric2D(x, y, z, Nfaces)
%metric2D. Function that computes for every face of every element the area
%doubled and the mappings from the 2D refence triangle to the faces of the
%3D reference tetrahedron.

K = size(x,2);
maps2D = zeros(3,3,Nfaces);

maps2D(:,:,1) = [0 1 0;
                 1 0 0;
                 0 0 0];
             
maps2D(:,:,2) = [1 0 0;
                 0 0 0;
                 0 1 0];
             
maps2D(:,:,3) = [0 0 0;
                 0 1 0;
                 1 0 0];
             
maps2D(:,:,4) = [1  0 0;
                 0  1 0;
                -1 -1 1];

areas = zeros(Nfaces,K); % sfruttando EToE e EToF potrei non calcolare le aree che così calcolo due volte

for k = 1:K        
        areas(1,k) = norm(cross([x(1,k)-x(2,k) y(1,k)-y(2,k) z(1,k)-z(2,k)],[x(1,k)-x(3,k) y(1,k)-y(3,k) z(1,k)-z(3,k)]));
        areas(2,k) = norm(cross([x(2,k)-x(1,k) y(2,k)-y(1,k) z(2,k)-z(1,k)],[x(2,k)-x(4,k) y(2,k)-y(4,k) z(2,k)-z(4,k)]));
        areas(3,k) = norm(cross([x(3,k)-x(1,k) y(3,k)-y(1,k) z(3,k)-z(1,k)],[x(3,k)-x(4,k) y(3,k)-y(4,k) z(3,k)-z(4,k)]));
        areas(4,k) = norm(cross([x(4,k)-x(2,k) y(4,k)-y(2,k) z(4,k)-z(2,k)],[x(4,k)-x(3,k) y(4,k)-y(3,k) z(4,k)-z(3,k)]));    
end

end

