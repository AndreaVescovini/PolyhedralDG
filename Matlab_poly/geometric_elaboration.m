function [geom] = geometric_elaboration(mesh)
%function [geom] = geometric_elaboration(mesh)
% This function computes geometric quantities needed to build the linear
% system.
% The output struct contains:
% E2P: Vector containing for each tetrahedron the number of the polyhedron it belongs to
% K: Number of polyhedral elements
% Ntet: Total number of tetrahedra of sub-triangulations
% Nfaces: Total numberof faces of the mesh
% faces: Matrix containing faces to vertices connectivity
% faces_neig: Matrix containing for each face the number of the tetrahedra
%             that share it and the number of that face in both the tetrahedra
% area: Vector containing face areas
% normal: Vector containing face normal vector directed from the first face_neig
%         to the second face_neig
% bb: Matrix containing coordinates of the cartesian bounding box for each polyhedron
% hk: Vector containing the diameter of each polyhedron
% hmax: maximum diameter
% Fk: Vector of matrices containing the mappings from the reference tetrahedron
%     to each tetrahedron
% Jinv: Vector containing the inverse of the jacobian of each Fk
% Jdet: Vector containing the determinant of each jacobian of each Fk
%
% Author: Andrea Vescovini

% set faces of polyhedra and tretrahedra that share them
[faces, faces_neig] = read_faces(mesh);

% set area and normal E1->E2 for every face
Nfaces = size(faces,1);
area = zeros(1,Nfaces);
normal = zeros(3,Nfaces);
for e = 1:Nfaces
    [area(e), normal(:,e)] = metric2D(mesh.VX(faces(e,:)), mesh.VY(faces(e,:)), mesh.VZ(faces(e,:)), faces_neig(e,2));
end

% set vertices of tetrahedra
[x, y, z] = set_vertices(mesh.E2V, mesh.VX, mesh.VY, mesh.VZ);

% set bounding box and diameter of each polyhedron
[bb, hk] = bbox(x,y,z, mesh.E2P, mesh.K);
hmax = max(hk);

% set maps from the reference tetrahedron to physical tetrahedra
[Fk, Jinv, Jdet] = jacobians(x,y,z);

geom = struct('E2P', mesh.E2P,...
              'K', mesh.K,...
              'Ntet', mesh.Ntet,...
              'Nfaces', Nfaces,...
              'faces', faces,...
              'faces_neig', faces_neig,...
              'area', area,...
              'normal', normal,...
              'bb', bb,...
              'hk', hk,...
              'hmax', hmax,...
              'Fk', Fk,...
              'Jinv', Jinv,...
              'Jdet', Jdet);
end
