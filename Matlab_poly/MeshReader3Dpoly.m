function [mesh] = MeshReader3Dpoly(FileName)
% [mesh] = MeshReader3D(FileName)
% Reads the mesh from the file FileName and build connectivity matrices.
% *.mesh format is assumed
% The output is a struct that contains:
% Nv: Number of vertices
% VX, VY, VZ: vectors containing coordinates of vertices
% Ntet: Total number of tetrahedra of sub-triangulations
% K: Number of polyhedral elements
% E2V: Matrix containing Tetrahedra to vertices connectivity
% E2P: Vector containing for each tetrahedron the number of the polyhedron it belongs to

Fid = fopen(FileName, 'rt');

% read intro 
for ii = 1:5
  line = fgetl(Fid);
end

% find number of vertices
Nv = fscanf(Fid, '%d', 1);

% read vertices coordinates
xyz = fscanf(Fid, '%lf', [4, Nv]);
xyz = xyz(1:3, :);
VX = xyz(1,:); VY = xyz(2,:); VZ = xyz(3,:);

% eventually translate the domain
% VX = (VX+ones(1,Nv))/2;
% VY = (VY+ones(1,Nv))/2;
% VZ = (VZ+ones(1,Nv))/2;

for ii = 1:3
  line = fgetl(Fid);
end

% read tetrahedron to vertex connectivity
Ntet = fscanf(Fid, '%d', 1);
E2V = fscanf(Fid, '%d', [5, Ntet]);
E2V = E2V(1:4,:)';

for ii = 1:3
  line = fgetl(Fid);
end

if strcmp(line,'Polyhedra')
    % read element to tetrahedron connectivity
    K = fscanf(Fid, '%d', 1);
    E2P = fscanf(Fid, '%d', Ntet);
    if min(E2P) == 0
        E2P = E2P+1;
    end
else
    K = Ntet;
    E2P = 1:Ntet;
end

fclose(Fid);

mesh = struct('Nv', Nv,...
              'VX', VX,...
              'VY', VY,...
              'VZ', VZ,...
              'Ntet', Ntet,...
              'K', K,...
              'E2V', E2V,...
              'E2P', E2P);

end
