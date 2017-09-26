function [mesh] = MeshReader3Dpoly(FileName)
% [mesh] = MeshReader3D(FileName)
% Read in basic grid information to build grid and build connectivity matrices.
% *.mesh format is assumed

Fid = fopen(FileName, 'rt');

% read intro 
for ii = 1:5
  line = fgetl(Fid);
end

% find number of vertices
Nv = fscanf(Fid, '%d', 1);

% read node coordinates
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
    E2P = E2P+1;
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
