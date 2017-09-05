function [mesh] = MeshReader3D(FileName)
% [mesh] = MeshReader3D(FileName)
% Read in basic grid information to build grid and build connectivity matrices.
% *.mesh format is assumed

Fid = fopen(FileName, 'rt');

% read intro 
for i=1:5
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

for i=1:3 
  line = fgetl(Fid);
end

% read element to node connectivity
K = fscanf(Fid, '%d', 1);
EToV = fscanf(Fid, '%d', [5, K]);
EToV = EToV(1:4,:)';

fclose(Fid);

% connectivity matrices
[EToE, EToF] = tiConnect3D(EToV);

mesh = struct('Nv', Nv,...
              'VX', VX,...
              'VY', VY,...
              'VZ', VZ,...
              'K', K,...
              'EToV', EToV,...
              'EToE', EToE,...
              'EToF', EToF);

end
