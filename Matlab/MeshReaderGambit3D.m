function [mesh] = MeshReaderGambit3D(FileName)

% function [mesh] = MeshReaderGambit3D(FileName)
% Purpose  : Read in basic grid information to build grid and build
% connectivity matrices.
% NOTE     : gambit *.neu format is assumed
%
% Author: Andrea Vescovini

Fid = fopen(FileName, 'rt');

% read intro
for i=1:6
  line = fgetl(Fid);
end

% find number of vertices and number of elements
dims = fscanf(Fid, '%d');

Nv = dims(1); K = dims(2);

for i=1:2
  line = fgetl(Fid);
end

% read node coordinates
xyz = fscanf(Fid, '%lf', [4, Nv]);
xyz = xyz(2:4, :);
VX = xyz(1,:); VY = xyz(2,:); VZ = xyz(3,:);

% eventually translate the domain
VX = (VX+ones(1,Nv))/2;
VY = (VY+ones(1,Nv))/2;
VZ = (VZ+ones(1,Nv))/2;

for i=1:3
  line = fgetl(Fid);
end

% read element to node connectivity
EToV = zeros(K, 4);
for k = 1:K
  line   = fgetl(Fid);
  tmpcon = sscanf(line, '%lf');
  EToV(k,1:4) = tmpcon(4:7);
end

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
