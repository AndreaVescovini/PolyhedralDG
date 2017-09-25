function [mesh] = MeshReaderGambit3D(FileName)

% function [mesh] = MeshReaderGambit3D(FileName)
% Purpose  : Read in basic grid information to build grid and build
% connectivity matrices.
% NOTE     : gambit *.neu format is assumed

Fid = fopen(FileName, 'rt');

% read intro 
for ii = 1:6 
  line = fgetl(Fid);
end

% find number of vertices and number of elements
dims = fscanf(Fid, '%d');

Nv = dims(1); Ntet = dims(2);

for ii = 1:2 
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

for ii = 1:3 
  line = fgetl(Fid);
end

% read element to node connectivity
E2V = zeros(Ntet, 4);
for ii = 1:Ntet
  line   = fgetl(Fid);
  tmpcon = sscanf(line, '%lf');
  E2V(ii,1:4) = tmpcon(4:7);
end

fclose(Fid);

mesh = struct('Nv', Nv,...
              'VX', VX,...
              'VY', VY,...
              'VZ', VZ,...
              'Ntet', Ntet,...
              'E2V', E2V);

end
