% Partition tetrahedral stuctured grid into hexahedral structured grid

% read input file and elaborate the mesh
file_input = '..\meshes\cube_str48.mesh';
mesh = MeshReader3Dpoly(file_input);
[x, y, z] = set_vertices(mesh.E2V, mesh.VX, mesh.VY, mesh.VZ);
[bb, ~] = bbox(x,y,z, mesh.E2P, mesh.K);
mesh.K = mesh.Ntet/6;
lato = round((mesh.K)^(1/3));

% open output file 
[path, name_out, ~] = fileparts(file_input);
file_output = strcat(path, '\',name_out, 'p.mesh');
fo = fopen(file_output, 'wt');
fi = fopen(file_input, 'rt');

line = fgetl(fi);
while strcmp(line, 'Triangles') == 0
    fprintf(fo, strcat(line, '\n'));
    line = fgetl(fi);
end

fclose(fi);
fprintf(fo, 'Polyhedra\n%d\n', mesh.K);

for ie = 1:mesh.Ntet
    posx = 0.5*(bb(1,1,ie)+bb(1,2,ie));
    posy = 0.5*(bb(2,1,ie)+bb(2,2,ie));
    posz = 0.5*(bb(3,1,ie)+bb(3,2,ie));
    
    pol = 1+floor(posx*lato)+floor(posy*lato)*lato+floor(posz*lato)*lato*lato;
    fprintf(fo, '%d\n', pol);
end
fprintf(fo, '\nEnd\n');

fclose(fo);