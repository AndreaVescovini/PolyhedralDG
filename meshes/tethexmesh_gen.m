% Partition hexahedral stuctured grid into hexahedral/tetrahedral structured grid
%
% Author: Andrea Vescovini

% read input file and elaborate the mesh
file_input = 'cube_str1296h.mesh';

% open output file
[~, name_out, ~] = fileparts(file_input);
file_output = strcat(name_out, 't.mesh');
fo = fopen(file_output, 'wt');
fi = fopen(file_input, 'rt');

line = fgetl(fi);
while strcmp(line, 'Tetrahedra') == 0
    fprintf(fo, strcat(line, '\n'));
    line = fgetl(fi);
end
fprintf(fo, strcat(line, '\n'));
line = fgetl(fi);
Ntet = str2double(line);
while strcmp(line, 'Polyhedra') == 0
    fprintf(fo, strcat(line, '\n'));
    line = fgetl(fi);
end
fprintf(fo, strcat(line, '\n'));
line = fgetl(fi);
NK = str2double(line);
fprintf(fo, '%d\n', NK*3.5);

ii = 1;
primo = zeros(NK/2);
for ie = 1:Ntet

    pol = fscanf(fi, '%d', 1);

    if mod(pol,2) == 0
        fprintf(fo, '%d\n', pol);
    elseif primo((pol+1)/2) == 0
        fprintf(fo, '%d\n', pol);
        primo((pol+1)/2) = 1;
    else
        fprintf(fo, '%d\n', NK+ii);
        ii = ii+1;
    end
end
fprintf(fo, '\nEnd\n');

fclose(fi);
fclose(fo);
