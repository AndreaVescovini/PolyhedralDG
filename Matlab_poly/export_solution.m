function [] = export_solution(u, mesh, bb)
% This function exports data of the solution and of the mesh into a vtk
% legacy format file called output.vtk.
% Inputs are:
% - u is the solution of the linear system
% - mesh is the mesh struct
% - bb is the polyhedron bounding box stored in geom.bb
% The solution is computed at the tetrahedra nodes, so it will be displayed
% exactly only if you use lineare finite elements (N == 1).
% The grid is saved is an unstructured_grid of tetrahedra.
% The polyhedron numbers are used to visualize a different color in each
% polyhedron.
%
% Author: Andrea Vescovini

fprintf('Saving output.vtk...');

N = 1;
Np = (N+1)*(N+2)*(N+3)/6;
blist = basis_list(N, Np);

[X, Y, Z] = set_vertices(mesh.E2V, mesh.VX, mesh.VY, mesh.VZ);
u_nod = zeros(4, mesh.Ntet);

% Computing nodal values of the solution.
for ie = 1:mesh.Ntet
    for f = 1:Np
        [valx, ~] = LegendreP(X(:,ie), blist(f,1), bb(1,:,mesh.E2P(ie)));
        [valy, ~] = LegendreP(Y(:,ie), blist(f,2), bb(2,:,mesh.E2P(ie)));
        [valz, ~] = LegendreP(Z(:,ie), blist(f,3), bb(3,:,mesh.E2P(ie)));
        u_nod(:,ie) = u_nod(:,ie) + u(f,mesh.E2P(ie))*valx.*valy.*valz;
    end
end

out = fopen('output.vtk', 'wt');

% VTK version, header and file format
fprintf(out, '# vtk DataFile Version 2.0\nVTK from Matlab\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS %d double\n', mesh.Ntet*4);

% Dataset structure
for ie = 1:mesh.Ntet
    for v = 1:4
        fprintf(out, '%f %f %f\n', X(v,ie), Y(v, ie), Z(v, ie));
    end
end
fprintf(out, 'CELLS %d %d\n', mesh.Ntet, mesh.Ntet*5);
for ie = 1:mesh.Ntet
    fprintf(out, '4 %d %d %d %d\n', [(ie-1)*4:ie*4-1]);
end
fprintf(out, 'CELL_TYPES %d\n', mesh.Ntet);
for ie = 1:mesh.Ntet
    fprintf(out, '10\n');
end

% Dataset attributes
fprintf(out, 'POINT_DATA %d\nSCALARS solution double\nLOOKUP_TABLE default\n', mesh.Ntet*4);
for ie = 1:mesh.Ntet
    for v = 1:4
        fprintf(out, '%f\n', u_nod(v,ie));
    end
end
mapped_E2P = randperm(mesh.K);
fprintf(out, 'CELL_DATA %d\nSCALARS mesh int\nLOOKUP_TABLE default\n', mesh.Ntet);
for ie = 1:mesh.Ntet
    fprintf(out, '%d\n', mapped_E2P(mesh.E2P(ie)));
end

fclose(out);
fprintf('done!\n');

end
