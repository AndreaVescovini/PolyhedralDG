% This script prints over a file the results of convergence tests.
%
% Author: Andrea Vescovini

res = fopen('results_gen.txt', 'wt');
fprintf(res, 'EXPONENTIAL - POLYHEDRAL GRID\n');

uex = @(p) exp(p(1,:).*p(2,:).*p(3,:));
f = @(p) -uex(p).*((p(1,:).*p(2,:)).^2+(p(1,:).*p(3,:)).^2+(p(2,:).*p(3,:)).^2);
uex_grad = @(p) [uex(p).*p(2,:).*p(3,:); uex(p).*p(1,:).*p(3,:); uex(p).*p(1,:).*p(2,:)];

file_names = {'..\meshes\cube_str48p.mesh';...
              '..\meshes\cube_str384p.mesh';...
              '..\meshes\cube_str1296p.mesh';...
              '..\meshes\cube_str3072p.mesh'};
for N = 1:3
    [err_L2, err_H10, r_L2, r_H10] = hconvergence_ratio(uex, f, uex_grad, N, file_names);
    fprintf(res,'\nN = %d\n', N);
    fprintf(res, 'errH10: %.4e %.4e %.4e %.4e\n', err_H10);
    fprintf(res, 'rH10: %f %f %f\n', r_H10);
    fprintf(res, 'errl2: %.4e %.4e %.4e %.4e\n', err_L2);
    fprintf(res, 'rl2: %f %f %f\n', r_L2);
end

fclose(res);
