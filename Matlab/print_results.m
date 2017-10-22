res = fopen('results_new.txt', 'wt');
fprintf(res, 'EXPONENTIAL - STRUCTURED TETRAHEDRAL GRID\n');


uex = @(p) exp(p(1,:).*p(2,:).*p(3,:));
f = @(p) -uex(p).*((p(1,:).*p(2,:)).^2+(p(1,:).*p(3,:)).^2+(p(2,:).*p(3,:)).^2);
uex_grad = @(p) [uex(p).*p(2,:).*p(3,:); uex(p).*p(1,:).*p(3,:); uex(p).*p(1,:).*p(2,:)];
gd = uex;

sigma = 10; % penalty coefficient
epsilon = -1;
file_names = {'..\meshes\cube_str48t.mesh';...
              '..\meshes\cube_str384t.mesh';...
              '..\meshes\cube_str1296t.mesh';...
              '..\meshes\cube_str3072t.mesh'};
n_it = 4;


for N = 1:2
    err_L2 = zeros(1,n_it);
    err_H10 = zeros(1,n_it);
    hh = zeros(1,n_it);

    for i = 1:n_it
        mesh = MeshReader3D(file_names{i});
        hh(i) = sqrt(3)/(mesh.K/6)^(1/3);
        [u, x, y, z] = linsys(mesh, f, gd, N, sigma, epsilon);
        [err_L2(i), err_H10(i)] = errors(uex, uex_grad, u, x, y, z, N);
    end        
    r_L2 = log( err_L2(1:end-1)./err_L2(2:end) )./log( hh(1:end-1)./hh(2:end) );
    r_H10 = log( err_H10(1:end-1)./err_H10(2:end) )./log( hh(1:end-1)./hh(2:end) );
    
    fprintf(res,'\nN = %d, epsilon = %d\n', N, epsilon);
    fprintf(res, 'errH10: %.4e %.4e %.4e %.4e\n', err_H10);
    fprintf(res, 'rH10: %f %f %f\n', r_H10);
    fprintf(res, 'errl2: %.4e %.4e %.4e %.4e\n', err_L2);
    fprintf(res, 'rl2: %f %f %f\n', r_L2);

end
fclose(res);