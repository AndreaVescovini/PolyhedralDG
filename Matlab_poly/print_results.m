res = fopen('results.txt', 'w');
fprintf(res, 'POLYNOMIAL\n\n');


uex = @(p) p(1,:).^3+10*p(2,:).*p(3,:).^2;
f = @(p) -20*p(2,:)-6*p(1,:);
uex_grad = @(p) [3*p(1,:).^2; 10*p(3,:).^2; 20*p(2,:).*p(3,:)];
gd = uex;

sigma = 10; % penalty coefficient
file_names = {'..\meshes\cube_str6.mesh'; '..\meshes\cube_str48.mesh';...
              '..\meshes\cube_str384.mesh'; '..\meshes\cube_str3072.mesh'};
n_it = 4;


for N = 2:2
    for epsilon = -1:1

        err_L2 = zeros(1,n_it);
        err_H10 = zeros(1,n_it);
        hh = zeros(1,n_it);

        for i = 1:n_it
            mesh = MeshReader3D(file_names{i});
            hh(i) = sqrt(3)/(mesh.K/6)^(1/3);
            [u, x, y, z] = linsys(mesh, f, gd, N, sigma, epsilon);
            [err_L2(i), err_H10(i)] = errors(uex, uex_grad, u, x, y, z, N);
        end
        fprintf(res, 'N: %d, eps: %d, errl2: %.4e, errh10: %.4e \n', N, epsilon, err_L2(end), err_H10(end));
        
        r_L2 = log( err_L2(1:end-1)./err_L2(2:end) )./log( hh(1:end-1)./hh(2:end) );
        r_H10 = log( err_H10(1:end-1)./err_H10(2:end) )./log( hh(1:end-1)./hh(2:end) );
        fprintf(res, 'rL2: %f, rH10: %f\n', r_L2(end), r_H10(end));
  
    end
end

fprintf(res, '\nEXPONENTIAL\n\n');
uex = @(p) exp(p(1,:) - p(2,:) + p(3,:).^2);
f = @(p) -4*uex(p).*(1+p(3,:).^2);
uex_grad = @(p) [uex(p); -uex(p); uex(p).*p(3,:)*2];
gd = uex;

for N = 1:2
    for epsilon = -1:1

        err_L2 = zeros(1,n_it);
        err_H10 = zeros(1,n_it);
        hh = zeros(1,n_it);

        for i = 1:n_it
            mesh = MeshReader3D(file_names{i});
            hh(i) = sqrt(3)/(mesh.K/6)^(1/3);
            [u, x, y, z] = linsys(mesh, f, gd, N, sigma, epsilon);
            [err_L2(i), err_H10(i)] = errors(uex, uex_grad, u, x, y, z, N);
        end
        fprintf(res, 'N: %d, eps: %d, errl2: %.4e, errh10: %.4e \n', N, epsilon, err_L2(end), err_H10(end));
        
        r_L2 = log( err_L2(1:end-1)./err_L2(2:end) )./log( hh(1:end-1)./hh(2:end) );
        r_H10 = log( err_H10(1:end-1)./err_H10(2:end) )./log( hh(1:end-1)./hh(2:end) );
        fprintf(res, 'rL2: %f, rH10: %f\n', r_L2(end), r_H10(end));
  
    end
end

fclose(res);