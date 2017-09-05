% Test for the convergence ratio on the unitary cubic domain

% uex = @(p) exp(p(1,:).*p(2,:).*p(3,:));
% f = @(p) -uex(p).*((p(1,:).*p(2,:)).^2+(p(1,:).*p(3,:)).^2+(p(2,:).*p(3,:)).^2);
% uex_grad = @(p) [uex(p).*p(2,:).*p(3,:); uex(p).*p(1,:).*p(3,:); uex(p).*p(1,:).*p(2,:)];

uex = @(p) p(1,:).^3+10*p(2,:).*p(3,:).^2;
f = @(p) -20*p(2,:)-6*p(1,:);
uex_grad = @(p) [3*p(1,:).^2; 10*p(3,:).^2; 20*p(2,:).*p(3,:)];

% problem data
gd = uex;
N = 1; % degree of finite finite elements
sigma = 10; % penalty coefficient

file_names = {'meshes\cube_str6.mesh'; 'meshes\cube_str48.mesh';...
              'meshes\cube_str384.mesh'; 'meshes\cube_str3072.mesh'};

err_L2 = zeros(1,4);
err_H10 = zeros(1,4);
hh = zeros(1,4);

for i = 1:4
    mesh = MeshReader3D(file_names{i});
    hh(i) = sqrt(3)/(mesh.K/6)^(1/3);
    [u, x, y, z] = linsys(mesh, f, gd, N, sigma);
    [err_L2(i), err_H10(i)] = errors(uex, uex_grad, u, x, y, z, N);
end

r_L2 = log( err_L2(1:3)./err_L2(2:4) )./log( hh(1:3)./hh(2:4) ); % in theory 2
r_H10 = log( err_H10(1:3)./err_H10(2:4) )./log( hh(1:3)./hh(2:4) ); % in theory 1
