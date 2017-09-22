% Test for the convergence ratio on the unitary cubic domain

% uex = @(p) exp(p(1,:).*p(2,:).*p(3,:));
% f = @(p) -uex(p).*((p(1,:).*p(2,:)).^2+(p(1,:).*p(3,:)).^2+(p(2,:).*p(3,:)).^2);
% uex_grad = @(p) [uex(p).*p(2,:).*p(3,:); uex(p).*p(1,:).*p(3,:); uex(p).*p(1,:).*p(2,:)];

% uex = @(p) (p(1,:)).^1.5;
% f = @(p) -0.75./sqrt(p(1,:));
% uex_grad = @(p) [1.5*sqrt(p(1,:)); 0; 0];

uex = @(p) p(1,:).^3+10*p(2,:).*p(3,:).^2;
f = @(p) -20*p(2,:)-6*p(1,:);
uex_grad = @(p) [3*p(1,:).^2; 10*p(3,:).^2; 20*p(2,:).*p(3,:)];

% problem data
gd = uex;
N = 1; % degree of finite finite elements
sigma = 18; % penalty coefficient
epsilon = -1; % SIP

file_names = {'..\meshes\cube_str6.mesh'; '..\meshes\cube_str48.mesh';...
              '..\meshes\cube_str384.mesh'; '..\meshes\cube_str3072.mesh'};

n_it = 4;
err_L2 = zeros(1,n_it);
err_H10 = zeros(1,n_it);
hh = zeros(1,n_it);

for ii = 1:n_it
    mesh = MeshReader3D(file_names{ii});
    hh(ii) = sqrt(3)/(mesh.K/6)^(1/3);
    [u, x, y, z] = linsys(mesh, f, gd, N, sigma, epsilon);
    [err_L2(ii), err_H10(ii)] = errors(uex, uex_grad, u, x, y, z, mesh.E2P, mesh.K, N);
end

r_L2 = log( err_L2(1:end-1)./err_L2(2:end) )./log( hh(1:end-1)./hh(2:end) );
r_H10 = log( err_H10(1:end-1)./err_H10(2:end) )./log( hh(1:end-1)./hh(2:end) );