% Test for the convergence ratio on the unitary cubic domain
%
% Author: Andrea Vescovini

% uex = @(p) exp(p(1,:).*p(2,:).*p(3,:));
% f = @(p) -uex(p).*((p(1,:).*p(2,:)).^2+(p(1,:).*p(3,:)).^2+(p(2,:).*p(3,:)).^2);
% uex_grad = @(p) [uex(p).*p(2,:).*p(3,:); uex(p).*p(1,:).*p(3,:); uex(p).*p(1,:).*p(2,:)];

uex = @(p) p(1,:).^3+10*p(2,:).*p(3,:).^2;
f = @(p) -20*p(2,:)-6*p(1,:);
uex_grad = @(p) [3*p(1,:).^2; 10*p(3,:).^2; 20*p(2,:).*p(3,:)];

% uex = @(p) (p(1,:)-0.5).*heaviside(p(1,:)-0.5);
% f = @(p) 0;
% uex_grad = @(p) [heaviside(p(1,:) -0.5); 0; 0];

% problem data
gd = uex;
N = 2; % degree of finite finite elements
sigma = 18; % penalty coefficient
epsilon = -1; % SIP

file_names = {'..\meshes\cube_str6.mesh'; '..\meshes\cube_str48.mesh';...
              '..\meshes\cube_str384.mesh'; '..\meshes\cube_str3072.mesh'};

n_it = 3;
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
