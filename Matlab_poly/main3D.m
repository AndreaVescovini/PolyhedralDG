% 3D Poisson problem with Dirichlet conditions
% -lapl(u) = f    in the domain
%    u = gd   on the boundary
%
% Author: Andrea Vescovini

uex = @(p) exp(p(1,:).*p(2,:).*p(3,:));
f = @(p) -uex(p).*((p(1,:).*p(2,:)).^2+(p(1,:).*p(3,:)).^2+(p(2,:).*p(3,:)).^2);
uex_grad = @(p) [uex(p).*p(2,:).*p(3,:); uex(p).*p(1,:).*p(3,:); uex(p).*p(1,:).*p(2,:)];

% uex = @(p) sin(pi*p(1,:)).*sin(pi*p(2,:)).*sin(pi*p(3,:));
% f = @(p) 3*pi^2*sin(pi*p(1,:)).*sin(pi*p(2,:)).*sin(pi*p(3,:));
% uex_grad = @(p) [pi*cos(pi*p(1,:)).*sin(pi*p(2,:)).*sin(pi*p(3,:));
%                  pi*sin(pi*p(1,:)).*cos(pi*p(2,:)).*sin(pi*p(3,:));
%                  pi*sin(pi*p(1,:)).*sin(pi*p(2,:)).*cos(pi*p(3,:))];

% uex = @(p) p(1,:).^3+10*p(2,:).*p(3,:).^2;
% f = @(p) -20*p(2,:)-6*p(1,:);
% uex_grad = @(p) [3*p(1,:).^2; 10*p(3,:).^2; 20*p(2,:).*p(3,:)];

% uex = @(p) p(1,:).*p(2,:);
% f = @(p) 0;
% uex_grad = @(p) [p(2,:); p(1,:); 0];

% uex = @(p) p(1,:);
% f = @(p) 0;
% uex_grad = @(p) [1; 0; 0];

% problem data
gd = uex; % Dirichlet datum
N = 3; % degree of finite elements
sigma = 10; % penalty coefficient
epsilon = -1; % method

% read the mesh and generate connectivity matrices
mesh = MeshReader3Dpoly('..\meshes\cube_str48h.mesh');

% elaborate geometric information
geom = geometric_elaboration(mesh);

% build the linear system and solve it
u = linsys(geom, f, gd, N, sigma, epsilon);
export_solution(u, mesh, geom.bb);

% compute errors
[err_L2, err_H10] = errors(uex, uex_grad, u, geom, N);
