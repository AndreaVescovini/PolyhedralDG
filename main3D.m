% 3D Poisson problem with Dirichlet conditions
% -u'' = f    nel dominio
%    u = gd   sul bordo

% dati del problema
uex = @(p) exp(p(1,:).*p(2,:).*p(3,:));
f = @(p) -uex(p).*((p(1,:).*p(2,:)).^2+(p(1,:).*p(3,:)).^2+(p(2,:).*p(3,:)).^2);
gd = uex;
N = 1; % degree
Np = (N+1)*(N+2)*(N+3)/6; %number of points for every element
Nfaces = 4; % number of faces for every element
sigma = 10;

% generate the mesh
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGambit3D('meshes\cubeK5.neu');

% set the neighboors
[EToE, EToF] = tiConnect3D(EToV);

% compute dof
%[x,y,z] = Nodes3D(N); [r,s,t] = xyztorst(x,y,z);
[r, s, t] = set_dof_lin;
va = EToV(:,1)'; vb = EToV(:,2)'; vc = EToV(:,3)'; vd = EToV(:,4)';
x = (-(-1+r+s+t)*VX(va)+(r)*VX(vb)+(s)*VX(vc)+(t)*VX(vd));
y = (-(-1+r+s+t)*VY(va)+(r)*VY(vb)+(s)*VY(vc)+(t)*VY(vd));
z = (-(-1+r+s+t)*VZ(va)+(r)*VZ(vb)+(s)*VZ(vc)+(t)*VZ(vd));

% x = 0.5*(-(1+r+s+t)*VX(va)+(1+r)*VX(vb)+(1+s)*VX(vc)+(1+t)*VX(vd));
% y = 0.5*(-(1+r+s+t)*VY(va)+(1+r)*VY(vb)+(1+s)*VY(vc)+(1+t)*VY(vd));
% z = 0.5*(-(1+r+s+t)*VZ(va)+(1+r)*VZ(vb)+(1+s)*VZ(vc)+(1+t)*VZ(vd));

% assemble the linear system
[A, b] = linsys(f, gd, sigma, K, Np, Nfaces, x, y, z, r, s, t, EToE, EToF);

% solve the linear system
u = A\b
uex_vec = uex([x(:) y(:) z(:)]')'

err2=norm(u-uex_vec)
err_inf=norm(u-uex_vec, inf)