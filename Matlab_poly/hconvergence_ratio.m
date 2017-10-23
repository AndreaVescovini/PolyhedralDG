function [err_L2, err_H10, r_L2, r_H10] = hconvergence_ratio(uex, f, uex_grad, N, file_names)
%function [err_L2, err_H10, r_L2, r_H10] = convergence_ratio(uex, f, uex_grad, N, file_names)
% Test for the convergence ratio with respect to h-refinement on the unitary
% cubic domain, on the meshes stored in file_names.
%
% Author: Andrea Vescovini

% problem data
gd = uex;
sigma = 10; % penalty coefficient
epsilon = -1;

n_it = length(file_names);
err_L2 = zeros(1,n_it);
err_H10 = zeros(1,n_it);
hh = zeros(1,n_it);

for ii = 1:n_it
    mesh = MeshReader3Dpoly(file_names{ii});
    geom = geometric_elaboration(mesh);
    hh(ii) = geom.hmax;
    u = linsys(geom, f, gd, N, sigma, epsilon);
    [err_L2(ii), err_H10(ii)] = errors(uex, uex_grad, u, geom, N);
end

r_L2 = log( err_L2(1:end-1)./err_L2(2:end) )./log( hh(1:end-1)./hh(2:end) );
r_H10 = log( err_H10(1:end-1)./err_H10(2:end) )./log( hh(1:end-1)./hh(2:end) );

end
