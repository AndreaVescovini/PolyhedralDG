function [bb] = box(x, y, z, E2P, K)
%[bb] = box(x, y, z, E2P, K)
% For every polyhedral element determine the cartesian box in which it is
% contained.

bb = zeros(3,2, K);

for ie = 1:K
    
    index = find(E2P == ie);
    vx = x(:, index);
    vy = y(:, index);
    vz = z(:, index);    
    
    bb(1,1,ie) = min(vx(:));
    bb(1,2,ie) = max(vx(:));
    bb(2,1,ie) = min(vy(:));
    bb(2,2,ie) = max(vy(:));
    bb(3,1,ie) = min(vz(:));
    bb(3,2,ie) = max(vz(:));
end

end

