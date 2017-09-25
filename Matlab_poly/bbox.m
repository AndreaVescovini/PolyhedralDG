function [bb, hk] = bbox(x, y, z, E2P, K)
%[bb, hk] = bbox(x, y, z, E2P, K)
% For every polyhedral element determine the cartesian box in which it is
% contained. Computes also the diameter.

bb = zeros(3,2, K);
hk = zeros(K,1);

for ie = 1:K % loop over polyhedra
    
    index = find(E2P == ie); % find tetrahedra that compose ie
    vx = x(:, index);
    vy = y(:, index);
    vz = z(:, index);    
    vx = vx(:);
    vy = vy(:);
    vz = vz(:);
    
    bb(1,1,ie) = min(vx);
    bb(1,2,ie) = max(vx);
    bb(2,1,ie) = min(vy);
    bb(2,2,ie) = max(vy);
    bb(3,1,ie) = min(vz);
    bb(3,2,ie) = max(vz);
    
    for ii = 1:length(vx)
        for jj = ii+1:length(vx)
            dist = sqrt((vx(ii)-vx(jj))^2+(vy(ii)-vy(jj))^2+(vz(ii)-vz(jj))^2);
            if dist > hk(ie)
                hk(ie) = dist;
            end
        end
    end
end

end

