function [bb] = box(x, y, z)
%[bb] = box(x, y, z)
% For every element determine the cartesian box in which it is contained

K = size(x,2);
bb = zeros(3,2, K);

for ie = 1:K
    bb(1,1,ie) = min(x(:,ie));
    bb(1,2,ie) = max(x(:,ie));
    bb(2,1,ie) = min(y(:,ie));
    bb(2,2,ie) = max(y(:,ie));
    bb(3,1,ie) = min(z(:,ie));
    bb(3,2,ie) = max(z(:,ie));
end

end

