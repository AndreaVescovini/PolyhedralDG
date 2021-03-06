function [blist] = basis_list(N, Np)
%[blist] = basis_list(N, Np)
% Function that lists the degrees of monomials of the Np basis functions up to
% a total degree N
%
% Author: Andrea Vescovini

blist = zeros(Np, 3);
ii = 1;
q1 = N;

while q1 >= 0
    q2 = N-q1;
    while q2 >= 0
        q3 = N-q1-q2;
        while q3 >= 0
            blist(ii,:) =[q1 q2 q3];
            q3 = q3-1;
            ii = ii+1;
        end
        q2 = q2-1;
    end
    q1 = q1-1;
end

end
