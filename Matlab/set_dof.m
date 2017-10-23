function [r, s, t] = set_dof(N)
%[r, s, t] = set_dof(N)
%Set degrees of freedom in the reference simplex.
%
% Author: Andrea Vescovini

% vertices must be the first 4 points
if N == 1
    r = [0; 1; 0; 0];
    s = [0; 0; 1; 0];
    t = [0; 0; 0; 1];
elseif N == 2
    r = [0; 1; 0; 0; 0.5;  0;  0; 0.5; 0.5;   0];
    s = [0; 0; 1; 0;  0;  0.5; 0; 0.5;  0;  0.5];
    t = [0; 0; 0; 1;  0;   0; 0.5; 0;  0.5; 0.5];

else
    disp('Wrong order');
end
