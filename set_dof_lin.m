function [r, s, t] = set_dof_lin()
%[r, s, t] = set_dof_lin()
%Set degrees of freedom in the reference simplex for linear fe.

r = [0; 1; 0; 0];
s = [0; 0; 1; 0];
t = [0; 0; 0; 1];

end

