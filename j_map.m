function [j_star] = j_map(x,y,z,ie,j,ei,e)
%function [j_star] = j_map(x,y,z,ie,j,ei,e)
% function that given the function j of the element el return the number of
% the function j_star of the element ie that has the same trace over the
% face e.
 
% if any(index) % true if the function j of el has non null trace over e, so el and ie share the vertex relative to j.
%     j_star = find(index);
% else % if j has null trace over e.
%     j_star = 5-e;
% end

for i = 1:size(x,1)
   if x(j,ei) == x(i,ie) && y(j,ei) == y(i,ie) && z(j,ei) == z(i,ie)
      j_star = i;
      return;
   end
end

j_star = 5-e;

end

