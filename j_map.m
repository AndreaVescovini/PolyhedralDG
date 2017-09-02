function [j_star] = j_map(index,e)
%function [j_star] = j_map(index, e)
% function that given the function j of the element el return the number of
% the function j_star of the element ie that has the same trace over the
% face e.
 
if any(index) % true if the function j of el has non null trace over e, so el and ie share the vertex relative to j.
    j_star = find(index);
else % if j has null trace over e.
    j_star = 5-e;
end
end

