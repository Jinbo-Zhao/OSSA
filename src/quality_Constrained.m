function [quality] = quality_Constrained(returned_subset,standard_point)
%This is a function evaluating a set of returned system compared with a
%given constrained Optimal point. It will only be penalized by being
%greater than the constraints.

%The first element of the standard_point should be the true constrained 
% optimal system's primary response. The remianing elements of it should be
% the constrints's on all the secondary responses.
temp=returned_subset-standard_point;
temp=max(temp,0);
quality=sum(temp,2);
end