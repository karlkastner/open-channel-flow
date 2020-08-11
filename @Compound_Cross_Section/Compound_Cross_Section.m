% Tue 21 Apr 08:55:19 +08 2020
% discharge in a compound cross-section
% under uniform-stationary flow
classdef Compound_Cross_Section
	methods (Static)
		[A,h]   = area(zs,zb,width);
		[Q, qn] = discharge(zs,zb,W,C,S,ismanning);
		ret     = roughness(Q,zs,zb,W,S,ismanning);
	end
end

