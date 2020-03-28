% Fri 18 May 18:45:42 CEST 2018
%% compute the sediment transport
function [qsx,qsy] = sediment_transport(obj,h)
	u = obj.u;
	v = obj.v;

	[ubed, vbed] = obj.velocity_near_bed();

	% rigid lid
	if (nargin() < 2)
		zb = -obj.h;
	else
		zb = -h;
	end
	% bed slope
	dzb_dx = obj.mesh.Dx*zb;
	dzb_dy = obj.mesh.Dy*zb;
	% sediment transport
	[qsx,qsy] = sediment_transport_directed(u,v,ubed,vbed,obj.h,...
							dzb_dx,dzb_dy, ...
							obj.chezy,obj.d_mm);
end % sediment_transport

