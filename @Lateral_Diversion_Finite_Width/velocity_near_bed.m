% Mon 20 Jan 10:51:16 +08 2020
function [ub,vb] = velocity_near_bed(obj,x,y)
%	h        = obj.h;
	[u,v]    = obj.velocity(x,y);
	R        = obj.streamline_radius_of_curvature(x,y);
	[ub,vb]  = bend_velocity_near_bed(u,v,obj.beta/12,R);
end

