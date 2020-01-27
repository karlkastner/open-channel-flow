% Mon 20 Jan 10:51:02 +08 2020
function R = streamline_radius_of_curvature(obj,x,y)
	u = obj.u(x,y);
	v = obj.v(x,y);
	[du_dx,du_dy] = obj.grad_u(x,y);
	[dv_dx,dv_dy] = obj.grad_v(x,y);
	R = -streamline_radius_of_curvature(u,du_dx,du_dy,v,dv_dx,dv_dy);
end

