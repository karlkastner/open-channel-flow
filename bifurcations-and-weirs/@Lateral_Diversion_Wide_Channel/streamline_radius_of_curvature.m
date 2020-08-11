% Mon 20 Jan 10:51:02 +08 2020
function R = streamline_radius_of_curvature(obj,x,y,u,v,du_dx,du_dy,dv_dx,dv_dy)

	%function R = R(obj,x,y,u,v,du_dx,du_dy,dv_dx,dv_dy)
			if (nargin()<6)
			if (nargin()<4)
				u = obj.u(x,y);
				v = obj.v(x,y);
			end
			%J = obj.J(x,y);
			%du_dx = J(1,1);
			%du_dy = J(1,2);
			%dv_dx = J(2,1);
			%dv_dy = J(2,2);
	% TODO exploit symmetry
			du_dx = obj.du_dx(x,y);
			du_dy = obj.du_dy(x,y);
			dv_dx = obj.dv_dx(x,y);
			dv_dy = obj.dv_dy(x,y);
	%[du_dx,du_dy] = obj.grad_u(x,y);
	%[dv_dx,dv_dy] = obj.grad_v(x,y);
			end
			R = -streamline_radius_of_curvature( ...
						    u,du_dx,du_dy, ...
						    v,dv_dx,dv_dy);
			% R = obj.evalk('R',x,y);
			%R = -streamline_radius_of_curvature(u,du_dx,du_dy,v,dv_dx,dv_dy);
end

