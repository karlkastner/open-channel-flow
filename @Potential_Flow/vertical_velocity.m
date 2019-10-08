% Fri 18 May 14:08:52 CEST 2018
%% determine the vertical velocity from continuity
function v = vertical_velocity(obj,baru,C)
	h  = obj.h;
	r  = obj.streamline_radius_of_curvature();
	% TODO sign required?
	u  = obj.umag();
	v  = bend_vertical_velocity(r,h,u,baru,C);
end

