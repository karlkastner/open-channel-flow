% Wed 14 Mar 16:06:34 CET 2018
%%
%% determine the velocity near the bed
% TODO account for adaptation of secondary flow
function [ubed, vbed] = velocity_near_bed(obj,varargin)
	u    = obj.u(varargin{:});
	v    = obj.v(varargin{:});
	h    = obj.h(varargin{:});

	R = obj.streamline_radius_of_curvature(varargin{:});
	
	[ubed,vbed] = bend_velocity_near_bed(u,v,h,-R);
end % velocity_near_bed

