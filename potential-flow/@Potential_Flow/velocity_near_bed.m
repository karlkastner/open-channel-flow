% Wed 14 Mar 16:06:34 CET 2018
%%
%% determine the velocity near the bed
% TODO account for adaptation of secondary flow
function [ubed, vbed] = velocity_near_bed(obj,varargin)
	%umag = obj.mag('u',varargin{:});
	u    = obj.u(varargin{:});
	v    = obj.v(varargin{:});
%	umag =  
	h    = obj.h(varargin{:});

	R = obj.streamline_radius_of_curvature(varargin{:});
	
	% 'hack' to keep 1/R finite
	% R = sign(R)+R;

	[ubed,vbed] = bend_velocity_near_bed(u,v,h,-R);

	% ubmag = hypot(ubed,vbed);

	% quick fix
	if (0)
	fdx = ~isfinite(ubmag);
	if (sum(fdx(:))>0)
		warning('fixing non-finite bed velocities')
		ubed(fdx)  = u(fdx);
		vbed(fdx)  = v(fdx);
		ubmag(fdx) = umag(fdx);
	end
	end

	%obj.ubed_ = ubed;
	%obj.vbed_ = vbed;
end % velocity_near_bed

