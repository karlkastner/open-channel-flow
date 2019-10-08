% Wed 14 Mar 16:06:34 CET 2018
%
%% near-bed-velocity in a meander bend
function [ubed,vbed,ubmag,vsi] = bend_velocity_near_bed(u,v,h,R)
	kappa = Constant.Karman;
	% rozovskii 1.114 (1.5)
	% vriend 1977 119 (2.0)
	% kalkwijk 36
	scale = 2.0/kappa^2*h./R;

	% near bed velocity = offset of secondary flow + mean velocity
	% Thu 30 Aug 17:19:15 CEST 2018 sign of u and v swapped
	ubed = u + scale.*v; % == u - v*vsi
	vbed = v - scale.*u; % == v + u*vsi

	if (0)
	if (~isempty(umag) && nargout > 2)
		ubmag = (1+scale.^2).*umag;
		if (nargout() > 3)
			% secondary flow intensity
			% TODO adapation along streamline
			% 1.114 in rozovskii
			vsi = umag.*scale;
		end
	end
	end
end

