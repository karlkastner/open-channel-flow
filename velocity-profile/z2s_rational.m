% Tue  5 Jun 09:20:42 CEST 2018
% 
% z0(h1)     = h0
% z0(0)      = 0
% dz0/dz1|_0 = 1
%
function [z0,s0] = z2s_rational(z,h,h0)
	% correction factor
	% beta = 1./h1 - 1./h0
	beta = 1./h - 1./h0;

	% corrected distance to bed
	%z0 = z./(1-beta.*z);
	beta_z = bsxfun(@times,beta,z);
	z0 = z./(1- beta_z);

	% normalized
	%s0 = z0./h0;
	if (nargout() > 1)
		s0 = bsxfun(@times,z0,1./h0);
	end
end

