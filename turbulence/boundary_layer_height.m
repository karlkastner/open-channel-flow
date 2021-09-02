% Fri 17 Apr 17:06:01 +08 2020
% boundary layer height
% c.f eh 1967
function h_t = boundary_layer_height(u,h,S,z0_t,n)
	reltol = 1e-3;
	if (nargin()<5)
		n = 5;
	end
	g     = Constant.gravity;
	kappa = Constant.Karman;
	h_t   = h;
	% TODO gn instead of pade
	k       = 0;
	maxiter = 100;
	fdx = (h>0);
	while (1)
		k      = k+1;
		ht_old = h_t;
		h_t    = u.^2./(g.*S.*(1./kappa.*(log(h_t./z0_t)-1)).^2);
		if (max(abs(h_t(fdx)-ht_old(fdx)))<reltol*mean(abs(h_t(fdx))))
			break;
		end
		if (k>maxiter)
			warning('boundary_layer_height: no convergence');
			break;
		end
		if (any(isnan(h_t(fdx))))
			warning('boundary_layer_height: NaN');
			break;
		end
	end
end

