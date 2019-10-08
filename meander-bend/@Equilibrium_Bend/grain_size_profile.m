% Wed 28 Feb 12:04:05 CET 2018
%
%% transverse (across channel) profile of the bed material grain size
%% in a river meander
% TODO why do not start with d_max at rhs?                                             
function [r_D, D, obj] = grain_size_profile(obj)
	% TODO no magic numbers
	reltol = 1e-2;
	R = obj.Rc;
	W = obj.width;

	% TODO no magic numbers
	p   = 1e-2;
	tol = 1e-4; 

	% minimum grain size
	obj.D_min = fzero(@(D) quad(obj.phi,sqrt(eps),D)-p,1e-4);

	% median grain size
	D_50      = fzero(@(D) quad(obj.phi,sqrt(eps),D)-0.5,1e-4);

	% maximum grain size
	obj.D_max = fzero(@(D) quad(obj.phi,sqrt(eps),D)-(1-p),2*D_50);

	% TODO no magic numbers
	p = obj.relaxation.gsd;

	% solve towards outer bend
	D_c = D_50;
	k = 0;
	while (1)
		% TODO iterate to achieve d(R+w/2) = d_max
		[r, D] = ode23s(@obj.dD_dr,[R, R+W/2], D_c, obj.odeset);
%		[D_min D_c D_max D(end)]*1000
	
		% ikeda, eq 29
		if ( abs(D(end) - obj.D_max) < reltol*obj.D_max )
			break;
		end
		% D_c = D_50 + 0.1*(D_max - D(end))
		D_c = D_c*((1-p) + p*obj.D_max/D(end));
		k   = k+1;
		if (k>obj.maxiter)
			warning('grain size profile did not converge');
			break;
		end
	end
	% solve towards inner bend
	[r_l, D_l] = ode23(@obj.dD_dr,[R, R-W/2],D_c);
	r_D = [flipud(r_l(2:end)); r];
	D   = [flipud(D_l(2:end)); D];

%	obj.grain_size.r = r_D;
%	obj.grain_size.D = D;
end

