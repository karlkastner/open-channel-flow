% Wed 28 Feb 12:04:05 CET 2018
% ikeda
%% transverse profile of the bed level of an equilibrium meander bend
%% with uniform grain size
function [r, h, D, f] = bed_profile_uniform(obj)
	Hc  = obj.Hc;
	Rc  = obj.Rc;
	W   = obj.width;

	if (isnumeric(obj.phi))
		D    = obj.phi;
		Dfun = D;
	else
		[r_D, D, obj] = obj.grain_size_profile();
		Dfun = @(r) interp1(r_D,D,r);
	end

	% initial bed level
	r = obj.r();
	h = Hc*ones(obj.nx,1);

	% TODO no magic numbers
	% compute from area and discharge
	u_bar  = 1;
	
	% TODO use pade
	% TODO hc must be adapated for continuity of transport
	k = 0;
	while (1)
		k = k+1;

		% velocity profile
		f    = real(obj.f(h));
		ffun = @(ri) interp1(r,f,ri);

		% to outer bank
		[r_h_r, h_r] = ode23(@(r,h) obj.dh_dr_uniform(r,Rc,h,u_bar,Dfun,ffun),[Rc, Rc+W/2],Hc);

		% to inner bank
		[r_h_l, h_l] = ode23(@(r,h) obj.dh_dr_uniform(r,Rc,h,u_bar,Dfun,ffun),[Rc, Rc-W/2],Hc);

		% stack
		r_h = [flipud(r_h_l(2:end)); r_h_r];
		h   = [flipud(h_l(2:end)); h_r];

		% interpolate to grid points
		h = interp1(r_h,h,r);
		if (k>obj.maxiter)
			break;
		end
	end % while 1

	obj.h = h;
	
	if (nargout() > 2)
	if (~isnumeric(obj.phi))
		D   = interp1(r_D,D,r);
	else
		D = D*ones(nx,1);
	end
	end

	if (nargout() > 3)
		obj.f_ = ffun(r); %@(r) interp1(r_f,f,r);
	end
end % bed_profile_uniform

