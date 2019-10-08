% Wed 28 Feb 12:04:05 CET 2018
% ikeda
%% predict transverse bed profile of an equilibrium meander bend
%
function [r, h, D, f] = bed_profile(obj)
	Hc = obj.Hc;
	Rc = obj.Rc;
	W = obj.width;
	r = obj.r();

	if (isnumeric(obj.phi))
		D = obj.phi;
		Dfun = @(ri) D*ones(size(ri));
	else
		[r_D, D, obj] = obj.grain_size_profile();
		D = interp1(r_D,D,r);
		Dfun = @(ri) interp1(r,D,ri);
	end
	obj.D = D;

	% initial bed level
	h = Hc*ones(obj.nx,1);
	
	k = 0;
	% TODO use pade
	% TODO hc must be adapated for continuity of transport
	while (1)
		k = k+1;

		% velocity profile
		f = real(obj.f(h));

		ffun = @(ri) interp1(r,f,ri);
		u_bar = obj.u_bar(h);

		% TODO avoid ode solution, just integrate on grid
		% to outer bank
		[r_h_r, h_r] = ode23s(@(r,h) obj.dh_dr(r,h,Dfun,ffun,u_bar),...
						   [Rc, Rc+W/2],Hc,obj.odeset);

		% to inner bank
		[r_h_l, h_l] = ode23s(@(r,h) obj.dh_dr(r,h,Dfun,ffun,u_bar),...
						   [Rc, Rc-W/2],Hc,obj.odeset);

		% stack
		r_h = [flipud(r_h_l(2:end)); r_h_r];
		h   = [flipud(h_l(2:end)); h_r];

		% interpolate to grid points
		h = interp1(r_h,h,r);
		if (k>obj.maxiter)
			break;
		end
		% enforce given discharge by changing u_bar
%		u = f.*r./Rc*;
%		Q = 
		Qs = obj.Qs(h);
		Qs0 = obj.Qs0;
		p = (1 + max(-0.1,min(0.1,0.1*(Qs-Qs0)./Qs)));
		Hc = Hc*p;
		h  = h*p;
	end % while 1

	obj.h = h;
	obj.D = D;

%	limits(r_f)
%	obj.f_ = @(r) interp1(r_f,f,r);
	
	%h   = interp1(r_h,h,r_D);
%	if (~isnumeric(obj.phi))
%		D   = interp1(r_D,D,r);
%		%D   = interp1(r_D,D,r_h);
%	else
%		D = D*ones(obj.nx,1);
%	end
end % bed_profile

