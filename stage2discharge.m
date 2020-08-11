% Wed 26 Feb 21:01:47 +08 2020
% Tue 31 Mar 15:41:24 +08 2020
%
%function [Q,C,u] = stage_discharge(h,width,S,d50_mm,d90_mm,T_C,profile,method)
function [Q,C,u] = stage2discharge(h,width,S,d50_mm,d90_mm,T_C,profile,method)
	p       = 1;
	maxiter = 100;
	reltol  = 1e-3;

	siz    = size(h);
	h      = flat(h);
	%zbn   = flat(zbn);
	S      = flat(S);
	d50_mm = flat(d50_mm);
	d90_mm = flat(d90_mm);

	rho_w   = Constant.density.quartz;
	g       = Constant.gravity;
	kappa   = Constant.Karman;

	% discharge per unit width
	q = zeros(size(h));
	% initial value, has to be coarse for karim
        C = manning2chezy(1e-2,h);
%	C = 5*ones(size(h));
	% iteratively predict the roughness
	k = 0;
	% TODO a has to be adapted for stratification
	a = 1;
	% TODO use pade
	fdx = (h>0);
	while (true)
		k     = k+1;
		q_old = q;
		switch (profile)
		case {'log'}
			z0 = chezy2z0(C,h);
			us = sqrt(g.*h.*S);
			% note: here, z0 should not be rounded off, to
			% assure convergence
			%q  = h.*us./(a*kappa).*(log(h/z0)-1);
			q  = h.*us./(a*kappa).*(z0+(log(h./z0)-1));
		case {'power'}
			z0  = chezy2z0(C,h);
			kc  = z02ks(z0);
			q   = (8.32/a.*sqrt(g*S).*h.^(5/3).*kc.^(-1/6));
			% parkers relation yields the identical result
			u   = q./h;
			us  = a./8.32.*u.*(kc./h).^(1/6);
		end
		u    = q./h;
		rms_ = rms(q(fdx)-q_old(fdx));
		if ( rms_<p*reltol*rms(q(fdx)))
			break;
		end	
		if (~isfinite(rms_))
			warning(['function returned NaN at step ' num2str(k)]);
			break;
		end
		if (k>maxiter)
			warning('no convergence');
			break;
		end
		switch (method)
		case {'karim'}
			H_d         = dune_height_karim(C,u,h,d50_mm,T_C);
			C_          = total_roughness_karim2(H_d,h,d50_mm);
		case {'rijn'}
			T           = transport_stage_rijn(d50_mm,d90_mm,h,u,T_C);
			[d_h,d_l]   = bedform_dimension_rijn(h,d50_mm,T);
			rgh         = total_roughness_rijn(d_h,d_l,d90_mm,h);
			C_          = rgh.C;
		case {'engelund-hansen-1967'}
			C_          = total_roughness_engelund_fredsoe2(u,S,h,d50_mm);
		case {'wright-parker-2003','eh-wp'}
			if (strcmp(method,'eh-wp'))
				method_ = 'engelund-hansen-1967';
			else
				method_ = method;
			end
			tau_t   = rho_w*us.^2;
			tau_rat = total2skin_stress_ratio(tau_t,d50_mm,u,h,method_);
			tau_rat = min(tau_rat,1);
			tau_s   = tau_rat.*tau_t;
			ks      = nikuradse_roughness_length(d90_mm,'parker');
			%d90_m   = 1e-3*d90_mm;
			%ks      = 3*d90_m;
			kc      = ks.*(tau_t./tau_s).^4;
			z0      = ks2z0(kc);
			% limit z0
			z0      = min(z0,0.9*exp(-1)*h);
			C_      = z02chezy(z0,h);
		otherwise
			error('here')
		end % switch
                C = p*C_ + (1-p)*C;
	end % while
	Q = width*q;
	Q(~fdx) = 0;
	C(~fdx) = 0;
	u(~fdx) = 0;

	Q=reshape(Q,siz);
	C=reshape(C,siz);
	u=reshape(u,siz);
end % function stage_discharge

function derive()
	syms q g d50_m qs a h S kc
	%solve(q/(sqrt(g*d50_m)*d50_m) - 8.32./a.*(h/d50_m).^(5/3)*sqrt(S).*(d50_m./kc).^(1/6),q)
	%solve(q/(sqrt(g*d50_m)*d50_m) - 8.32./a.*(h/d50_m).^(5/3)*sqrt(S).*(d50_m./kc).^(1/6),q)
	q_=solve( h/d50_m - (a/8.32*(q/(sqrt(g*d50_m)*d50_m))/sqrt(S)*(kc/d50_m).^(1/6)).^(3/5),q)
	vpa(simplify(q_,'ignoreanalyticconstraints',true))
	% from eq 20:
	% q_ = solve(q/(h*sqrt(g*h*S)) - 8.32/a*(h/kc)^(1/6),q), vpa(simplify(q_,'ignoreanalyticconstraints',true))
end

