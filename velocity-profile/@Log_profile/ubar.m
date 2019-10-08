% Thu 10 Aug 11:20:55 CEST 2017
%% depth averaged velocity
function [ubar, obj] = ubar(obj,h)
	us    = obj.shear_velocity();
	ln_z0 = obj.ln_z0();
	if (issym(us) || issym(ln_z0) || issym(h))
		syms kappa
		ubar = us/kappa*(log(h) - ln_z0 - 1);
	else
		if (isrow(h))
			% not cvec(h)
			ubar  = us/Constant.Karman.*(bsxfun(@minus,log(h),(rvec(ln_z0))) - 1);
		else
			ubar  = cvec(us)/Constant.Karman.*(bsxfun(@minus,log(h),(cvec(ln_z0))) - 1);
		end % else of if row
		ubar(h<0) = 0;
	end % else of if sym
end % ubar

