% Wed  8 Jan 18:20:09 +08 2020
% Wed  8 Jan 18:20:16 +08 2020
function C = z02chezy(z0,h)
	if (~issym(z0))
		g = Constant.gravity;
		kappa = Constant.Karman;
	else
		syms g kappa positive
	end
%	z0 = h./exp((C*kappa/sqrt(g))+1);
	C  = sqrt(g)./kappa.*(log(h./z0) - 1);
	if (~issym(C))
		C(C<0) = NaN;
	end
end

