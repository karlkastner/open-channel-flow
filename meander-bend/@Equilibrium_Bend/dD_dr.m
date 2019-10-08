% Wed 28 Feb 12:04:05 CET 2018
% eq 27 in ikeda
%
% change of grain size across channel
%
function dD_dr = dD_dr(obj,r,D)
	R = obj.Rc;
	W = obj.width;
	if (D < 2*obj.D_max && D>obj.D_min)
		phi           = obj.phi(D);
		if (0 == phi)
			dD_dr = 0;
		else
			dD_dr = 1./obj.phi(D)*r/(W*R);
		end
	else
		dD_dr = 0;
	end
end

