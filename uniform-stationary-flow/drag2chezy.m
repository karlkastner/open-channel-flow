% 2017-09-14 13:51:14.099671825 +0200
%% convert drag coefficient to chezy coefficient
%% g dz_s/dx + cd w u^2/h   = 0 (swe formalism)
%%      - S  +  1/C^2 U^2/H = 0 (chezy formalism)
function C = cd_to_chezy(cd)
	if (issym(cd))
		syms g
	else
		g = Constant.gravity;
	end
	C = sqrt(g./cd);
end

