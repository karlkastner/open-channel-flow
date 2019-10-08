% Sun 20 Aug 11:15:35 CEST 2017
%% sensitivity of profile with respect to depth
function df = df_dh_(ln_z0,h,s)
	if (issym(ln_z0) || issym(h) || issym(s))
		df = (ln_z0 - log(s*h))/(h*(log(h) - 1 - ln_z0)^2);
		%df = (ln_z0 - log(hi))/(H*(log(H) - 1 - ln_z0)^2);
	else
		df = bsxfun(@times,bsxfun(@minus,ln_z0, - log(s.*h)), ...
                                   1./(h.*(log(h) - 1 - ln_z0).^2));
	end
end

