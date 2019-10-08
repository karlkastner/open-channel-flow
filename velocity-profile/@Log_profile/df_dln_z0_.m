% So 21. Jun 11:31:35 CEST 2015
% 2017-08-15 14:33:03.768513930 +0200
%% sensitivity of profile with respect to roughness length
function df = df_dln_z0_(ln_z0,h,s)
%function df = df_dln_z0_(H,hi,ln_z0)
	% sensitivity of f_v w/r to ln_z0
	% df_dln_z0 = bsxfun(@times, (1-ln_h_minus_ln_ha), 1./(ln_h_minus_ln_z0-1).^2 );
	% hi = h*s;

	if (issym(ln_z0) || issym(h) || issym(s))
		df = (log(s)+1)/(log(h)-1-ln_z0)^2;
		%df = (log(hi)+1-log(H))/(log(H)-1-ln_z0)^2;
	else
		df = bsxfun(@times,log(s)+1,1./(log(h)-1-ln_z0).^2);
		%df = bsxfun(@times,bsxfun(@minus,log(hi)+1,log(H)),1./(log(H)-1-ln_z0).^2);
	end
end

