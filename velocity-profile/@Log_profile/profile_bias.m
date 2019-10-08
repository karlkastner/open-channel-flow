% So 21. Jun 11:31:35 CEST 2015
	% TODO how to treat the bias with respect to prediction / confidence interval?
	if (~isempty(s))
		fz_bias = bsxfun(@times, s.^2, bsxfun(@plus, ...
			bsxfun(@times, 2*(ln_h_minus_ln_z0),1./ln_ha_minus_ln_z0.^3), ...
			-2./ln_ha_minus_ln_z0.^2));
	end

