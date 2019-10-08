% wo apr  1 11:52:36 CEST 2015
% Karl Kastner, Berlin
%
%% fit the vertical velocity profile
%
function obj = fit(obj,U,S,mask,h)

	nens = size(U,2);
	if (isscalar(h))
		h = h*ones(nens,1);
	end
	if (isempty(mask))
		mask = true(size(U));
	else
		if (~islogical(mask))
			error('mask is not logical');
		end
	end
	o = ones(size(U,1),1);

	% allocate memory
	param = NaN(obj.np,nens);
	% diagonals of covariance matrices
	s2param = NaN(obj.np,nens);
	% covariances of covariance matrices
	cparam = NaN(1/2*obj.np*(obj.np-1),nens);
	% error variance
	serr2 = NaN(nens,1);
	% number of valid samples per profile
	n      = sum(mask);
	Z      = bsxfun(@times,S,rvec(h));
	ln_Z   = log(Z);
	ln_1mS = log(1-S);
	tmask  = triu(true(obj.np),+1);

        % regress profile for each ensemble
	for idx=1:nens
		if (n(idx) >= obj.np)
			% log law with dip: u = us/k ln(z) - us/k ln(z0) + a*us/k(1-y/h)
			A     = [ln_Z(mask(:,idx),idx), ...
	                                o(mask(:,idx)), ...
				 ln_1mS(mask(:,idx),idx)];
			% regress linear parameters
			param(:,idx) = A \ U(mask(:,idx),idx);
			% regression residual
			res        = A*param(:,idx) - U(mask(:,idx),idx);
			% error variance
			serr2(idx) = res'*res/(length(res)-obj.np);
			% error covariance matrix
			S2         = serr2(idx)*inv(A'*A);
			% error variance of linear parameters
			s2param(:,idx) = diag(S2);
			% error covariance of linear parameters
			cparam(:,idx) = flat(S2(tmask));
		end % if n >= np
	end % for idx

	% write back
	obj.param   = param;
	obj.s2param = s2param;
	obj.cparam  = cparam;

	% standard error
	obj.serr = sqrt(serr2);
end % fit (Log_profile_with_dip)

