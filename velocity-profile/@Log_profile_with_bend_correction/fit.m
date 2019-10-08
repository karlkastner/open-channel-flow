% Thu Jul 17 13:16:00 WIB 2014
% Sat  7 Jan 14:24:42 CET 2017
% Karl Kastner, Berlin
%
%% fit the vertical velocity profile
% function obj = fit(obj,U,S,h,binmask)
%
function obj = fit(obj,U,S,h,binmask,ensmask)
	nbin = size(U,1);
	nens = size(U,2);
	if (isscalar(h))
		h = h*ones(nens,1);
	end
	if (nargin() < 5 || isempty(binmask))
		binmask = true(size(U));
	else
		if (~islogical(binmask))
			error('mask is not logical');
		end
	end
	if (nargin() < 6 || isempty(ensmask))
		ensmask=1:size(U,2);
	else
	if (islogical(ensmask))
		ensmask = find(ensmask);
	end
	end

	% allocate memory
	if (isempty(obj.param))
		% profile parameter
		param   = NaN(obj.np,nens);
		% covariance matrices
		S2     = NaN(obj.np,obj.np,nens);
		% covariances of covariance matrices
		%cparam  = NaN(1/2*obj.np*(obj.np-1),nens);
		% error variance
		serr  = NaN(nens,1);
	else
		param  = obj.param;
		%s2param = obj.s2param;
		S2     = obj.S2;
		%cparam = obj.s2param;
		serr   = obj.serr;
	end
	%o    = ones(nbin,1);

	% mean of ther log Z
	% mu_ln_Z  = NaN(nens,1);

	% apply maximum relative range above bottom range
	binmask = binmask & (S <= obj.smax);

	% number of valid samples per profile
	n       = sum(binmask);
	Z       = bsxfun(@times,S,rvec(h));
	%ln_Z    = log(Z);
	% mu_ln_Z = (sum(ln_Z.*mask)./n)';
	tmask  = triu(true(obj.np),+1);

        % regress profile for each ensemble
        for idx=rvec(ensmask)
		if (n(idx) >= obj.np)
			c = [0.1*mean(U(binmask(:,idx),idx));0;0];
			[c,resn,res] = lsqnonlin(@(c) obj.u_(S(binmask(:,idx)),h(idx),c) - U(binmask(:,idx),idx), c);
			%A = obj.regmtx(Z(binmask(:,idx),idx),h(idx));
			% regress linear parameters
			param(:,idx) = c;
			%A \ U(binmask(:,idx),idx);
			% regression residual
			%res = A*param(:,idx) - U(binmask(:,idx),idx);
			% root mean square residual = error variance
			n_ = length(res);
			serr2 = (res'*res)/(n_-obj.np);
			serr(idx) = sqrt(serr2);
			% error covariance matrix
    	        	%S2(:,:,idx)     = serr2*inv(A'*A);
		end % if n >= np
      	end % for idx

	% write back
	obj.param   = param;
	obj.S2      = S2;

  	% standard error
	obj.serr  =  serr;
end % Vertical_profile::fit

