% Tue Feb  3 18:26:50 CET 2015
% Karl Kastner, Berlin
%% fit the log profile to the vertical profile of the streamwise velocity
function [us, ln_z0, d, serr, aic] = regress_log_profile_(Hens,U,Z,scut,mode)
	nens = size(U,2);

	us    = NaN(nens,1,class(U));
	ln_z0 = NaN(nens,1,class(U));
	d     = NaN(nens,1,class(U));
	serr  = NaN(nens,1,class(U));
	aic   = NaN(nens,1,class(U));

	ln_Z = log(Z);
	switch (mode)
	case {'linear'}
	% for each ensemble
	for idx=1:nens
		% get valid samples
		fdx = find(isfinite(U(:,idx)) & (Z(:,idx) < scut*Hens(idx)) );
		n = length(fdx);
		if (n > 1)		
			% regress
			A     = [ln_Z(fdx,idx) ones(n,1)];
			param = A \ U(fdx,idx);
			% transform parameter
			us(idx)    =  param(1)*Constant.KAPPA;
			ln_z0(idx) = -param(2)/param(1);
			% estimate error
			res = U(fdx,idx) - A*param;
			k = length(param);
			serr(idx) = sqrt(res'*res/(n-k));
			aic(idx) = akaike_information_criterion(serr(idx),n,k);
		end % if (n > 1)
	end % for idx
	case {'displacement'}
	for idx=1:nens
		% get valid samples
		fdx = find(isfinite(U(:,idx)) & (Z(:,idx) < scut*Hens(idx)) );
		n = length(fdx);
		if (n > 1)		
			% regress
			A     = [ln_Z(fdx,idx) ones(n,1)];
			% linear estimate
			param = A \ U(fdx,idx);
			% add displacement height
			param(end+1) = 0;
			u = double(U(fdx,idx));
			z = double(Z(fdx,idx));
			param = double(param);
			% TODO limit, such that z-param(3) always positive
			f = @(param) (param(1)*log(z-param(3)) - param(2)) - u;
			param = lsqnonlin(f,param);
			% transform parameter
			us(idx)    =  param(1)*Constant.KAPPA;
			ln_z0(idx) = -param(2)/param(1);
			d(idx)     = param(3);
			% estimate error
			res = U(fdx,idx) - f(param);
			k = length(param);
			serr(idx) = sqrt(res'*res/(n-k));
			aic(idx) = akaike_information_criterion(serr(idx),n,k);
		end % if (n > 1)
	end % for idx
	otherwise
		error('unimplemented regression method');
	end
end % regress_log_profile

