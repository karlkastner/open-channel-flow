% 2015-03-11 10:08:52.449062265 +0100
% Tue  4 Jul 15:45:38 CEST 2017
% Karl Kastner, Berlin
%
% prediction error of the vertical velocity profile
%
%% input :
%% U     : [nbin x nens]
%%         - values for each bin (or across section) and ensemble (or reference measurement)
%%         this are estimates estimates of the discharge or the cross sectional averaged
%%         velocity from the raw values
%%         - the profile should be limited to the effective profiling range,
%%         abobj 75-100m for a 600kHz ADCP
%%      
%% dn    : distance between HADCP bins
%% width : cross section width
%%
%% objput:
%%       sd_n : expected standard deviation for increasing profiling range
%
%%function [s_rel s_err s_dat rho res m2 u_pred fdx] = velocity_variation(U)
%% hadcp_prediction_error
%
%% TODO take scales and unscaled velocity to do combine with harmmean estimate
%
%% note: previus versions:
%%	residual was computed with respect to the predicted local velocity
%%	mse was not upscaled to cs, as profile was expected to cover entire cs
%%	finite width of cs was not considered
%% parametric estimate from moments, objliers should be filtered beforehand
%% Note that the median absolute deviation is not a good estimate,
%% because it may excludes rare events like reverse flow of floods
%% thus, the only acceptible more robust estimate would be mean absolute deviation
function [ff, f] = profile_prediction_error(f,ff,isfin,Q)
	n = size(ff.res);
%	if (isempty(obj))
%		obj = struct();
%	end
	if (isempty(f.width))
		f.width = f.dw*n(1);
	end
	if (isempty(f.width_predict'))
		f.width_predict = f.width;
	end
	if (isempty(f.inner))
		f.inner = (1:n(1))';
	end
	if (~isfield(ff,'scale') || isempty(ff.scale))
		ff.scale = 1;
	end
	if (~isfield(ff,'weight') || isempty(ff.weight))
		ff.weight = ones(n);
	end
	if (nargin() < 3 || isempty(isfin))
		isfin = true;
	end
	if (islogical(f.inner))
		ni = sum(f.inner);
	else
		ni = length(f.inner);
	end

	nw      = round(f.width_predict/f.dw);
	f.range = f.dw*(1:nw)';
	k       = ni;

	if (isfin)
		k_ = k;
		nw_ = nw;
	else
		k_ = 1e7;
		nw_ = 1e7;
	end
	
	% 1) autocorrelation function of the residual
	z      = zeros(ni,n(2));
	resflag = true;
	[acorr avar] = wautocorr(flat([ff.weight(f.inner,:);z]), ...
				 flat([ff.res(f.inner,:);z]),k,[],resflag,f.acfbias);
	if (0)
		[acorr avar] = wautocorr(ff.weight(f.inner,:), ...
				 ff.res(f.inner,:),k,[],resflag,f.acfbias);
		%avar  = mean(avar,2);
		%w = sqrt(sum(weight(inner,:).^2));
		w = rvec(Q.^2);
		avar  = wmean(ones(ni,1)*w,avar,2);
		acorr = avar/avar(1);
		%s2     = wvar(flat((weight(inner,:))),flat((res(inner,:))));
		%avar   = s2*acorr;    
	end

	ff.acorr    = acorr;

	% scale
	ff.avar     = ff.scale*avar;

	% 2) correlation coefficient
	if (f.acfbias)
		n = length(acorr);
		w = n./(n:-1:n-k+1)';
	else
		w = 1;
		%w   = (n:-1:n-k+1)'/n;
	end
	ff.rho = lsqnonlin(@(rho) w.*(acfar1(rho,k_,(0:k-1)',f.acfbias)-acorr),acorr(2));

	ff.acorr_model  = acfar1(ff.rho,k_,(0:k-1)',f.acfbias);
	% f_v.rho.^(0:ni-1)';
	ff.avar_model   = ff.avar(1)*ff.acorr_model;

	% associated correlation length
	ff.L        = -f.dw/log(ff.rho);

	% 3) error variance along the HADCP profile
	%     scale the variance of the profile to the variance of
	%     the underlying process comprising of the entire cross section
%	obj.rmse1    = sqrt(mse/ar1_var_range2(1,obj.rho,[],ni));
	ff.rmse0    = sqrt(ff.avar(1)/ar1_mse_range(1,ff.rho,nw_,1));

	% 4.) estimate the error variance for an increasing profiling range
	% s_rel = s_err/s_dat*f_correlation(rho,idx).*f_finite(n,idx).*sqrt(1./(idx-1));
	% rmse ~ f_correlation(rho,idx)*sqrt(s2_err*1./(N-1));
	% rmse ~ ( s_0 2 (W-l)/W L/l )^(1/2)
	% note : if finite sample is not applicable, set width to inf

	ff.range_rmse  = sqrt(ar1_mse_range(ff.rmse0,ff.rho,nw_,1:nw));

	% variance of the mean over the range
	ff.range_s2    = ar1_var_range2(ff.rmse0,ff.rho,nw,1:nw);

	% coefficient of determination
%	rval.R2    = 1-obj.rmse.^2/s2_dat;
	ff.range_R2  = 1-mean((cvec(ff.range_rmse)*rvec(Q)).^2,2)/var(Q);
end % profile_prediction_error

