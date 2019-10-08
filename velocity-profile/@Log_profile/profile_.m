% So 21. Jun 11:31:35 CEST 2015
% Karl Kastner, Berlin
%
%% scale of velocity at instrument depth to depth average velocity
%% roughness length and associated standard error can change in time,
%% i.e. may be passed as vectors
%%
%% zs    : [1xn] water surface level
%% zb    : [1x1] bottom level
%% za    : [1xn] or [1x1]
%%         level of velocity measurement,
%%         i.e. level of HADCP beam bin centre, coincides with instrument level,
%%         if the HADCP is horizontally aligned
%%         only needs to be passed as vector if instrument is redeployed or
%%         becomes misaligned
%% ln_z0 : [1xn] or [1x1]
%%         natural logarithm of the roughness length
%% s     : [1xn] or [1x1]
%%         standard error of ln_z0
%% function [fz_mu fz_s fz_sp fz_bias fz_eps] = log_profile(zs,zb,za,ln_z0,s,sp,e)
function [f_v] = profile_(ln_z0,h,s)
	% zs,zb,za,ln_z0)

	% logarithm of local flow depth
	% ln_h   = log(bsxfun(@minus,zs,zb));

	% logarithm of instrument elevation above the bottom
	% this remains constant in time
	% ln_ha  = log(bsxfun(@minus,za,zb));
	
	%ln_ha_minus_ln_z0 = bsxfun(@minus,ln_ha,ln_z0);
	%ln_h_minus_ln_z0  = bsxfun(@minus,ln_h,ln_z0);
	%ln_h_minus_ln_ha  = bsxfun(@minus,ln_h,ln_ha);

	% expected value of velocity scale
	f_v  = bsxfun(@times, log(h.*s) - ln_z0, ...
			      1./(log(h) - ln_z0 - 1));

%	% error varinace of fz_mu
%	if (~isempty(s))
%		fz_s    = bsxfun(@times,  s, abs(df_dln_z0));
%	end
%	if (~isempty(sp))
%		fz_sp   = bsxfun(@times, sp, abs(df_dln_z0));
%	end
%	% perturbation, if eps_ln_z0 is known
%	% no absolute value here, e is signed
%	if (~isempty(e))
%		fz_eps = bsxfun(@times, e, df_dln_z0);
%	else
%		fz_eps = [];
%	end
end % function profile

