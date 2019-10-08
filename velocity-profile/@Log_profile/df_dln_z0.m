% 2015-08-10 11:05:07.767262814 +0200
%% sensitivity of velocity profile with respect to roughness length
function [df, obj] = df_dln_z0(obj,H,s)
	ln_z0 = cvec(obj.ln_z0);
	df = obj.df_dln_z0_(ln_z0,H,s);
end

