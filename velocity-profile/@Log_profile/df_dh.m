% Sun 20 Aug 11:14:31 CEST 2017
%% sensitivity of profile with respect to depth
function [df, obj] = df_dh(obj,H,hi)
	ln_z0 = cvec(obj.ln_z0);
	
	df = obj.df_dh_(H,hi,ln_z0);
end

