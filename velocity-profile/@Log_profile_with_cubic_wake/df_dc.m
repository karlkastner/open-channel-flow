% Tue 22 Aug 11:24:22 CEST 2017
%% sensitivity of profile with respect to wave parameter
function [df_dc, obj] = df_dc(obj,h,s)
	ln_z0 = obj.ln_z0;
	df_dc = obj.df_dc_(ln_z0,h,s);
end

