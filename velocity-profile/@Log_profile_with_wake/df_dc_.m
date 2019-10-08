% Tue 22 Aug 11:23:53 CEST 2017
%% sensitivity of velocity profile with respect to wake parameter
function df_dc = df_dc_(ln_z0,c,h,s)
	df_dc = s.*(3*s - 2)./(ln_z0 - log(h) + 1);
end


