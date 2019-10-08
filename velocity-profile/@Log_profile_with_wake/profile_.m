% 2017-08-22 11:53:43.124189515 +0200
%% predict velocity profile
function  f_v = profile_(ln_z0,c,h,s)
	f_v = (log(h.*s) - ln_z0 + 2*c.*h - 3*c*h.^2)./(log(h) - ln_z0 - 1);
end

