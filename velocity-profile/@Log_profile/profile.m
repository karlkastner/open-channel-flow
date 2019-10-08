% Sun 20 Aug 10:46:33 CEST 2017
%% vertical profile of the streamwise velocity
function [f_v, obj] = profile(obj,S,h)
	u    = obj.u(S,h);
	ubar = obj.ubar(h);
	f_v  = u./ubar;
end

