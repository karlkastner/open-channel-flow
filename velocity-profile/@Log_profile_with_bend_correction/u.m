% 2018-09-21 14:56:17.531756026 +0800
%% streamwise velocity
function u = u(obj,S,h)
	c = obj.param;
	u = obj.u_(S,h,c);
%	u = (  c(1)*log(z) + c(2) ...
%		    + c(3)*z.*log(z) ...
%		    + c(2)*c(3)/c(1)*z);
end

