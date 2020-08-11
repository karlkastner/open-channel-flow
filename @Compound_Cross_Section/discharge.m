% Tue 21 Apr 09:10:43 +08 2020
function [Q, qn] = discharge(zs,zb,W,C,S,ismanning);
	nn = length(zb);
	if (isvector(zb))
		zb = rvec(zb);
	end
	zs = cvec(zs);
	H  = max(zs-zb,0);
	qn = normal_flow_discharge(H,W,C,S,ismanning);
	Q  = W.*mean(qn,2);
end

