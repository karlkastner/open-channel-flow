% Wed 22 Apr 16:13:47 +08 2020
% Karl Kastner
function [A,h] = area(zs,zb,width)
	if (isvector(zb))
		zb = rvec(zb);
	end
	if (isvector(zs))
		zs = cvec(zs);
	end
	h   = zs-zb;
	h   = max(h,0);
	A   = width.*mean(h,2);
end

