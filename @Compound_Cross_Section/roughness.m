% Tue 21 Apr 08:53:36 +08 2020
function ret = roughness(Q,zs,zb,W,S,ismanning)
	Q  = cvec(Q);
	zs = cvec(zs);
%	if (isvector(zb))
%		zb = rvec(zb);
%	end
	h  = max(0,zs-zb);
	if (ismanning)
		% Q = 1/n w mean(h^(5/3)) S^(1/2)
		n = (W./Q.*mean(h.^(5/3),2).*sqrt(S));
		ret = n;
	else
		% Q = C w mean(h^(3/2)) S^(1/2)
		C = Q./(W.*mean(h.^1.5,2).*sqrt(S));
		ret = C;
	end
end

