% 2020-04-16 16:47:02.694662220 +0800
function [d50i, d90i] = grain_size(obj,xi,zbi)

	if (nargin()<3)
		% pi   = 1 - 2*xi;
		pi   = 1 - xi;
	else
		zblim = obj.zblim;
		pi = 1 - (zbi-zblim(1))./(zblim(2)-zblim(1));
	end

	if (0)
		sd   = obj.sd;
		d50 = obj.d50;
		d90 = obj.d90;
		d50i = exp(norminv(pi,log(d50),log(sd)));
	else
		param = obj.dparam;
		d50i = logninv(pi,param(1),param(2));
	end

	% this is not quite true for rivers, as the distribution becomes narrower where sediment is finer/coarser
	% but it is also skewed
	% (local sd is a function of local d50)
	% i.e. d90 > d90/d50 * d50i, if d50i < d50
	% and  d90 < d90/d50 * d50i, if d50i > d50
	d90i = (obj.d90./obj.d50).*d50i;
end

