% Wed 26 Feb 21:01:47 +08 2020
function [zbi, d50n, d90i, weight, width_, xn] = wetted_cross_section(obj)
	
	% oversampled polynomial transverse cross section profile
	nn = obj.nn;
	zs = obj.zs;
	nt = length(zs);
%	xn = obj.n.x;

	width  = obj.width;

	% expansion location
	xn  = mid((0:nn)/nn);
%	xn  = 0.5*mid((0:nn)/nn);

	if (obj.widthadapt)
		width_ = obj.wetted_width(zs);
		xn     = (width_./width).*rvec(xn);
	else
		width_ = repmat(width,nt,1);
		xn     = repmat(rvec(xn),nt,1);
	end

	zbi = obj.expand(xn);

	% averaging weights
	nn = obj.nn;
	weight = 1/nn;


	% (width_./width)*(ones(1,nn)./nn);
	%weight = (width_./width)*(ones(1,nn)./nn);

	% assign the bed material, coarser grains in deeper parts of the section
	[d50n, d90i]  = obj.grain_size(xn);
end


