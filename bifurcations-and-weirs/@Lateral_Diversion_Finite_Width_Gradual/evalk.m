% Mon 20 Jan 20:08:31 +08 2020
function val = evalk(obj,fname,x,y)
	n    = obj.n;
	m    = obj.m;
	k    = (-n:n);

	c    = obj.coefficients();

	% matching locations
	xmatch  = -1/2+(1:m)'/(m+1);
	% segment end-points
	xlr = [-1/2;mid(xmatch);1/2];
	% segment mid-points
	xc   = mid(xlr);
	% segment lengths
	dw  = diff(xlr);

	val  = 0;
	for l=1:m
		val = val + c(l)*sum(obj.fun.k.(fname)((x-xc(l))./dw(l), y./dw(l) ...
			                    , obj.alpha,obj.beta,obj.gamma,k ...
					   )...
			             ,2);
	end
	%val = val/(m+1);
end % evalk

