% Mon 10 Feb 18:45:01 +08 2020
function val = evalk_(obj,fname,x,y)
	[xp,dw] = obj.xp;
	n = obj.n;
	k = (-n:n);

	alpha_dummy = 1;
	xx   = (x-xp')./dw.';
	yy   = repmat(y,[1,length(xp)])./dw.';
	val1 = sum(obj.fun.k.(fname)(flat(xx), flat(yy) ...
			        , alpha_dummy,[],obj.gamma ...
				, k ...
			       ),2);
	val = reshape(val1,size(xx));
%	val = sum(val);
end % evalk_


