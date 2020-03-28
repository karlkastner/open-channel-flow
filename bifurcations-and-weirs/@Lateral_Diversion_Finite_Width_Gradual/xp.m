% Mon 10 Feb 18:59:11 +08 2020

function [xp,dw,np] = xp(obj)
	m = obj.m;
	switch (obj.mode)
	case {0} % constant outflow velocity along segments
		xmatch  = -1/2+(1:m)'/(m+1);
		% left and right end points
		xlr = [-1/2;mid(xmatch);1/2];
		% segment mid-points, coincide with matching locations
		% except in first and last segment
		xp  = mid(xlr);
		% segment length
		dw  = diff(xlr);
		% number of potential support points
		%np = ns;
		np = length(xp);
	case {1} % linearly varying strength along segments
if(1)
		xmatch  = -1/2+(1:m)'/(m+1);
		% left and right end points
		xlr = [-1/2;mid(xmatch);1/2];
		% segment mid-points, coincide with matching locations
		% except in first and last segment
		xp  = mid(xlr);
else
		xlr = -1/2+(0:m+1)'/(m+1);
		% segment mid-points
		xp  = mid(xlr);
		%xp  = xlr;
end
		% segment length
		dw  = diff(xlr);
		% number of potential support points
		np = length(xlr);
	end
end

