% Mon 10 Feb 18:49:38 +08 2020
function [u, v] = uv1(obj,x,y)
	u0      = obj.evalk_('u',x,y);
	v0      = obj.evalk_('v',x,y);

	cp      = obj.cp;
	[xp,dw] = obj.xp;	
	
	yp = 0;	
	dx = (x-xp')./dw';
	dy = (y-yp')./dw';
	u1 = ( dy.*v0 + dx.*u0 + 1/(pi));
	v1 = ( dx.*v0 - dy.*u0);

	[xp,dw,np]= obj.xp;

	IA = (spdiags(ones(np,1)*[0.5,0.5],0:1,np-1,np));
	IB = (spdiags((1./dw.^0)*[-1,1],0:1,np-1,np));

	u = u0*IA*cp + u1*IB*cp;
	v = v0*IA*cp + v1*IB*cp;

	u = obj.u0 + sum(u,2);
	v = sum(v,2);
end % uv1

