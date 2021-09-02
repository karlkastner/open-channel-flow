% Sat 23 May 21:31:44 +08 2020
% alpha : initial angle
function [x, y, z, R] = idealized_meander_bend(Rcm,alpha,T,w,n,h0,ah)
	% curvature
	ccm = 1./Rcm;

	if (nargin()<6)
		h0 = 0;
	end
	if (nargin()<7)
		ah = 0;
	end

	% TODO resample
	[xy,uv,cc] = meander_centreline(Rcm,alpha,T,n(1));

	umag = hypot(uv(:,1),uv(:,2));
	yn   = linspace(-1/2,1/2,n(2));

	% normal direction v,-u
	dx  =  uv(:,2)./umag.*w.*yn;
	dy  = -uv(:,1)./umag.*w.*yn;
	dR  = hypot(dx,dy).*sign(yn);

	x  = repmat(xy(:,1),1,n(2)) + dx;
	y  = repmat(xy(:,2),1,n(2)) + dy;

	R  = 1./cc + dR;

	% series expansion
	c  = cc - cc.^2.*dR;

	%z  = -h0*(1 + ah*R./Rc);
	z  = -h0 - ah*h0*(cc./c-1);
	% fix for c = 0;
	z(abs(c)<sqrt(eps)) = -h0;
end


