% Sat 23 May 21:31:44 +08 2020
function [x, y, z, R] = idealized_meander_bend(Rcm,alpha,w,h0,ah,n)
	ccm = 1./Rcm;

	% TODO resample
	[xy,uv,cc] = idealized_meander_centreline(Rcm,alpha,w,n(1));

	umag=hypot(uv(:,1),uv(:,2));
	yn  = linspace(-1/2,1/2,n(2));

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

% Sat 23 May 21:31:44 +08 2020
function [xy,uv,cc] = idealized_meander_centreline(Rcm,alpha,w,nt);
	t   = linspace(0,2,nt)';
	cm  = 1./Rcm;
	cc  = (cm*sin(pi*t));
	um  = -pi*(alpha+pi/2)/cm;
	t_  = (cm*um*cos(pi*t))/pi;
	uv  = [um*cos(t_), um*sin(t_)];
	xy=[0,0;
            cumsum(mid(uv).*(t(2)-t(1)))];
end

