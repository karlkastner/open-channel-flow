% 2015-02-06 20:53:40.766971402 +0100
% Karl Kastner, Berlin

	% acceleration by gravity [m/s]
	g = 9.81;

	% length of bend (2pi is unrealistic, but the plot in rosovsky goes that far as the non dimensional theta scales with other bend properties)
	Tlim = [0 2*pi];
	% steps along cross section (17 in rosovsky)
	nn     = 170;
	% step width along bend
	dtheta = 0.004;
	% rosovski uses a different van Karman constatn in bends
	kappa = 0.5;

	% Chezy's coefficient (constant)
	C = 68;
	% radius of curvature
	Rc = 5;
	% width of the bend
	B  = 1.64;
	% grid along cross section
	y = B*((0:nn)'/nn-0.5);
	r  = Rc+y;
	hmax = 0.140;
	vmax = 0.4
	% longitudinal water level slope
	It   = 0.0015;

	% local depth
	yt = 2*y/B;
	h  = hmax*(1-yt.^2);
	% initial velocity
	v = vmax*(h/hmax).^0.4;

	% solve pde
        [theta v] = rozovskii(v,h,r,It,C,kappa,g,Tlim,dtheta);

%	theta = 1.5*sqrt(g)/(kappa^3*C)*hmax/B*theta0;
	% interpolate to fixed steps
	thetat0 = (0:4.2:21)/100;
	theta0 = thetat0 / (1.5*sqrt(g)/(kappa^3*C)*hmax/B)
	v = interp1(theta,v',theta0)';
	v
	figure(1)
	clf
	plot(2*y/B,v/vmax)
	axis equal;
	grid on;
	set(gca,'ytick',0:0.1:2)	
	figure(2);
	plot(nansum(v.*repmat(h,1,size(v,2))));



	
