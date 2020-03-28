% Fri 17 Nov 10:58:09 CET 2017
%% predefined functions to set up channel geometry
function y = swe_geometry(x,X,type,x0,y0,a,b)
	x  = x-x0;
	L  = X(2)-x0;
	switch (type)
	case {'flat'}
		y = y0*ones(size(x));
	case {'step'}
		a  = exp(a*L);
		xc = 0.5*(0+X(2));
		y  = y0*(1 + (a-1)*(x>xc));
	case {'gstep'}
		xc = 0.5*(0+X(2));
		a = exp(a*L);
		%y = y0*(1 + (a-1)*normcdf(x,x0,b));
		y = y0*a.^normcdf(x,xc,b);
	case {'step2'}
		a = exp(a*L);
		q = 1/4;
		x1 = 0+q*L;
		x2 = 0+(1-q)*L;
		y  = y0*(1 - a*((x>x1) + a*(x>x2)));
	case {'gstep2'}
		a  = exp(a*L);
		q  = 1/4;
		x1 = x0+q*L;
		x2 = x0+(1-q)*L;
		y  = y0*(a.^normcdf(x,x1,b)).*(a.^normcdf(x,x2,b));
		%y  = y0*(1 - a*((x>x1) + a*(x>x2)));
	case {'exp'}
		y = exponent(x,0,y0,a);
	case {'linear'}
		% scale at end
		a = exp(a*L);
		y = y0*(1-(1-a)*(x-0)/L);
%		y = y0+a*(x-0);
%		y = y0+y0*a/L*(x-X(1));
		%y = y0+y0*(a-1)/L*(x-X(1));
	otherwise
		error('here')
	end

	% straight reach in fromt of estuary
	y(x<0) = y0;

end

