function [u,v] = lateral_outflow_finite_width(x,y,w0,n)
	y0  = w0*[ (0:-2:-(n-2))'; (2:2:n)'];
%	y00 = 0.5; R=sqrt(0+(y0 - y00).^2); A=(y0-y00)./R.^2; A=-sum(A,2); E=(cumsum(A)./sum(A)-1); loglog(abs(E(1:end/2))); 
%	x0 = -2:

	% outflow at x = 0
	dx = x(:);
	dy = bsxfun(@minus,y(:),(y0(:))');
	R  = bsxfun(@hypot,dx,dy);
	
	u = -1/pi*dx./R.^2;
	v = -1/pi*dy./R.^2;
	u = sum(u,2);
	v = sum(v,2);
	u = reshape(u,size(x));
	v = reshape(v,size(x));
end

