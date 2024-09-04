% 2017-03-15 13:45:16.321974154 +0100
%
%% time passes and phase shifts
%% transmission and reflection coefficient depend on direction !
%% iterative (recursive) reflection and transmission
% TODO consider phase shift
% f : frequency
% L : length
function [ar_r, at_r, kdx] = sw_reflection_stepwise(h1,w1,h2,w2,L,f,n)
	% TODO use tol
	tol = 1e-4;

if (0)
	%h2_ = h1*(h2/h1)^(1/n);
	%w2_ = w1*(w2/w1)^(1/n);
	[ar_ at_] = sw_reflection(h1,w1,h2_,w2_);

	% coefficients for waves travelling into the opposit direction
	[ar_b_ at_b_] = sw_reflection(h2_,w2_,h1,w1);
%	ar_b_ = sign(ar_b)*abs(ar_b).^(1/n);
%	at_b_ = at_b.^(1/n);

else
	% forward coefficient
%	[ar at] = sw_reflection(h1,w1,h2,w2);
	% by design
%	at_ = at^(1/n)
%	at_ = 1 + 0.5*(at-1)
	% by definition
%	ar_ = 1-at_
%	s = sign(ar);
%	ar_ = s*(1-(1-abs(ar))^(1/n))
%	wc = w1*(1+ar_)/(1-ar_);
%	wc = (w1*w2)^(1/n);
%	hc = (h1*h2)^(1/n);
	wc = w1*(w2/w1)^(1/n);
	hc = h1*(h2/h1)^(1/n);

	[ar_, at_] = sw_reflection(h1,w1,hc,wc);
%	[ar_ at_] = sw_reflection(hc,wc,h2,w2)

	% backward coefficient
	[ar_b_, at_b_] = sw_reflection(hc,wc,h1,w1);
%	[ar_b_ at_b_] = sw_reflection(h2,w2,hc,wc)



%	[ar_b at_b] = sw_reflection(h2,w2,h1,w1)
%	at_b_ = at_b^(1/n)
%	ar_b_ = 1-at_b_
%pause
	
%	ar_b_ = sign(ar_b)*abs(ar_b).^(1/n);
%	s     = sign(ar_b);
%	ar_b_ = s*(1-(1-abs(ar_b))^(1/n))
%	at_b_ = at_b.^(1/n)
%pause
%	ar_ = sign(ar)*abs(ar).^(1/n);
%	at_ = at.^(1/n);

end

	% TODO use local wave number if depth differs
	g = 9.81;
	c = sqrt(g*sqrt(h1*h2));
	k = 2*pi*f/c;

	% time for passing over one segment
	dT = 1/c*L/(n-1);
	if (1==n)
		dT = 0;
	end

%	% transmission solution for n==2
	if (2==n)
		at_r = sum(at_^n * bsxfun(@times, ...
					bsxfun(@power,ar_*ar_b_,[0:5]), ...
					cos(2*pi*f*(2*dT)*(0:5))));
		at_r
	end
	at_
	at_^n
	1-(1-ar_)^n
n
pause


	% attenuation by phase shift for backward travelling waves
	if (1==n)
		p = 1;
	else
		p = cos(k*L/(n-1));
	end

	ar_r = 0;
	at_r = 0;
	kdx =0;
	rr(1,1,1,0);
	% correct transmitted wave
%	at_r = at_r/p^(n-1);
	p = 1;

	% recursive reflection
	function rr(a,ndx,d,T)
		if (T*f > 1/4)
			%warning('funnel too long');
			%return;
		end
		if (ndx > 0 && ndx <= n && abs(a) > tol)
		kdx=kdx+1;
			% do not use switch case here, as n can be 1
			if (1 == ndx)
				% at the first step
				if (d<0) % backward (out)
					% one part is transmitted out
					ar_r = ar_r + a*at_b_*cos(2*pi*f*T);
					% one part is reflected back in
					rr(ar_b_*a, ndx+1, +1, T+dT);
				else % forward (in)
					% one part is reflected out
					% this only happens at T=0, when the wave enters
					ar_r = ar_r + a*ar_*cos(2*pi*f*T);
					% one part is transmitted in
					rr(at_*a, ndx+1, +1, T+dT);
				end
				% one part is reflected out and not any more considered	
				% ar_r = ar_r + a*ar_;
				% one part is transmitted
				% rr(at_*a, ndx+1, d);
			end % if
			if (n == ndx)
				% at last step
				if (d<0) % backward
					% cannot happen, wave enters from left
					error('here');
				else % forward
					% at the last step,
					% one part is transmitted out of the domain and not any more considered
					at_r = at_r + a*at_*cos(2*pi*f*(T-(n-1)*dT));
					% one part is reflected
					rr(ar_*a, ndx-1, -1, T+dT);
				end
			end % if
			if (ndx > 1 && ndx < n)
				% at inbetween steps, nothing is relfected or transmitted out of the domain,
				% reflection and transmission is internal
				if (d<0)
					rr(at_b_*a, ndx+d,  d, T+dT);
					rr(ar_b_*a, ndx-d, -d, T+dT);
				else
					rr(at_*a, ndx+d,  d, T+dT);
					rr(ar_*a, ndx-d, -d, T+dT);
				end
			end % if
		end % no stop of recursion
	end % function rr(a,ndx)
end % sw_reflection_recursive

