% Mon Feb 28 01:33:50 MSK 2011
% Karl Kästner
%
%% linearised analytic solution of the swe
%
% TODO, velocity
% TODO, circular bc, non-linear (local) characteristics
function H = hw2_analytic(T, X, h0, H0)
        g = 9.81;
        nx = length(X);
	L  = X(end)-X(1);

	% celerity (characteristic velocity)
        c = nx/X(end)*sqrt(g*h0);

	% allocate space for water depth
        H = zeros(length(H0),length(T));
	l = zeros(size(X));
	r = zeros(size(X));
        for idx=1:length(T)
                t = T(idx);
		% distance the wave traveled since t0
		dx = c*t;
		% distance in indices
		d = round(nx/L*dx);
                if (d < nx)
			% shift (interpolated) wave by travelled distance
			% characteristics tavelling to the left
                        l(1:end-d) = H0(1+d:end);
			% set piece entering the domain at right end to h0
                        l(end-d+1:end) = h0;
			% characteristics travelling to the right
                        r(1+d:end) = H0(1:end-d);
			% set pience entering the domain at left end to h0
                        r(1:d) = h0;
			% update water depth: average of coming from right and coming from left
                        H(:,idx)= 0.5*(l+r);        
                else
			% distance of wave travelled further than domain
                        H(:,idx) = h0;
                end
        end
end % hw2_analytic

%		jdx_0l = 1 + round(1:length(X) - c*t);
%		jdx_0r = 1 + round(1:length(X) + c*t);
%		jdx_0l = max(1, jdx_0l);
%		jdx_0r = min(1, jdx_0r);
%		H0 = [0; H0; 0];

%		l = H0(jdx_0l)) + sum(1 == jdx_0l);
%		r = H0(jdx_0r)) + sum(nx+1 == jdx_0r);


% {		for jdx=1:nx
%			if (jdx_0l >= 1)
%				l = H0(jdx_0l,1);
%			else
%				l = 1;
%			end
%			if (jdx_0r <= length(H0))
%				r = H0(jdx_0r,1);
%			else
%				r = 1;
%			end
%			H(jdx,idx) = 0.5*(l + r);
%		end
% }
	
