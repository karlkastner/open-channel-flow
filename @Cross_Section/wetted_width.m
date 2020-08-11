% 2020-04-23 20:58:26.348779277 +0800

function [width_, obj] = wetted_width(obj,zs)
%	flag  = obj.flag;
	c     = obj.c;
	width = obj.width;
	nt    = length(zs);

	% points to discretize the cross section
	% here the specific discharge is computed
	% TODO allow for gauss integration
	if (obj.widthadapt)
		% dynamically limit cross section width during low flow
		% and locate integration points dynamically along the wetted cross-section
		% TODO this works only, when the cross section consists only of
		% a simple transverse slope
		switch (obj.mode)
		case ('piecewise-linear')
			cc = repmat(rvec(c),nt,1);
			cc = cc-zs;
			n  = length(c);
			x  = 0.5*(0:n-1)/(n-1);
			r  = roots_piecewise_linear(x,cc);
			r  = min(r,[],2);
			r(~isfinite(r)) = 0.5;
			width_ = 2*width*r;
		case {'polynomial'}
		switch(obj.order)
		case {'quadratic'}
			cc = repmat(rvec(c),nt,1);
			%cc(:,1) = cc(:,1) - cvec(zs);
			cc(:,1) = cc(:,1)-zs;
			r = roots2([cc(:,2),zeros(nt,1),cc(:,1)]);
			%r = real(r);
			r(abs(imag(r))>0) = NaN;
			r(r<=0) = NaN;
			% the roots are identical with swapped sign,
			% due to absence of the linear term
			r = min(r,[],2); 
			%r = max(r,[],2); 
			r = min(r,0.5);
			width_ = 2*width.*r;
		case {'quartic'}
			cc = repmat(rvec(c),nt,1);
			cc(:,1) = cc(:,1)-zs;
			%cc(:,1) = cc(:,1) - cvec(zs);
			r = roots4_reduced([cc(:,3),cc(:,2),cc(:,1)]);
			r(abs(imag(r))>0) = NaN;
			r(r<=0) = NaN;
			%r = min(r,1);
			r = min(r,[],2); 
			r = min(r,0.5);
			width_ = 2*width.*r; 
		case {'sextic'}
			cc = repmat(rvec(c),nt,1);
			%cc(:,1) = cc(:,1) - cvec(zs);
			cc(:,1) = cc(:,1)-zs;
			r = roots6_reduced([cc(:,4),cc(:,3),cc(:,2),cc(:,1)]);
			r(abs(imag(r))>0) = NaN;
			r(r<=0) = NaN;
			%r = min(r,1);
			r = min(r,[],2); 
			%r = max(r,[],2); 
			r = min(r,0.5);
			width_ = 2*width.*r;
		otherwise
			% TODO this is superfluos since the piecewise linear
			if (order < 2)
				zb_min = min(zb);
				zb_max = max(zb);
				if (abs(zb_max-zb_min)>sqrt(eps))
					width_ = (min(zb_max,zs)-zb_min)/(zb_max-zb_min)*width;
				else
					width_ = width*ones(nt,1);
				end
			else
				cc = repmat(rvec(c),nt,1);
				cc(:,1) = cc(:,1)-cvec(zs);
				r = real(roots2(fliplr(cc)));
				% choose positive root
				% TODO, what, if there are 2 positive roots?
				r(r<0) = NaN;
				% limit cross section to maximum width
				r = min(r,1);
				%r(r>(1+sqrt(eps))) = NaN;
				width_ = width*min(r,[],2);
			end
		end
		end
	else
		% do not adapt width and let points dry during low flow
		width_ = width*ones(nt,1);
	end
end

