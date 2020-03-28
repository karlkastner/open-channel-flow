% Sun  9 Feb 19:17:08 +08 2020
function [u,v] = uv_side_branch(obj,x,y)
	v0 = -obj.alpha;
	
	ns = obj.m;
	
	cs_ = obj.cs;


	ns_ = length(cs_);
if (obj.lmode)
	l = 2*pi*(1:ns_);
	% caesaro weights
	if (obj.cflag)
	w_ = (1-l/(pi*(2*ns_)))
	else
	w_ = ones(size(l));
	end
else
	l = pi*(1:ns_);
	w_ = ones(size(l));
end

if (0)
	cs = reshape(cs_,[],2);
	u = 0 + sum( ...
		    (    l.*cos(x*l).*rvec(cs(:,1)) ... 
                      + -l.*sin(x*l).*rvec(cs(:,2))).*exp(y*l), 2);
	v = v0+ sum( ...
		    (   -l.*sin(x*l).*rvec(cs(:,1)) ...
                      + -l.*cos(x*l).*rvec(cs(:,2))).*exp(y*l), 2);
else
	cs = cs_;
	u = 0 + sum( ...
		    (    w_.*l.*cos((x)*l).*rvec(cs(:,1))).*exp(y*l), 2);
		    %(    w_.*l.*cos((x+ox)*l).*rvec(cs(:,1))).*exp(y*l), 2);
	v = v0+ sum( ...
		    (    w_.*l.*sin((x)*l).*rvec(cs(:,1)).*exp(y*l)), 2);
		    %(    w_.*l.*sin((x+ox)*l).*rvec(cs(:,1)).*exp(y*l)), 2);
	
end
	fdx = abs(x)>1/2;
	u(fdx) = NaN;
	v(fdx) = NaN;	
end
