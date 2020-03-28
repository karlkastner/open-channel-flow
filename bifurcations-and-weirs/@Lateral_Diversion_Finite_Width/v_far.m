% Thu 16 Jan 20:06:07 +08 2020
function vf = v_far(obj,x,y)
	alpha = obj.alpha;
	gamma = obj.gamma;
	% pseudo width
	ws = 1;
	w0 = 1./gamma*ws;
	% roots	
	r1 =  (y + x*1i)/(2*w0);
	r2 =  (y - x*1i)/(2*w0);
	%c  = -a/(4*w0*pi);
	%vf =  pi*c.*(cot(pi*r1) + cot(pi*r2))
	vf =  -alpha./(4*w0).*(cot(pi*r1) + cot(pi*r2));
%	end
end
