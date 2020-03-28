% Wed 22 Jan 19:39:02 +08 2020
function Jb   = Jb(obj,x,y)
	u     = obj.u(x,y);
	v     = obj.v(x,y);

	du_dx = obj.du_dx(x,y);
	du_dy = obj.du_dy(x,y);
	dv_dx = obj.dv_dx(x,y);
	dv_dy = obj.dv_dy(x,y);

	R     = obj.R(x,y,u,v,du_dx,du_dy,dv_dx,dv_dy);

	[dR_dx,dR_dy]    = obj.dR(x,y,u,v,du_dx,du_dy,dv_dx,dv_dy);

	%h     = obj.h;
	
	% TODO, no magic numbers
	% fs = 1681./20000;
%	fs = obj.fs;

%	dub_dx = du_dx + 1./fs*h.*(dv_dx - dR_dx.*v./R)./R;
%	dub_dy = du_dy + 1./fs*h.*(dv_dy - dR_dy.*v./R)./R;
%	dvb_dx = dv_dx - 1./fs*h.*(du_dx - dR_dx.*u./R)./R; % sign swapped!
%	dvb_dy = dv_dy - 1./fs*h.*(du_dy - dR_dy.*u./R)./R;

	beta = obj.beta;

	dub_dx = du_dx + beta.*(dv_dx - dR_dx.*v./R)./R;
	dub_dy = du_dy + beta.*(dv_dy - dR_dy.*v./R)./R;
	dvb_dx = dv_dx - beta.*(du_dx - dR_dx.*u./R)./R; % sign swapped!
	dvb_dy = dv_dy - beta.*(du_dy - dR_dy.*u./R)./R;

	if (isscalar(dub_dx))
		Jb = [dub_dx, dub_dy;
		      dvb_dx, dvb_dy];
	else
		Jb = [diag(sparse(dub_dx)),diag(sparse(dub_dy));
       	              diag(sparse(dvb_dx)),diag(sparse(dvb_dy))];
	end
end

