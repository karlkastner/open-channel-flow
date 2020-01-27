% 21:22:43.789246429 +0800
function [dR_dx,dR_dy] = dR(obj,x,y,u,v,du_dx,du_dy,dv_dx,dv_dy)
	if (nargin()<4)
		u = obj.u(x,y);
		v = obj.v(x,y);
	%J = obj.J(x,y);
	end
%		du_dx = J(1,1);
%		du_dy = J(1,2);
%		dv_dx = J(2,1);
%		dv_dy = J(2,2);
	% TODO accept R as input
	d2u_dxdx = obj.d2u_dx2(x,y);	
	d2u_dxdy = obj.d2u_dxdy(x,y);	
	d2u_dydy = obj.d2u_dy2(x,y);	
	d2v_dxdx = obj.d2v_dx2(x,y);	
	d2v_dxdy = obj.d2v_dxdy(x,y);	
	d2v_dydy = obj.d2v_dy2(x,y);	

	dR_dx = (3.*(2.*du_dx.*u + 2.*dv_dx.*v).*(u.^2 + v.^2).^(1./2))./(2.*(du_dy.*v.^2 - dv_dx.*u.^2 + (du_dx - dv_dy).*u.*v)) ...
		 - ((u.^2 + v.^2).^(3./2).*(v.^2.*d2u_dxdy - u.^2.*d2v_dxdx - 2.*du_dx.*dv_dx.*u ...
		 + 2.*du_dy.*dv_dx.*v + du_dx.*(du_dx - dv_dy).*v + dv_dx.*(du_dx - dv_dy).*u + (d2u_dxdx - d2v_dxdy).*u.*v)) ...
  			   ./(du_dy.*v.^2 - dv_dx.*u.^2 + (du_dx - dv_dy).*u.*v).^2;
	dR_dy = (3.*(2.*du_dy.*u + 2.*dv_dy.*v).*(u.^2 + v.^2).^(1./2))./(2.*(du_dy.*v.^2 - dv_dx.*u.^2 + (du_dx - dv_dy).*u.*v)) ...
		 - ((u.^2 + v.^2).^(3./2).*(v.^2.*d2u_dydy - u.^2.*d2v_dxdy ...
		- 2.*du_dy.*dv_dx.*u + 2.*du_dy.*dv_dy.*v + du_dy.*(du_dx - dv_dy).*v + dv_dy.*(du_dx - dv_dy).*u + (d2u_dxdy - d2v_dydy).*u.*v)) ...
		   ./(du_dy.*v.^2 - dv_dx.*u.^2 + (du_dx - dv_dy).*u.*v).^2;
end
