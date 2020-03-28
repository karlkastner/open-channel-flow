% Sun 20 May 11:57:21 CEST 2018
%% TODO, g should be replaced by gx,gy,gz, see chaudhri
%% assemble discretisation matrix for stationary flow
function [A,b,obj] = assemble_stationary(obj,h,qx,qy)
% 	continuity
%	Dx*qx + Dy*qy = 0
%	momentum
%	Dx*(qx.^2./h) + Dx*(qx.*qy./h) + g*h*Dx*h + g/C^2*qmax.*qx./h^2 = g*h*Dx*zb
%	Dy*(qy.^2./h) + Dy*(qx.*qy./h) + g*h*Dx*h + g/C^2*qmag.*qy./h^2 = g*h*Dy*zb

	g = obj.g;
	n = prod(obj.mesh.n);
	zb = obj.zb;
	C  = obj.C;
	Dx = obj.mesh.Dx;
	Dy = obj.mesh.Dy;

	Z = spzeros(n);
	H = diag(sparse(h));
	ux = qx./h;
	Ux = diag(sparse(ux));
	uy = qy./h;
	Uy = diag(sparse(uy));
	% qmag = hypot(qx,qy);
	umag = hypot(ux,uy);
	% bed shear stress, to be scaled by qx and qy
	tau_b = -diag(sparse(g./C.^2.*umag./h));
	% e = 1e-2;
	e = 0;
	
	if (1)
	A = [      Z,                  Dx,                    Dy;
	      g*H*Dx, Dx*Ux+Dy*Uy - tau_b + e*(Dx*Dx),                     Z;
              g*H*Dy,                   Z, Dy*Uy + Dx*Ux - tau_b + e*(Dy*Dy)]; 	   
	else
		Ux2 = diag(sparse(ux.^2));
		Uy2 = diag(sparse(uy.^2));
	A = [      Z,                  Dx,                    Dy;
	      g*H*Dx-2*Ux2*Dx, 2*Ux*Dx + Dy*Uy - tau_b + e*(Dx*Dx),                   Z;
              g*H*Dy-2*Uy2*Dy,                       Z, 2*Uy*Dy + Dx*Ux - tau_b + e*(Dy*Dy)]; 	   
	end
	M = averaging_matrix_2(obj.mesh.n);
	% lax-friedrich like stabilization
	MM = [M,Z,Z; Z,M,Z; Z,Z,M];
	p = 0.125;
	A = (1-p)*MM + p*A;
	b = [ zeros(n,1);
              g*H*Dx*zb;
              g*H*Dy*zb];
end % assemble_stationary

