% Thu 17 May 14:31:02 CEST 2018
%% apply boundary conditions
%% p*phi + (1-p)*d/db phi = rhs
%%
%% TODO, make this return the bc-struct
%%
function [p,rhs] = boundary_condition_side_outflow(x,y,h,id,Lx,Ly,wsb,qin,qout,qside)
	p   = zeros(size(x));
	rhs = zeros(size(x));
	h   = h(id);

	tol = sqrt(eps);

	% left, no flow across boundary
	fdx      = abs(x) < tol;
	p(fdx)   = 0; % set phi' = v = 0
	rhs(fdx) = 0;

	% right, no flow across boundary
	fdx      = abs(x-Lx) < tol;
	p(fdx)   = 0; % set phi' = v = 0
	rhs(fdx) = 0;
	
	% bottom, flow into main channel
	fdx      = abs(y) < tol & (x < (Lx-wsb));
	p(fdx)   = 0; % set phi' = u = u0
	rhs(fdx) = qin./h(fdx);

	% bottom, flow out of side channel 
	% TODO dx
	fdx      = abs(y) < tol & (x > (Lx-wsb));
	p(fdx)   = 0; % set phi' = u = u0
	rhs(fdx) = qside./h(fdx);

	% top, flow out of main channel
	fdx      = abs(y-Ly) < tol & (x < (Lx-wsb));
	% p(fdx)   = 0; % set phi' = u = u0;
	% rhs(fdx) = 
	p(fdx)   = 1;
	rhs(fdx) = 1;
end % boundary condition

