% 2015-02-06 17:48:41.315560290 +0100
% Karl Kastner, Berlin
%
%% transversal velocity distribution in a bend
%% Rososkii,
%% as in the book central differences along the radius and euler forward in space
%% are used, note that since the advent of the computer more advanced schemes
%% could be used (see build in solvers)
%% cfl condition is not explicitely checked
%% Rosovsky assumes a constant water level, e.g. does not consider superelevation
%%
%% I_theta = -1/r dz/d_theta	(p. 22)
%% d_theta = 1/R ds|_R
%% => I_theta = -R/r dz/ds = -R/r I_0
%
%% It : (1.32) drop of level per unit angle, identical across section
function [theta, v] = transverse_velocity_profile_rozovskii(v,h,r,It,C,T,dtheta);
	kappa = 0.41;
	g     = 9.81;


	% initial distribution
	% ODE to determine initial velocity distribution
	% dv^2/dr = d v0^2/dr - 2 v^2/r (preceeding 2.32)
	% for v0 == condt: law of areas, rozovskii, 2.30
	vt = sqrt(h./r);

	% scale to discharge
	vt = Q0./sum(centre(v).*diff(r))*vt;

	% differential operator along the radius (assume constant step with)
	nr = length(r);
	dr = (r(2)-r(1));
	Dr = spdiags(0.5/dr*ones(nr,1)*[-1 0 1],-1:1,nr,nr);

	% neumann boundary condition
	if (1)
		Dr(1,:) = 0;
		Dr(1,1:2) = [-1,1]/dr;
		Dr(end,:) = 0;
		Dr(end,end-1:end) = [-1,1]/dr;
	else
		% dirichlet bc
		Dr(1,:) = 0;
		Dr(1,1) = [1]/dr;
		Dr(end,:) = 0;
		Dr(end,end) = [-1]/dr;
	end
	
	% steps along theta
	theta = cvec(T(1):dtheta:T(2));

	% pseudo variable
	u = r.*v.^2.*h.^2;

	% initial momentum (discharge)
	Q0 = sum(v.*h);

	% sove PDE
	for idx=2:length(theta)
		du_dtheta_i(:,idx) = du_dtheta(theta,u(:,idx-1));

		% step euler forward along theta
		u(:,idx) = u(:,idx-1)+dtheta*du_dtheta_i(:,idx);
		v(:,idx) = sqrt(u(:,idx)./r)./h;
		v = real(v);

		% apply boundary conditions
		% u(1,idx)   = 0;
		% u(end,idx) = 0;

%		v(1,idx) = 0;
%		v(end,idx) = 0;
%		v(1,idx)   = 2*v(2,idx) - v(3,idx);
%		v(end,idx) = 2*v(end-1,idx)-v(end-2,idx);

		% apply conservation of momentum
		% TODO this is a hack, one should instead implement a momentum presevering scheme
		% or itheta should be iterated
		s        = Q0/sum(v(:,idx).*h);
		v(:,idx) = s*v(:,idx);
		u(:,idx) = r.*v(:,idx).^2.*h.^2;
	end % for idx

	function du_dtheta = du_dtheta(theta,u)
		% eq. 2.37' in rozovskii
		du_dtheta = - 0.75/(kappa^3*C)*sqrt(g)*h.*(Dr*u) ...
			    - g/C^2*r.*u./h ...
			    + g*It.*r.*h.^2;
	end % du_dtheta

end % transverse_velocity_profile_rosovskii

