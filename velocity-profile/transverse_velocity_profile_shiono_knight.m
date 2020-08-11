% Sun 23 Jul 11:47:45 CEST 2017
%% transverse profile of the streamwise velocity, determined analytically
%% by the method of shiono and knight
%% shape of velocity profile only dependent on lambda, f, H, not slope
% y : across channel coordinate
% S : slope
% H : depth
% w : width
function [U, dU_dy, y0, app] = transverse_velocity_profile_shiono_knight(w,H,S,f,lambda,y)
	if (nargin()<5 || isempty(lambda))
		lambda = 0.5;
	end
	if (nargin()<6 || isempty(y))
		y = linspace(-w/2,w/2);
	end


	g = Constant.gravity;
	
	% free stream velocity
	U0 = sign(S).*sqrt(8*g*abs(S)*H/f)
	
	% constant depth (eq 10)
	% for derivation of A1 and A2 see separate script 
	gamma = sqrt(2/lambda*sqrt(f/8))*1./H;
	
	A1 = -1/(exp(gamma*w/2) + exp(-gamma*w/2))

	% velocity
	U  = U0.*sqrt(1 + A1*(exp(gamma*y) + exp(-gamma*y)));
	
	% higher order leads to 1% difference, so not necessary
	% y0 = log((exp(-(gamma*w)/2)*(exp(gamma*w) + (2*exp(gamma*w) + exp(2*gamma*w) - gamma^2*w^2*exp(gamma*w) + 1)^(1/2) + 1))/(gamma*w))/gamma
	y0 = -(w/2 + log(2/(gamma*w))/gamma);
	
	if (nargout() > 1)
		% TODO this has a singularity at 0
		dU_dy = -0.5*U0 ...
			   .* A1*gamma*(exp(gamma*y) - exp(-gamma*y)) ...
			   ./  sqrt(1 + A1*(exp(gamma*y) + exp(-gamma*y)));
	end
	
	%Ubara = U0*(gamma*w - exp(gamma*w) + gamma*w*exp(gamma*w) + 1)/(gamma*w*(exp(gamma*w) + 1));
	
	if (nargout() > 3)
		app.U    = sqrt(8*g*S*H/f)*(1 + 1/2*A1*(exp(gamma*y) + exp(-gamma*y)));
		app.dU   = sqrt(8*g*S*H/f)*(1/2*gamma*(exp(gamma*(-w/2-y))));
		app.Ubar = U0*(1 - 1/(w*gamma));
	end
	
	if (0)
	% du/dy = 0 at y = 0
	% u(+/-w/2) = 0
	
	% with side slop (eq 11)
	a1 = -1/2 + 1/2*sqrt(1 + s*sqrt(1+s^2)/lambda*sqrt(8*f));
	a1 = +1/2 + 1/2*sqrt(1 + s*sqrt(1+s^2)/lambda*sqrt(8*f));
	omega = g*S/(sqrt(1+s^2)/s*f/8 - lambda/s^2*sqrt(f/8));
	
	U = sqrt(A3*Y^a1 + A4*Y^-a2 + omega*Y);
	
	end % if 0
	
end

