% Fri 30 Oct 20:48:28 +08 2020
% Fri 30 Oct 19:11:03 +08 2020
function [resc,resm,rescf,resmf] = residual_swe(t,x,A,Q,g,w,zb,Cd,nf)
	omega = 2*pi/86400;
	nx = length(x);
	nt = length(t);
%	Dx = obj.D1_dx(1,x);
	Dx = derivative_matrix_1_1d(x,2);
	Dt = derivative_matrix_1_1d(t,2);
	Dx = kron(Dx,speye(nt));
	Dt = kron(speye(nx),Dt);

	% continuity
	resc = (Dt*A + Dx*Q);
	 
	% momentum
	%resm = (      Dt*Q + Dx*(Q.^2./A) + 0.5*g*Dx*(A.^2) ...
	resm = (      Dt*Q + Dx*(Q.^2./A) + g*A.*(Dx*A) ...
                  - ( -g*A.*(Dx*zb) + g*(A./w).^2.*(Dx*w) - Cd.*w.*abs(Q).*Q./A.^2 ) ...
	       );

	% normalize by width
	resc = resc./w;
	resm = resm./w;

	resc = reshape(resc,nx,nt);
	resm = reshape(resm,nx,nt);

	% TODO fourier decomposition of the residual
	o = zeros(1,1+2*nf);
	o(1) =1;
	for idx=1:nf
		o(2*idx)   = omega*idx;
		o(2*idx+1) = omega*idx;
	end
	if (nf>0)
		F = fourier_matrix_exp(1./(1:nf),t);
		rescf = (F'*resc.');
		resmf = (F'*resm.');
		% normalize
		rescf = rescf ./o';
		resmf = resmf ./o';
	end
end % residual swe

