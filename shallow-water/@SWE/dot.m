% Thu Apr 28 04:50:17 MSD 2011
% Mi 3. Feb 12:26:04 CET 2016
% Karl Kästner, KTH Stockholm
%
%% time derivative
%% (only for matlab internal ode-solver)
%% TODO this is not swe specific
%
%% continuity
%% dA/dt + dQ/dx = I
%% 
%% momentum
%% dQ/dt + d/dx( Qu + 1/2 gh^2) = gA(S_f - S_b)
%% S_b = dz_b/dx
%% S_f = tau_x/rho_w = C_f u|u|
%
% TODO : eddy viscosity
% TODO : velocity profile correction (vertical and transversal)
%
% SWE::dot
function u_dot = dot(t, u, flux, bc)
        % apply boundary condition
        u = feval(bc, t, u);
        f = feval(flux, t, u);
        u_dot = -0.5/dx*Dc*f;
end % u_dot

