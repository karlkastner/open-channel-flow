% 2011-02-28 01:35:34.000000000 +0100
%% determine total energy as sump of potential and kinetic energy
%% this is preserved for fricitionless flows
% SWE::energy
function [E, Epot, Ekin] = energy(H,h0,V)
	g = Constant.gravity;
	% ok for tidal wave, but the same for water flowing down a slope?
	Epot = sum(1/2*g*(H-h0).^2);
	Ekin = sum(1/2*H.*V.^2); % or hbar here?
	E = Epot+Ekin;
end

