% Tue 18 Feb 14:48:13 +08 2020
% ponce 1991
% dQ/dt + (dQ/dA)dQ/dx = nu dQ^2/dx
function nu = flood_wave_diffusion_coefficient(Q,w,h,S)
	U  = Q./(w.*h);
	nu = Q./(2*w*S).*(1-U.^2);
end

